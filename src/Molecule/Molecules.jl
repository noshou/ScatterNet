"""
A molecule: per-atom coordinates centered at the centroid, held in both
cartesian and spherical form (computed eagerly), with `radii`, `vols` and
`r_max` lazily.

Lazy memoization (`Lazy`/`make`/`force`) and the default radii backend
(`Ion`/`resolve_one`/`AtomicRadiiSource`/...) are spliced directly into this
module from `Cache.jl`/`AtomicRadii.jl` rather than living in their own
submodules -- nothing outside `Molecules` uses either, so there's no
separate module identity worth keeping (Julia has no true private scoping;
not exporting these names is as close as it gets).
"""
module  Molecules

import  ...Interfaces
using   ...Interfaces: RadiiSource, lookup
using   SQLite: SQLite
using   DBInterface: DBInterface

include("Cache.jl")
include("AtomicRadii.jl")

export  Molecule, MoleculeError, create, coords_cartesian, coords_spherical,
        radii, vols, r_max, elms, name

"Raised for malformed molecule input (empty or mismatched coords, missing radii)."
struct MoleculeError <: Exception; msg::String end
Base.showerror(io::IO, e::MoleculeError) = print(io, "MoleculeError: ", e.msg)

"""
Per-atom centered coordinates in both frames, with lazy `radii`/`vols`/`r_max`.

Both coordinate frames are `(3, n)` matrices sharing a column index (the atom),
so a single atom's data is one contiguous column in either frame; the spherical
rows are `r`, `theta`, `phi` in that order.
"""
struct Molecule
    _name   :: String
    _elms   :: Vector{String}
    _cart   :: Matrix{Float64}       # (3, n) centered (x, y, z)
    _sph    :: Matrix{Float64}       # (3, n) (r, theta, phi)
    _radii  :: Lazy{Vector{Float64}}
    _vols   :: Lazy{Vector{Float64}}
    _r_max  :: Lazy{Float64}         # largest per-atom radius
end

"Volume of a sphere of radius `rad`."
sphere_volume(rad::Float64)::Float64 = (4.0 / 3.0) * π * rad^3

"""
    _center(cs::Vector{NTuple{3,Float64}}) -> Matrix{Float64}

Stack coordinates into a `(3, n)` matrix translated to the centroid.

# Arguments
- `cs`: per-atom `(x, y, z)` tuples; must be non-empty.
"""
function _center(cs::Vector{NTuple{3,Float64}})::Matrix{Float64}
    n = length(cs)
    n == 0 && throw(MoleculeError("Empty coordinates"))
    sx = 0.0; sy = 0.0; sz = 0.0
    @inbounds for c in cs
        sx += c[1]; sy += c[2]; sz += c[3]
    end
    nf = Float64(n)
    mx = sx / nf; my = sy / nf; mz = sz / nf
    out = Matrix{Float64}(undef, 3, n)
    @inbounds for j in 1:n
        c = cs[j]
        out[1, j] = c[1] - mx; out[2, j] = c[2] - my; out[3, j] = c[3] - mz
    end
    return out
end

# theta = acos(z/r) in [0, π], phi = atan2(y, x). r = 0 (single atom) would give
# 0/0, so clamp with rsafe; j_l(0) = 0 for l > 0 makes the angle irrelevant there.
"""
    _geometry(c::Matrix{Float64}) -> Matrix{Float64}

Spherical `(r, theta, phi)` per column of `c` as a `(3, n)` matrix, with
`theta = acos(z/r)` in `[0, π]` and `phi = atan(y, x)`. `r = 0` is handled
without a `0/0`.

# Arguments
- `c`: `(3, n)` centered coordinate matrix.
"""
function _geometry(c::Matrix{Float64})::Matrix{Float64}
    n = size(c, 2)
    out = Matrix{Float64}(undef, 3, n)
    @inbounds for j in 1:n
        x = c[1, j]; y = c[2, j]; z = c[3, j]
        rj = sqrt(x * x + y * y + z * z)
        rsafe = rj > 0.0 ? rj : 1.0
        out[1, j] = rj
        out[2, j] = acos(clamp(z / rsafe, -1.0, 1.0))
        out[3, j] = atan(y, x)
    end
    return out
end

"""
    _compute_radii(src::RadiiSource, es::Vector{String}) -> Vector{Float64}

Resolve per-element radii through `src`; throws `MoleculeError` on an empty list
or any element with no radius data.

A negative radius is clamped to `0.0`. Shannon's tables carry a handful of
these (`h1+`, `c4+`, `n5+`) as extrapolation artifacts of fitting to
coordination-number trends, not as physical sizes. Clamping keeps the
downstream invariants that actually matter -- non-negative volumes, and an
expanded SASA radius `r + probe` that never inverts -- and a bare proton with
no electron density around it is, for scattering purposes, exactly the
zero-radius object the clamp makes it. These ions are rare enough in practice
that anything relying on them is not stable input for scattering analysis
regardless.

# Arguments
- `src`: radii backend to query.
- `es`: element/ion strings, one per atom.
"""
function _compute_radii(src::S, es::Vector{String})::Vector{Float64} where {S<:RadiiSource}
    isempty(es) && throw(MoleculeError("Empty elements"))
    pairs = lookup(src, es)
    out = Vector{Float64}(undef, length(pairs))
    @inbounds for i in eachindex(pairs)
        el, rad = pairs[i]
        rad === nothing && throw(MoleculeError("no radius data for element \"$el\""))
        out[i] = max(0.0, rad)   # negative table entries are artifacts; see above
    end
    return out
end

"""
    _to_tuples(cs) -> Vector{NTuple{3,Float64}}

Normalize any iterable of 3-component coordinates to `Float64` tuples (identity
when already `Vector{NTuple{3,Float64}}`). Throws `MoleculeError` if an entry
lacks 3 components.

# Arguments
- `cs`: iterable of per-atom coordinates.
"""
_to_tuples(cs::Vector{NTuple{3,Float64}}) = cs
function _to_tuples(cs)
    out = Vector{NTuple{3,Float64}}(undef, length(cs))
    @inbounds for (i, c) in enumerate(cs)
        length(c) == 3 || throw(MoleculeError("each coordinate needs 3 components"))
        out[i] = (Float64(c[1]), Float64(c[2]), Float64(c[3]))
    end
    return out
end

"""
    create(name, elms, coords; radii_source::RadiiSource = AtomicRadiiSource()) -> Molecule

Build a `Molecule`: `coords` are centered at the centroid, both coordinate
frames computed now, `radii`/`vols`/`r_max` on first access.

# Arguments
- `name`: molecule label.
- `elms`: element/ion string per atom.
- `coords`: per-atom `(x, y, z)` in any frame; length must match `elms`.

# Keywords
- `radii_source`: radii backend (defaults to `AtomicRadiiSource()`).
"""
function create(name::AbstractString, elms::AbstractVector{<:AbstractString}, coords;
                radii_source::RadiiSource = AtomicRadiiSource())
    cs = _to_tuples(coords)
    length(cs) == length(elms) || throw(MoleculeError("coords and elms length mismatch"))
    es = collect(String, elms)
    cart = _center(cs)
    sph  = _geometry(cart)
    rad  = Lazy{Vector{Float64}}(() -> _compute_radii(radii_source, es))
    vol  = Lazy{Vector{Float64}}(() -> sphere_volume.(force(rad)))
    rmax = Lazy{Float64}(() -> maximum(force(rad)))
    return Molecule(String(name), es, cart, sph, rad, vol, rmax)
end

"`(3, n)` centroid-centered cartesian coordinates; rows are `x`, `y`, `z`."
coords_cartesian(m::Molecule)::Matrix{Float64} = m._cart

"`(3, n)` spherical coordinates about the centroid; rows are `r`, `theta`, `phi`."
coords_spherical(m::Molecule)::Matrix{Float64} = m._sph

"Per-atom radius; resolved and cached on first call."
radii(m::Molecule)::Vector{Float64}  = force(m._radii)

"Per-atom sphere volume; computed and cached on first call."
vols(m::Molecule)::Vector{Float64}   = force(m._vols)

"""
Largest per-atom radius in the molecule; forces (and caches) `radii`.

SASA's coarse neighbour filter needs this to bound how far away an atom can
still occlude another, before any individual radius is known.
"""
r_max(m::Molecule)::Float64          = force(m._r_max)

"Element/ion string per atom."
elms(m::Molecule)::Vector{String}    = m._elms

"Molecule label."
name(m::Molecule)::String            = m._name

end # module
