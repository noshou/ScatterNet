"""
A molecule.
"""
module  Molecules

import  ...Interfaces
using   ...Interfaces: RadiiSource, lookup
using   ...Interfaces: AtomicRadiiSource

include("Cache.jl")

export  Molecule, MoleculeError, create, coords_cartesian, coords_spherical,
        to_spherical, radii, vols, r_max, elms, name, n_atoms

"Raised for malformed molecule input (empty or mismatched coords, missing radii)."
struct MoleculeError <: Exception; msg::String end
Base.showerror(io::IO, e::MoleculeError) = print(io, "MoleculeError: ", e.msg)

"""
Per-atom centered coordinates in both frames, with lazy `radii`/`vols`/`r_max`.

Both coordinate frames are `(3, n)` matrices sharing a column index (the atom),
so a single atom's data is one contiguous column in either frame; the spherical
rows are `r`, `theta`, `phi` in that order. `_n` is the atom count, captured once
at construction from the coordinate pass rather than recomputed on demand.
"""
struct Molecule
    _name   :: String
    _elms   :: Vector{String}
    _n      :: Int                   # atom count; set at construction, never recomputed
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

"""
    to_spherical(c::AbstractMatrix{<:Real}) -> Matrix{Float64}

Spherical `(r, theta, phi)` per column of the `(3, n)` cartesian matrix `c`,
returned as a `(3, n)` matrix sharing `c`'s column index (the atom/point).
`theta = acos(z/r)` lies in `[0, π]` and `phi = atan(y, x)` in `(-π, π]`.

`r = 0` would make `theta` a `0/0`; it is handled without one. The angle is arbitrary there
and unobservable downstream, since `j_l(0) = 0` for every `l > 0`.

# Arguments
- `c`: `(3, n)` cartesian coordinates; rows are `x`, `y`, `z`.
"""
function to_spherical(c::AbstractMatrix{<:Real})::Matrix{Float64}
    size(c, 1) == 3 || throw(MoleculeError(
        "to_spherical: expected a (3, n) matrix with rows (x, y, z); got $(size(c, 1)) rows"))
    n = size(c, 2)
    out = Matrix{Float64}(undef, 3, n)
    @inbounds for j in 1:n
        x = Float64(c[1, j]); y = Float64(c[2, j]); z = Float64(c[3, j])
        rj = sqrt(x * x + y * y + z * z)
        rsafe = rj > 0.0 ? rj : 1.0   # see the r = 0 note above
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
coordination-number trends, not as physical sizes.

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
- `elms`: element/ion string per atom, e.g. `"c"`, `"Fe"`, `"o2-"`. Case is
        normalized to lowercase (the radii and form-factor tables are lowercase-keyed),
        so `elms(m)` returns the lowercased strings.
- `coords`: per-atom `(x, y, z)` in any frame; length must match `elms`.

# Keywords
- `radii_source`: radii backend (defaults to `AtomicRadiiSource()`).
"""
function create(name::AbstractString, elms::AbstractVector{<:AbstractString}, coords;
                radii_source::RadiiSource = AtomicRadiiSource())
    cs = _to_tuples(coords)
    n  = length(cs)
    n == length(elms) || throw(MoleculeError("coords and elms length mismatch"))
    es = String[lowercase(e) for e in elms]   # radii/form-factor tables are lowercase-keyed
    cart = _center(cs)
    sph  = to_spherical(cart)
    rad  = Lazy{Vector{Float64}}(() -> _compute_radii(radii_source, es))
    vol  = Lazy{Vector{Float64}}(() -> sphere_volume.(force(rad)))
    rmax = Lazy{Float64}(() -> maximum(force(rad)))
    return Molecule(String(name), es, n, cart, sph, rad, vol, rmax)
end

"`(3, n)` centroid-centered cartesian coordinates; rows are `x`, `y`, `z`."
coords_cartesian(m::Molecule)::Matrix{Float64} = m._cart

"`(3, n)` spherical coordinates about the centroid; rows are `r`, `theta`, `phi`."
coords_spherical(m::Molecule)::Matrix{Float64} = m._sph

"Number of atoms; `O(1)`, captured at construction."
n_atoms(m::Molecule)::Int = m._n

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
