"""
A molecule: per-atom coordinates centered at the centroid, with `r`, `theta`,
`phi` computed eagerly and `radii`, `vols` lazily.
"""
module  Molecules

using   ..Cache: Lazy, force
using   ...Interfaces: RadiiSource, lookup
using   ..AtomicRadii: AtomicRadiiSource

export  Molecule, MoleculeError, create, r, theta, phi, 
        coords, radii, vols, elms, name

"Raised for malformed molecule input (empty or mismatched coords, missing radii)."
struct MoleculeError <: Exception; msg::String end
Base.showerror(io::IO, e::MoleculeError) = print(io, "MoleculeError: ", e.msg)

"Per-atom centered coords with eager `r`/`theta`/`phi` and lazy `radii`/`vols`."
struct Molecule
    _name::String
    _elms::Vector{String}
    _coords::Matrix{Float64}          # (3, n) centered
    _r::Vector{Float64}               # (n,)
    _theta::Vector{Float64}           # (n,)
    _phi::Vector{Float64}             # (n,)
    _radii::Lazy{Vector{Float64}}
    _vols::Lazy{Vector{Float64}}
end

"Volume of a sphere of radius `rad`."
sphere_volume(rad::Float64)::Float64 = (4.0 / 3.0) * pi * rad^3

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
    _geometry(c::Matrix{Float64}) -> (Vector{Float64}, Vector{Float64}, Vector{Float64})

Spherical `(r, theta, phi)` per column of `c`, with `theta = acos(z/r)` in
`[0, π]` and `phi = atan(y, x)`. `r = 0` is handled without a `0/0`.

# Arguments
- `c`: `(3, n)` centered coordinate matrix.
"""
function _geometry(c::Matrix{Float64})
    n = size(c, 2)
    r = Vector{Float64}(undef, n); th = Vector{Float64}(undef, n); ph = Vector{Float64}(undef, n)
    @inbounds for j in 1:n
        x = c[1, j]; y = c[2, j]; z = c[3, j]
        rj = sqrt(x * x + y * y + z * z)
        r[j] = rj
        rsafe = rj > 0.0 ? rj : 1.0
        th[j] = acos(clamp(z / rsafe, -1.0, 1.0))
        ph[j] = atan(y, x)
    end
    return r, th, ph
end

"""
    _compute_radii(src::RadiiSource, es::Vector{String}) -> Vector{Float64}

Resolve per-element radii through `src`; throws `MoleculeError` on an empty list
or any element with no radius data.

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
        out[i] = rad
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

Build a `Molecule`: `coords` are centered at the centroid, `r`/`theta`/`phi`
computed now, `radii`/`vols` on first access.

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
    c = _center(cs)
    rr, th, ph = _geometry(c)
    rad = Lazy{Vector{Float64}}(() -> _compute_radii(radii_source, es))
    vol = Lazy{Vector{Float64}}(() -> sphere_volume.(force(rad)))
    return Molecule(String(name), es, c, rr, th, ph, rad, vol)
end

"Per-atom radial distance from the centroid."
r(m::Molecule)::Vector{Float64}      = m._r
"Per-atom polar angle `acos(z/r)` in `[0, π]`."
theta(m::Molecule)::Vector{Float64}  = m._theta
"Per-atom azimuth `atan(y, x)`."
phi(m::Molecule)::Vector{Float64}    = m._phi
"`(3, n)` centroid-centered coordinate matrix."
coords(m::Molecule)::Matrix{Float64} = m._coords
"Per-atom radius; resolved and cached on first call."
radii(m::Molecule)::Vector{Float64}  = force(m._radii)
"Per-atom sphere volume; computed and cached on first call."
vols(m::Molecule)::Vector{Float64}   = force(m._vols)
"Element/ion string per atom."
elms(m::Molecule)::Vector{String}    = m._elms
"Molecule label."
name(m::Molecule)::String            = m._name

end # module
