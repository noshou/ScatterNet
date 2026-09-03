"""
A molecule: per-atom coordinates centered at the centroid, with `r`, `theta`,
`phi` computed eagerly and `radii`, `vols` lazily.
"""
module  Molecules

using   ..Cache: Lazy, force
using   ..Interfaces: RadiiSource, lookup
using   ..AtomicRadii: AtomicRadiiSource

export  Molecule, MoleculeError, create, r, theta, phi, 
        coords, radii, vols, elms, name

struct MoleculeError <: Exception; msg::String end
Base.showerror(io::IO, e::MoleculeError) = print(io, "MoleculeError: ", e.msg)

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

sphere_volume(rad::Float64)::Float64 = (4.0 / 3.0) * pi * rad^3

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
    create(name, elms, coords; radii_source = AtomicRadiiSource()) -> Molecule

`coords` is per-atom `(x, y, z)` (any frame); it is centered at the centroid.
`r`, `theta`, `phi` are computed now; `radii`, `vols` on first access.
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

r(m::Molecule)::Vector{Float64}      = m._r
theta(m::Molecule)::Vector{Float64}  = m._theta
phi(m::Molecule)::Vector{Float64}    = m._phi
coords(m::Molecule)::Matrix{Float64} = m._coords
radii(m::Molecule)::Vector{Float64}  = force(m._radii)
vols(m::Molecule)::Vector{Float64}   = force(m._vols)
elms(m::Molecule)::Vector{String}    = m._elms
name(m::Molecule)::String            = m._name

end # module
