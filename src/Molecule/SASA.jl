"""
Solvent-accessible surface area estimation.
"""
module SASA

"""
Even point set on the unit sphere, drawn from the 
2-D plastic (R₂) low-discrepancy sequence.
"""
module PlasticMap

using ....Constants: PLASTIC_RATIO, PLASTIC_RATIO_SQR, _frac

export Vec3, plastic_points, PLASTIC_RATIO

"A unit vector on the sphere, `(x, y, z)`."
const Vec3 = NTuple{3,Float64}

"""
    _plastic_point(i::Int) -> Vec3

Unit-sphere point for 1-based plastic-sequence term `i`. The 2-D term
`(frac(i/ρ), frac(i/ρ²))` is read as (azimuth, height) and lifted to the sphere
through the equal-area cylindrical projection, so the points are uniform in area
rather than clustered at the poles.

Division by `ρ`/`ρ²` (rather than multiplication by reciprocals) keeps the
fractional part accurate; it degrades only once `i` nears the mantissa limit
(~1e15), far above any SASA point count.
"""
@inline function _plastic_point(i::Int)::Vec3
    φ = 2.0 * π * _frac(i / PLASTIC_RATIO)
    z = 2.0 * _frac(i / PLASTIC_RATIO_SQR) - 1.0
    r = sqrt(max(0.0, 1.0 - z * z))
    sinφ, cosφ = sincos(φ)
    return (r * cosφ, r * sinφ, z)
end

"""
    plastic_points(n::Int) -> Vector{Vec3}

The first `n` plastic-sequence points distributed uniformly on a 3D unit sphere 
using Lambert's cylindrical  equal-area projection. The longitudinal angle φ maps 
horizontally, while the vertical component `z` scales uniformly between [-1.0, 1.0], 
representing the sine of the latitude.

# Arguments
- `n`: number of points to generate; `n >= 0`.
"""
function plastic_points(n::Int)::Vector{Vec3}
    n < 0 && throw(DomainError(n, "n must be >= 0"))
    pts = Vector{Vec3}(undef, n)
    @inbounds for i in 1:n
        pts[i] = _plastic_point(i)
    end
    return pts
end

end # module PlasticMap

using .PlasticMap: PlasticMap, Vec3
using NearestNeighbors: KDTree, inrange
using ..Molecules: Molecules, Molecule

"""
    _occluded(p, candidates, crds, rads, probe, self) -> Bool

Whether point `p` lies inside the expanded sphere
(`radius + probe`) of any candidate atom other than `self`.

`self` is skipped because `p` is generated *on* atom `self`'s own expanded
sphere: its distance to that center is exactly `rads[self] + probe`, so an
unguarded `dst <= ρ_c` would report every point of every atom as occluded,
and `sasa` would return `0.0` for every molecule.

# Arguments
- `p`: point to test.
- `candidates`: indices of atoms to test against.
- `crds`: `(3, n)` coordinate matrix.
- `rads`: per-atom radius, indexed like `crds`'s columns.
- `probe`: solvent probe radius.
- `self`: index of the atom `p` was sampled on; never occludes `p`.
"""
function _occluded(
    p::Vec3,
    candidates::Vector{Int},
    crds::Matrix{Float64},
    rads::Vector{Float64},
    probe::Float64,
    self::Int
)::Bool
    x_p, y_p, z_p = p
    @inbounds for c in candidates
        c == self && continue
        ρ_c = rads[c] + probe
        x_c = crds[1, c]; y_c = crds[2, c]; z_c = crds[3, c]
        dst² = (x_c - x_p)^2 + (y_c - y_p)^2 + (z_c - z_p)^2
        if dst² <= ρ_c * ρ_c
            return true
        end
    end
    return false
end

"""
    Coverage

How much of an atom's expanded sphere its neighbours cover.
- `ALL_EXPOSED`:    no neighbour reaches the sphere, so the exposed fraction is
                    exactly 1 and the area is exactly `4π(r+probe)²`.
- `ALL_BURIED`:     a single neighbour swallows the whole sphere, so the exposed
                    fraction is exactly 0.
- `AMBIGUOUS`:      neighbours cut caps but no single one settles it; only point
                    sampling can estimate the fraction.
"""
@enum Coverage ALL_EXPOSED ALL_BURIED AMBIGUOUS

"""
    _classify(i, candidates, crds, rads, probe) -> Coverage

Each neighbour `j` cuts a spherical cap out of `i`'s expanded sphere. Writing
`ρᵢ = rads[i] + probe`, `ρⱼ = rads[j] + probe` and `d = |cᵢ - cⱼ|`, three cases
are decidable by comparing scalars, with no point sampling at all:

- `d + ρᵢ <= ρⱼ`:   `j` engulfs `i` entirely, so *every* point is occluded.
- `d >= ρᵢ + ρⱼ`:   `j` is too far to reach `i`'s surface, so it cuts nothing.
- `d + ρⱼ <= ρᵢ`:   `j`'s ball sits strictly inside `i`'s surface, so it also
                    cuts nothing (`i` encloses `j`).

If no neighbour cuts a cap the atom is fully exposed. Everything else is a
union-of-caps question that this predicate deliberately does not answer.

Sampling can only prove an atom is not fully covered, it can never prove burial, 
since nothing is stoping the next point from beign covered.

# Arguments
- `i`: atom to classify.
- `candidates`: neighbour indices from the coarse range query; may include `i`.
- `crds`: `(3, n)` coordinate matrix.
- `rads`: per-atom radius, indexed like `crds`'s columns.
- `probe`: solvent probe radius.
"""
function _classify(
    i::Int,
    candidates::Vector{Int},
    crds::Matrix{Float64},
    rads::Vector{Float64},
    probe::Float64
)::Coverage
    ρ_i = rads[i] + probe
    x_i = crds[1, i]; y_i = crds[2, i]; z_i = crds[3, i]
    cuts = false
    @inbounds for j in candidates
        j == i && continue
        ρ_j = rads[j] + probe
        d² = (crds[1, j] - x_i)^2 + (crds[2, j] - y_i)^2 + (crds[3, j] - z_i)^2

        # a neighbour that swallows i settles it outright
        ρ_j >= ρ_i && d² <= (ρ_j - ρ_i)^2 && return ALL_BURIED

        # neighbours that never reach i's surface cut nothing
        (d² >= (ρ_i + ρ_j)^2 || (ρ_i >= ρ_j && d² <= (ρ_i - ρ_j)^2)) && continue
        cuts = true
    end
    return cuts ? AMBIGUOUS : ALL_EXPOSED
end

"""
    sasa_atoms(mol; probe, n_occ, n_exp, area_tol) -> Vector{Float64}

Per-atom solvent-accessible surface area of `mol`, via Shrake–Rupley point
sampling over a plastic-sequence point set. Every atom first goes through 
[`_classify`](@ref). If there are no neighbours  reaching the surface, it will return 
`4π(r+probe)²`, and if there exists a single  engulfing neighbour it will return `0` 
(runs in `O(k)`, where k is # neighbours). Ambigous cases are tested with `n_occ` points, 
and if there exists a single non-occluded  point (a witness) we calculate its SASA 
contribution as `4π(r+probe)² × (exposed / n_exp)`. Non-existence of a witness is *not* 
proof of burial, so by the rule of 3 up to 3/n_occ points could be exposed. Code will keep
tracking unless 3/n is less than the tolerance.  

# Arguments
- `mol`: molecule to score.

# Keywords
- `probe`:      solvent probe radius; `probe >= 0`. Default 1.4 Å is the radius of a 
                water molecule and is the convention for SASA calculations.
- `n_occ`:      points for the witness pass; `n_occ > 0`. 512 is the smallest round count
                measured to lose no area on a dense lattice. Catches any atom exposed by more 
                than ~1/512 of its sphere, about 0.2 Å²
- `n_exp`:      Points per atom for the exposed-fraction pass. Default measured relative error 
                against the analytic two-sphere cap: 1.3 % at 64 points, 0.36 % at 1024, **0.065 %
                at 4096**, 0.02 % at 16384. Costs ~0.09 ms/atom.
- `area_tol`:   if no exposed point is found in `n_occ` samples, the atom might still
                have a tiny exposed patch (≤ 3/n_occ of its sphere). If that worst‑case area is
                below `area_tol`, we skip the full `n_exp` pass and treat it as buried.
                Default `0.8` Å² (safe for most uses); set `0.0` for exact (slower).
# Returns
Area per atom, indexed like `coords_cartesian(mol)`'s columns.
"""
function sasa_atoms(
    mol::Molecule;                        
    probe::Float64    = 1.4,    
    n_occ::Int        = 512,              
    n_exp::Int        = 4096,             
    area_tol::Float64 = 0.8
)::Vector{Float64}

    # Assertion checks
    if (n_occ <= 0)
        throw(DomainError(n_occ, "n_occ must be > 0"))
    elseif (n_exp <= 0)
        throw(DomainError(n_exp, "n_exp must be > 0"))
    elseif (n_exp < n_occ)
        throw(DomainError((n_occ, n_exp), "n_occ must be <= n_exp"))
    elseif (probe < 0.0)
        throw(DomainError(probe, "probe must be >= 0"))
    elseif (area_tol < 0.0)
        throw(DomainError(area_tol, "area_tol must be >= 0"))
    end

    # Calculate set of mappings for plastic points
    pmap = PlasticMap.plastic_points(n_exp)

    # Radii, coords, and maximum radius
    rads = Molecules.radii(mol)
    rmax = Molecules.r_max(mol)
    crds = Molecules.coords_cartesian(mol)

    tree = KDTree(crds)

    # SASA_i = 4π(rᵢ + probe)² × (n_exposed / n_exp)
    areas = zeros(Float64, size(crds, 2))

    # `tree` is not concretely typed, `KDTree(crds)` cannot infer from a Matrix, because NearestNeighbors 
    # keys the tree type on the point dimension, which is a runtime property of the array. Calling out to 
    # a separate function lets Julia specialise on the concrete tree  type at the call, so the `inrange` 
    # insidethe per-atom loop is statically dispatched instead of once per atom.
    return _sasa_loop!( areas, tree, crds, rads, rmax, pmap, probe, n_occ, n_exp, area_tol)
end

"""
    _sasa_loop!(areas, tree, crds, rads, rmax, pmap, probe, n_occ, n_exp, area_tol) -> areas

Per-atom scoring loop behind [`sasa_atoms`](@ref)'s function barrier; `tree` is
specialised to its concrete type here. Fills and returns `areas`.
"""
function _sasa_loop!(
    areas::Vector{Float64},
    tree::T,
    crds::Matrix{Float64},
    rads::Vector{Float64},
    rmax::Float64,
    pmap::Vector{Vec3},
    probe::Float64,
    n_occ::Int,
    n_exp::Int,
    area_tol::Float64
)::Vector{Float64} where {T}

    pts = Vector{Vec3}(undef, n_exp)

    for i in axes(crds, 2)

        # expanded radius of atom_i plus probe
        ρ = rads[i] + probe
        full = 4 * π * ρ^2

        # Do coarse filter to filter out atoms which are too far away to ever
        # occlude atom_i. Since  we don't know every radius at first, by treating
        # the closest  possible relevant neighbor as if it could be as large as the
        # single biggest atom in the molecule. Therefore, after this coarse filter we
        # have an upper bound of all atoms which could be *potentially* occluding atom_i
        # (a superset of occulding atoms).
        candidates = inrange(tree, @view(crds[:,i]), ρ + rmax + probe)

        # Settle the two exactly-decidable regimes straight from the neighbour
        # list, before spending a single point on them.
        status = _classify(i, candidates, crds, rads, probe)
        if status == ALL_BURIED
            continue                    # areas[i] stays 0.0, exactly
        elseif status == ALL_EXPOSED
            areas[i] = full             # exactly 4π(r+probe)²
            continue
        end

        # surface_pt = (cᵢ[1] + ρᵢ * p[1], cᵢ[2] + ρᵢ * p[2], cᵢ[3] + ρᵢ * p[3])
        # map minimum num of points onto surface of i
        x = crds[1, i]; y = crds[2, i]; z = crds[3, i]
        @inbounds for j in 1:n_occ
            ux, uy, uz = pmap[j]
            pts[j] = Vec3((x + ρ*ux, y + ρ*uy, z + ρ*uz))
        end

        # Hunt for a single unoccluded witness; one is enough to prove exposure,
        # so the coarse pass stops at the first rather than counting them.
        first_exp = 0
        @inbounds for j in 1:n_occ
            if !_occluded(pts[j], candidates, crds, rads, probe, i)
                first_exp = j
                break
            end
        end

        # No witness, rule of three: 0 exposed of n_occ still admits a true fraction up to 
        # ~3/n_occ, so only write the atom off when even that much area is negligible. Otherwise
        # fall through and confirm against the full n_exp set.
        if first_exp == 0 && 3.0 / n_occ * full < area_tol
            continue
        end

        # extend to the higher-res mesh, reusing the prefix already built
        @inbounds for j in n_occ+1:n_exp
            ux, uy, uz = pmap[j]
            pts[j] = Vec3((x + ρ*ux, y + ρ*uy, z + ρ*uz))
        end

        # Carry the coarse pass's verdicts forward instead of re-testing them.
        # Either it stopped at a witness or it ran to n_occ with every point occluded.
        p_exp, from = first_exp == 0 ? (0, n_occ) : (1, first_exp)
        @inbounds for j in from+1:n_exp
            if !_occluded(pts[j], candidates, crds, rads, probe, i)
                p_exp += 1
            end
        end
        areas[i] = full * (p_exp / n_exp)
    end
    return areas
end

"""
    sasa(mol; probe, n_occ, n_exp, area_tol) -> Float64

Total solvent-accessible surface area of `mol`: the sum of
[`sasa_atoms`](@ref). See that function for the sampling scheme and the
argument contract.
"""
sasa(mol::Molecule; kwargs...)::Float64 = sum(sasa_atoms(mol; kwargs...))

end # module SASA
