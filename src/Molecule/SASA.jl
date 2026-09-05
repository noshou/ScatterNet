"""
Solvent-accessible surface area estimation.
"""
module SASA

"""
Even point set on the unit sphere, drawn from the 
2-D plastic (R₂) low-discrepancy sequence.
"""
module PlasticMap

using Roots: find_zero

export Vec3, plastic_points, PLASTIC_RATIO

"A unit vector on the sphere, `(x, y, z)`."
const Vec3 = NTuple{3,Float64}

"Cubic `x³ − x − 1`; its real root in `(1, 2)` is the plastic ratio."
_plastic_poly(x) = x^3 - x - 1

"Plastic ratio ρ ≈ 1.324718, the real root of `x³ = x + 1`."
const PLASTIC_RATIO     = find_zero(_plastic_poly, (1.0, 2.0))
const PLASTIC_RATIO_SQR = PLASTIC_RATIO^2

"Fractional part of `x`, i.e. `x - floor(x)`, in `[0, 1)`."
_frac(x) = x - floor(x)

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
    _check_sasa_args(probe, n_occ, n_exp, area_tol) -> Nothing

Shared argument contract for [`sasa`](@ref); throws `DomainError` on the
first violated constraint.
"""
function _check_sasa_args(probe::Float64, n_occ::Int, n_exp::Int, area_tol::Float64)::Nothing
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
    return nothing
end

"""
    _unit_patch_moments(pmap, n_exp) -> (Vec3, Float64)

Centroid and mean-squared spread of `n_exp`-point plastic set on
the unit sphere: `(ubar, rg2_unit)` where `ubar = mean(pmap[1:n_exp])` and
`rg2_unit = mean(|u|²) - |ubar|² = 1 - |ubar|²` (every `u` is a unit vector, 
so `mean(|u|²) = 1`).
"""
function _unit_patch_moments(pmap::Vector{Vec3}, n_exp::Int)::Tuple{Vec3,Float64}
    sx = 0.0; sy = 0.0; sz = 0.0
    @inbounds for j in 1:n_exp
        ux, uy, uz = pmap[j]
        sx += ux; sy += uy; sz += uz
    end
    nf = Float64(n_exp)
    ubar = (sx / nf, sy / nf, sz / nf)
    rg2_unit = 1.0 - (ubar[1]^2 + ubar[2]^2 + ubar[3]^2)
    return ubar, rg2_unit
end

"""
    sasa(mol; probe, n_occ, n_exp, area_tol) -> (area, centroid, patch_rg2, exposed)

Per-atom solvent-accessible surface area of `mol`, via Shrake–Rupley point
sampling over a plastic-sequence point set, plus the accessible surface's
local geometry: the area-weighted centroid of each atom's accessible sample
points, and their mean squared spread about that centroid (`patch_rg2`). 

If there are no neighbours reaching the surface, area is exactly `4π(r+probe)²`; 
if a single neighbour engulfs it, area is exactly `0`; either way this costs 
`O(k)` (`k` = # of neighbours) with no sampling and, for the free-exposed case,
Only the ambiguous case must, tested with `n_occ` points first, and if there exists
a single non-occluded point (a witness) or the worst-case remaining exposure exceeds
`area_tol`, it is confirmed against the full `n_exp` set. Non-existence of a
witness is *not* proof of burial, so by the rule of three up to `3/n_occ` of
the sphere could still be exposed.

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
                Default `2.0` Å².

# Returns
-   `area`: `(n,)`, indexed like `coords_cartesian(mol)`'s columns.
-   `centroid`: `(3, n)`, indexed the same way. Equal to the atom's own
    position where `exposed` is `false` (unused there).
-   `patch_rg2`: `(n,)`, mean squared distance of accessible points from
    `centroid`, Å². Zero where `exposed` is `false`.
-   `exposed`: `(n,)`, `true` where the atom has at least one accessible point.
"""
function sasa(
    mol::Molecule;
    probe::Float64    = 1.4,
    n_occ::Int        = 512,
    n_exp::Int        = 4096,
    area_tol::Float64 = 2.0
)::Tuple{Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Bool}}

    _check_sasa_args(probe, n_occ, n_exp, area_tol)

    pmap = PlasticMap.plastic_points(n_exp)
    ubar, rg2_unit = _unit_patch_moments(pmap, n_exp)

    rads = Molecules.radii(mol)
    rmax = Molecules.r_max(mol)
    crds = Molecules.coords_cartesian(mol)

    tree = KDTree(crds)
    n = size(crds, 2)

    areas     = zeros(Float64, n)
    centroid  = Matrix{Float64}(undef, 3, n)
    patch_rg2 = zeros(Float64, n)
    exposed   = falses(n)

    return _sasa_loop!(
        areas, centroid, patch_rg2, exposed, tree, crds,
        rads, rmax, pmap, ubar, rg2_unit, probe, n_occ,
        n_exp, area_tol
    )
end

"""
    _sasa_loop!(areas, centroid, patch_rg2, exposed, tree, crds, rads,
                rmax, pmap, ubar, rg2_unit, probe, n_occ, n_exp,
                area_tol) -> (areas, centroid, patch_rg2, exposed)

Per-atom loop behind [`sasa`](@ref). `KDTree(crds)` can't infer a concrete tree 
type from a bare `Matrix` (NearestNeighbors keys the tree type on point dimension, 
a runtime property of the array), so calling out to a separate function lets Julia 
specialise on the concrete type at the call, statically dispatching `inrange` inside 
the per-atom loop instead of once per atom. 
"""
function _sasa_loop!(
    areas::Vector{Float64},
    centroid::Matrix{Float64},
    patch_rg2::Vector{Float64},
    exposed::BitVector,
    tree::T,
    crds::Matrix{Float64},
    rads::Vector{Float64},
    rmax::Float64,
    pmap::Vector{Vec3},
    ubar::Vec3,
    rg2_unit::Float64,
    probe::Float64,
    n_occ::Int,
    n_exp::Int,
    area_tol::Float64
)::Tuple{Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Bool}} where {T}

    for i in axes(crds, 2)
        x = crds[1, i]; y = crds[2, i]; z = crds[3, i]
        ρ = rads[i] + probe
        full = 4 * π * ρ^2

        candidates = inrange(tree, @view(crds[:, i]), ρ + rmax + probe)
        status = _classify(i, candidates, crds, rads, probe)

        if status == ALL_BURIED
            centroid[1, i] = x; centroid[2, i] = y; centroid[3, i] = z
            continue                                # areas/patch_rg2/exposed stay 0/0/false
        elseif status == ALL_EXPOSED
            areas[i] = full
            exposed[i] = true
            centroid[1, i] = x + ρ * ubar[1]
            centroid[2, i] = y + ρ * ubar[2]
            centroid[3, i] = z + ρ * ubar[3]
            patch_rg2[i] = ρ^2 * rg2_unit
            continue
        end

        # AMBIGUOUS: no shortcut exists for a centroid, so every accepted point
        # must be folded into a running sum. Accumulate over the first n_occ
        # points, in local (atom-centred) coordinates `v = ρ*u`; the parallel-
        # axis theorem (E[|v|²] - |E[v]|²) then gives patch_rg2 from the same
        # single pass that gives the centroid, no second loop over the points.
        sx = 0.0; sy = 0.0; sz = 0.0; s2 = 0.0; cnt = 0
        @inbounds for j in 1:n_occ
            ux, uy, uz = pmap[j]
            p = (x + ρ * ux, y + ρ * uy, z + ρ * uz)
            if !_occluded(p, candidates, crds, rads, probe, i)
                vx = ρ * ux; vy = ρ * uy; vz = ρ * uz
                sx += vx; sy += vy; sz += vz
                s2 += vx * vx + vy * vy + vz * vz
                cnt += 1
            end
        end

        # No witness among n_occ, which can still admit a true fraction up to 
        # ~3/n_occ, so only give up on it when even that much area is negligible.
        if cnt == 0 && 3.0 / n_occ * full < area_tol
            centroid[1, i] = x; centroid[2, i] = y; centroid[3, i] = z
            continue
        end

        @inbounds for j in n_occ+1:n_exp
            ux, uy, uz = pmap[j]
            p = (x + ρ * ux, y + ρ * uy, z + ρ * uz)
            if !_occluded(p, candidates, crds, rads, probe, i)
                vx = ρ * ux; vy = ρ * uy; vz = ρ * uz
                sx += vx; sy += vy; sz += vz
                s2 += vx * vx + vy * vy + vz * vz
                cnt += 1
            end
        end

        areas[i] = full * (cnt / n_exp)
        if cnt > 0
            exposed[i] = true
            cx = sx / cnt; cy = sy / cnt; cz = sz / cnt
            centroid[1, i] = x + cx; centroid[2, i] = y + cy; centroid[3, i] = z + cz
            patch_rg2[i] = s2 / cnt - (cx^2 + cy^2 + cz^2)
        else
            centroid[1, i] = x; centroid[2, i] = y; centroid[3, i] = z
        end
    end
    return areas, centroid, patch_rg2, exposed
end

end # module SASA
