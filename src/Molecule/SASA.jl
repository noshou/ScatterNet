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
Å² of accessible surface each shell point stands for; sets the cloud's spacing
at `≈ √SHELL_AREA_PER_POINT ≈ 2 Å`.

Calibrated against CRYSOL: its default `--fb 17` puts `F(17) = 1597` points on a
typical globular protein (~6500 Å²), i.e. ~4 Å² each. Budgeting by *area* rather
than by a fixed count keeps that spacing at every size, and since surface area
grows as `N^(2/3)`, the point count is sub-linear in atom count rather than flat.
CRYSOL instead caps `--fb` at `F(18) = 2584` for any structure, which
under-resolves large complexes; pass `n_target` explicitly to reproduce that.
"""
const SHELL_AREA_PER_POINT = 4.0

"Floor on the derived point budget: `F(10)`, the smallest Fibonacci grid CRYSOL's `--fb` accepts."
const SHELL_MIN_POINTS = 55

"Sample directions per atom, before occlusion and before thinning. Internal: only fine enough to resolve one atom's patch."
const _SHELL_SAMPLE = 256

"Range (Å) over which [`_bead_class`](@ref) casts escape rays. A void whose wall is further than this in every direction is bulk solvent, not a cavity."
const _BEAD_RAY_RANGE = 12.0

"Directions sampled by [`_bead_class`](@ref); about half fall in the outward hemisphere and are used."
const _BEAD_RAY_DIRS = 64

"Escaping fraction at or above which a bead is [`CONVEX`](@ref); below it (but nonzero) [`CONCAVE`](@ref)."
const _BEAD_CONVEX_ESCAPE = 0.5

"""
    BeadClass

Where a hydration-shell bead sits, matching CRYSOL 3's three border-layer
populations (each carries its own fitted contrast; CRYSOL's defaults are
`1, 1, 0` in units of `0.03 e/Å³`).

- `CONVEX`:  outer surface, open solvent ahead of it.
- `CONCAVE`: outer surface but recessed -- a groove or pocket.
- `CAVITY`:  enclosed interior void, unreachable from outside.
"""
@enum BeadClass CONVEX CONCAVE CAVITY

"""
    _ray_blocked(p, d, nb, crds, rads, probe) -> Bool

Does the ray from `p` along unit `d` hit any expanded sphere in `nb`?

Standard ray/sphere test: project each centre onto the ray, reject anything
behind `p`, and compare the perpendicular offset against `r + probe`. The
bead's own atom always projects backwards, so it never self-blocks.
"""
function _ray_blocked(
    p::NTuple{3,Float64}, d::Vec3, nb::Vector{Int},
    crds::Matrix{Float64}, rads::Vector{Float64}, probe::Float64
)::Bool
    @inbounds for j in nb
        wx = crds[1, j] - p[1]; wy = crds[2, j] - p[2]; wz = crds[3, j] - p[3]
        t = wx * d[1] + wy * d[2] + wz * d[3]
        t <= 0.0 && continue
        ρ = rads[j] + probe
        (wx * wx + wy * wy + wz * wz) - t * t < ρ * ρ && return true
    end
    return false
end

"""
    _bead_class(p, n̂, nb, dirs, crds, rads, probe) -> BeadClass

Classify one shell bead by what fraction of its outward hemisphere escapes the
molecule. Rays are cast only over [`_BEAD_RAY_RANGE`](@ref), and the outward normal is
tried first so open surface -- the common case -- costs one ray.
"""
function _bead_class(
    p::NTuple{3,Float64}, n̂::NTuple{3,Float64}, nb::Vector{Int},
    dirs::Vector{Vec3}, crds::Matrix{Float64}, rads::Vector{Float64},
    probe::Float64
)::BeadClass
    isempty(nb) && return CONVEX
    _ray_blocked(p, n̂, nb, crds, rads, probe) || return CONVEX

    esc = 0; tot = 0
    @inbounds for d in dirs
        d[1] * n̂[1] + d[2] * n̂[2] + d[3] * n̂[3] > 0.0 || continue
        tot += 1
        _ray_blocked(p, d, nb, crds, rads, probe) || (esc += 1)
    end
    tot == 0 && return CONVEX

    esc == 0 && return CAVITY
    return esc / tot >= _BEAD_CONVEX_ESCAPE ? CONVEX : CONCAVE
end

"""
    _prefix_thin(counts, m, budget) -> Vector{Int}

Indices selecting a prefix of each atom's block, sized in proportion to that
block, totalling exactly `budget` out of `m` points.

Allocation is cumulative-floor (Bresenham): atom `i` gets
`floor(C_i·budget/m) - floor(C_{i-1}·budget/m)` where `C_i` is the running
accepted count. The differences sum to `budget` with no rounding drift, and each
is within one point of the proportional share.
"""
function _prefix_thin(counts::Vector{Int}, m::Int, budget::Int)::Vector{Int}
    idx = Vector{Int}(undef, 0); sizehint!(idx, budget)
    base = 0; prev = 0
    @inbounds for c in counts
        c == 0 && continue
        cum = base + c
        take = (cum * budget) ÷ m - prev
        for t in 1:take
            push!(idx, base + t)          # prefix of this atom's block
        end
        base = cum; prev = (cum * budget) ÷ m
    end
    return idx
end

"""
    _class_loop(tree, pts, nrm, crds, rads, probe, dirs) -> Vector{BeadClass}

Classification pass behind [`shell_points`](@ref). A second function barrier for
the same reason [`_shell_loop`](@ref) is one: the `inrange` query per bead would
otherwise dispatch dynamically on the non-concrete tree type, once per bead.
"""
function _class_loop(
    tree::T, pts::Matrix{Float64}, nrm::Matrix{Float64}, crds::Matrix{Float64},
    rads::Vector{Float64}, probe::Float64, dirs::Vector{Vec3}
)::Vector{BeadClass} where {T}
    out = Vector{BeadClass}(undef, size(pts, 2))
    @inbounds for k in axes(pts, 2)
        nb = inrange(tree, view(pts, :, k), _BEAD_RAY_RANGE)
        out[k] = _bead_class(
            (pts[1, k], pts[2, k], pts[3, k]), (nrm[1, k], nrm[2, k], nrm[3, k]),
            nb, dirs, crds, rads, probe)
    end
    return out
end

"""
    shell_points(mol; probe, n_target) -> (pts, areas, class)

The solvent-accessible surface of `mol` as a point cloud: the area each point
stands for, and which CRYSOL border-layer population it belongs to.

Where [`sasa`](@ref) reduces each patch to scalars, this keeps the patch itself
-- what a hydration-shell model needs, since a patch's spatial extent is the
thing that scatters.

Each atom is sampled at [`_SHELL_SAMPLE`](@ref) directions and the cloud is then
thinned to `n_target` points for the whole molecule, by keeping a *prefix* of
each atom's accepted points sized in proportion to that atom's count -- hence to
its area. Every survivor then carries an equal `sum(area)/M`, so `sum(areas)`
still matches `sum(sasa(mol)[1])`.

Prefixes, not strides. The plastic sequence is progressive, so any prefix is
itself low-discrepancy, while an evenly strided subset is a *different* Kronecker
sequence whose quality depends on how near-rational the stride makes `k/ρ`.
Measured on a 256-point set thinned to 64: prefix gives first-moment `|mean|`
`0.04` and worst cap discrepancy `0.078`, striding gives `0.48` and `0.44` --
the strided points bunch on one side of the sphere, since `frac(4/ρ) ≈ 0.0195`
advances azimuth only ~7° per kept point.

# Arguments
- `mol`: molecule whose surface to sample.

# Keywords
- `probe`: solvent probe radius in Å; `probe >= 0`. Default `1.4` (water).
- `n_target`: total points to keep; `> 0`. Default `nothing`, deriving it from
the accessible area via [`SHELL_AREA_PER_POINT`](@ref) so spacing stays fixed
as the molecule grows. A cloud already smaller than the budget is kept whole.

# Returns
-   `pts::Matrix{Float64}`, `(3, M)`: accessible points in `mol`'s centred
    cartesian frame, sharing the origin `coords_cartesian` uses.
-   `areas::Vector{Float64}`, `(M,)`: Å² per point, equal across all `M`.
-   `class::Vector{BeadClass}`, `(M,)`: per-point [`BeadClass`](@ref). Cavity
    detection is exact for voids up to [`_BEAD_RAY_RANGE`](@ref) across;
    anything larger classifies as open surface, which is the intent -- a void
    that wide holds bulk-like water, not ordered shell water.

`M` is `0` for a molecule with no accessible surface; `pts` is then `(3, 0)`.
"""
function shell_points(
    mol::Molecule;
    probe::Float64                  = 1.4,
    n_target::Union{Nothing,Int}    = nothing
)::Tuple{Matrix{Float64},Vector{Float64},Vector{BeadClass}}
    n_target === nothing || n_target > 0 ||
        throw(DomainError(n_target, "n_target must be > 0"))
    probe >= 0.0 || throw(DomainError(probe, "probe must be >= 0"))

    pmap = PlasticMap.plastic_points(_SHELL_SAMPLE)
    rads = Molecules.radii(mol)
    rmax = Molecules.r_max(mol)
    crds = Molecules.coords_cartesian(mol)
    tree = KDTree(crds)
    pts, areas, nrm, counts =
        _shell_loop(tree, crds, rads, rmax, pmap, probe, _SHELL_SAMPLE)

    total = sum(areas)
    m = length(areas)
    budget = n_target === nothing ?
        max(SHELL_MIN_POINTS, round(Int, total / SHELL_AREA_PER_POINT)) : n_target

    if m > budget
        idx = _prefix_thin(counts, m, budget)
        pts, nrm = pts[:, idx], nrm[:, idx]
        areas = fill(total / length(idx), length(idx))
    end

    class = _class_loop(tree, pts, nrm, crds, rads, probe,
                        PlasticMap.plastic_points(_BEAD_RAY_DIRS))
    return pts, areas, class
end

"""
    _shell_loop(tree, crds, rads, rmax, pmap, probe, n_pts) -> (pts, areas, nrm)

Per-atom loop behind [`shell_points`](@ref), split out for the same reason
[`_sasa_loop!`](@ref) is: `KDTree` over a bare `Matrix` has no concrete type at
the call site, so the barrier lets Julia specialize.

Reuses [`_classify`](@ref)'s shortcuts. Unlike `sasa` there is no
witness/confirm split -- the points are the output, so nothing can be settled
early. `nrm` carries each point's outward unit normal, which only this loop
knows and [`_bead_class`](@ref) needs; `counts` gives each atom's block length,
which [`_prefix_thin`](@ref) needs.
"""
function _shell_loop(
    tree::T,
    crds::Matrix{Float64},
    rads::Vector{Float64},
    rmax::Float64,
    pmap::Vector{Vec3},
    probe::Float64,
    n_pts::Int
)::Tuple{Matrix{Float64},Vector{Float64},Matrix{Float64},Vector{Int}} where {T}

    xs = Float64[]; ys = Float64[]; zs = Float64[]; areas = Float64[]
    nx = Float64[]; ny = Float64[]; nz = Float64[]
    counts = zeros(Int, size(crds, 2))

    for i in axes(crds, 2)
        x = crds[1, i]; y = crds[2, i]; z = crds[3, i]
        ρ = rads[i] + probe
        per_pt = 4 * π * ρ^2 / n_pts   # every direction stands for this much

        candidates = inrange(tree, @view(crds[:, i]), ρ + rmax + probe)
        status = _classify(i, candidates, crds, rads, probe)
        status == ALL_BURIED && continue
        keep_all = status == ALL_EXPOSED

        got = 0
        @inbounds for j in 1:n_pts
            ux, uy, uz = pmap[j]
            p = (x + ρ * ux, y + ρ * uy, z + ρ * uz)
            (keep_all || !_occluded(p, candidates, crds, rads, probe, i)) || continue
            push!(xs, p[1]); push!(ys, p[2]); push!(zs, p[3])
            push!(nx, ux); push!(ny, uy); push!(nz, uz)
            push!(areas, per_pt)
            got += 1
        end
        counts[i] = got
    end

    pts = Matrix{Float64}(undef, 3, length(areas))
    nrm = Matrix{Float64}(undef, 3, length(areas))
    @inbounds for k in eachindex(areas)
        pts[1, k] = xs[k]; pts[2, k] = ys[k]; pts[3, k] = zs[k]
        nrm[1, k] = nx[k]; nrm[2, k] = ny[k]; nrm[3, k] = nz[k]
    end
    return pts, areas, nrm, counts
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
    sasa(mol; probe, n_occ, n_exp, area_tol) -> (area, exposed)

Per-atom solvent-accessible surface area of `mol`, via Shrake–Rupley point
sampling over a plastic-sequence point set.

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
-   `exposed`: `(n,)`, `true` where the atom has at least one accessible point.
"""
function sasa(
    mol::Molecule;
    probe::Float64    = 1.4,
    n_occ::Int        = 512,
    n_exp::Int        = 4096,
    area_tol::Float64 = 2.0
)::Tuple{Vector{Float64},Vector{Bool}}

    _check_sasa_args(probe, n_occ, n_exp, area_tol)

    pmap = PlasticMap.plastic_points(n_exp)
    rads = Molecules.radii(mol)
    rmax = Molecules.r_max(mol)
    crds = Molecules.coords_cartesian(mol)

    tree = KDTree(crds)
    n = size(crds, 2)

    areas   = zeros(Float64, n)
    exposed = falses(n)

    return _sasa_loop!(
        areas, exposed, tree, crds, rads, rmax, pmap, probe, n_occ, n_exp, area_tol
    )
end

"""
    _sasa_loop!(areas, exposed, tree, crds, rads, rmax, pmap, probe,
                n_occ, n_exp, area_tol) -> (areas, exposed)

Per-atom loop behind [`sasa`](@ref). `KDTree(crds)` can't infer a concrete tree 
type from a bare `Matrix` (NearestNeighbors keys the tree type on point dimension, 
a runtime property of the array), so calling out to a separate function lets Julia
specialize the whole loop on it once instead of dispatching `inrange` per atom.
"""
function _sasa_loop!(
    areas::Vector{Float64},
    exposed::BitVector,
    tree::T,
    crds::Matrix{Float64},
    rads::Vector{Float64},
    rmax::Float64,
    pmap::Vector{Vec3},
    probe::Float64,
    n_occ::Int,
    n_exp::Int,
    area_tol::Float64
)::Tuple{Vector{Float64},Vector{Bool}} where {T}

    for i in axes(crds, 2)
        x = crds[1, i]; y = crds[2, i]; z = crds[3, i]
        ρ = rads[i] + probe
        full = 4 * π * ρ^2

        candidates = inrange(tree, @view(crds[:, i]), ρ + rmax + probe)
        status = _classify(i, candidates, crds, rads, probe)

        if status == ALL_BURIED
            continue                                # areas/exposed stay 0/false
        elseif status == ALL_EXPOSED
            areas[i] = full
            exposed[i] = true
            continue
        end

        # AMBIGUOUS: only point sampling can settle the fraction. Count the
        # first n_occ points, and bail early if even the worst case they leave
        # open is negligible.
        cnt = 0
        @inbounds for j in 1:n_occ
            ux, uy, uz = pmap[j]
            p = (x + ρ * ux, y + ρ * uy, z + ρ * uz)
            _occluded(p, candidates, crds, rads, probe, i) || (cnt += 1)
        end

        # No witness among n_occ, which can still admit a true fraction up to 
        # ~3/n_occ, so only give up on it when even that much area is negligible.
        cnt == 0 && 3.0 / n_occ * full < area_tol && continue

        @inbounds for j in n_occ+1:n_exp
            ux, uy, uz = pmap[j]
            p = (x + ρ * ux, y + ρ * uy, z + ρ * uz)
            _occluded(p, candidates, crds, rads, probe, i) || (cnt += 1)
        end

        areas[i] = full * (cnt / n_exp)
        cnt > 0 && (exposed[i] = true)
    end
    return areas, exposed
end

end # module SASA
