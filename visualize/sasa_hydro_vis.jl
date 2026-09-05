# Visual check for the hydration-shell dummy placement `SASA.sasa` derives
# from its Shrake-Rupley patch geometry (centroid + patch_rg2 per exposed
# atom). Independent of `visualize/sasa_vis.jl` -- that file is about the
# raw occlusion machinery (`SASA._occluded`/`SASA._classify`); this one is
# about what a hydration-shell dummy generator would actually place on top
# of it, so it defines its own scene and drawing helpers rather than reusing
# that file's.
#
# Every scene uses `HydroRadii`, a test-only `RadiiSource` defined below, so
# the geometry is exact and reproducible instead of depending on the
# atomic-radii database.
#
# Run with:
#   julia --project=visualize -e 'include("visualize/sasa_hydro_vis.jl"); vis_sasa_hydro()'
#
# Numbers only, no window (works headless):
#   julia --project=visualize -e 'include("visualize/sasa_hydro_vis.jl"); sasa_hydro_report()'

using ScatterNet
using ScatterNet.Interfaces: Interfaces, RadiiSource
using ScatterNet.Molecule.Molecules: Molecules, Molecule
using ScatterNet.Molecule.SASA: SASA
using ScatterNet.Molecule.SASA.PlasticMap: plastic_points
using Printf: @printf, @sprintf
using GLMakie

# --------------------------------------------------------------------------
# radii source
# --------------------------------------------------------------------------

"""
Test-only [`RadiiSource`](@ref) mapping element labels to hand-picked radii,
so the packed-sphere scene has exact, reproducible geometry instead of
whatever the atomic-radii database happens to hold.
"""
struct HydroRadii <: RadiiSource
    table::Dict{String,Float64}
end

"""
    lookup(src::HydroRadii, ions) -> Vector{Tuple{String,Union{Float64,Nothing}}}

Resolve each label against `src.table`, `nothing` for anything absent.
"""
function Interfaces.lookup(src::HydroRadii, ions)
    out = Vector{Tuple{String,Union{Float64,Nothing}}}(undef, length(ions))
    for (i, ion) in enumerate(ions)
        s = String(ion)
        out[i] = (s, get(src.table, s, nothing))
    end
    return out
end

# --------------------------------------------------------------------------
# the scene: a real packed cluster (FCC lattice carved to a sphere)
# --------------------------------------------------------------------------

"`0:n-1` grid over 3 axes, as a flat vector of `(i,j,k)`."
_fcc_grid(n::Int) = vec([(i, j, k) for i in 0:n-1, j in 0:n-1, k in 0:n-1])

"""
    fcc_lattice(n, a) -> Vector{NTuple{3,Float64}}

Face-centred-cubic lattice sites over an `n×n×n` block of conventional cells
of side `a` (4 atoms/cell). Nearest-neighbour distance is `a/√2`.
"""
function fcc_lattice(n::Int, a::Float64)
    g = _fcc_grid(n)
    pts = [(a * i, a * j, a * k) for (i, j, k) in g]
    for (dx, dy, dz) in ((0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))
        append!(pts, [(a * (i + dx), a * (j + dy), a * (k + dz)) for (i, j, k) in g])
    end
    return pts
end

"""
    packed_cluster_scene(; probe = 1.4, r = 1.5, n = 4, a = 2.2, frac = 0.42) -> NamedTuple

A real packed object: `n³` FCC conventional cells of equal spheres (`r = 1.5`
by default), trimmed to the sites within `a*n*frac` of the centroid 
At the defaults this is 80 atoms, ~2/3 exposed and ~1/3 buried.
"""
function packed_cluster_scene(; probe::Float64 = 1.4, r::Float64 = 1.5, n::Int = 4, a::Float64 = 2.2, frac::Float64 = 0.42)
    pts = fcc_lattice(n, a)
    ctr = ntuple(t -> sum(p[t] for p in pts) / length(pts), 3)
    keep = [p for p in pts if sqrt(sum((p[t] - ctr[t])^2 for t in 1:3)) <= a * n * frac]
    src = HydroRadii(Dict("A" => r))
    mol = Molecules.create("packed cluster", fill("A", length(keep)), keep;
                            radii_source = src)
    return (; mol, probe, title = "Packed cluster ($(length(keep)) spheres, FCC)")
end

# --------------------------------------------------------------------------
# pure geometry: plastic sample points, no plotting
# --------------------------------------------------------------------------

const EXPOSED_COLOR  = RGBf(1.00, 0.62, 0.13)   # warm/bright
const OCCLUDED_COLOR = RGBf(0.24, 0.26, 0.32)   # dark/desaturated
const DUMMY_COLOR    = RGBf(0.90, 0.15, 0.55)   # distinct from atoms and points
const ATOM_COLOR     = RGBf(0.55, 0.62, 0.75)

"""
    atom_sample_points(mol, i, n, probe) -> NamedTuple

The `n` plastic-sequence directions mapped onto atom `i`'s expanded sphere
(radius `radii(mol)[i] + probe`), each tested with `SASA._occluded` against
every other atom.

# Returns
`(; pts, exposed)`: the `n` sample points and which of them are solvent
accessible.
"""
function atom_sample_points(mol::Molecule, i::Int, n::Int, probe::Float64)
    crds = Molecules.coords_cartesian(mol)
    rads = Molecules.radii(mol)
    cands = collect(1:size(crds, 2))

    ρ = rads[i] + probe
    cx = crds[1, i]; cy = crds[2, i]; cz = crds[3, i]
    dirs = plastic_points(n)

    pts = Vector{NTuple{3,Float64}}(undef, n)
    exposed = Vector{Bool}(undef, n)
    @inbounds for j in 1:n
        ux, uy, uz = dirs[j]
        p = (cx + ρ * ux, cy + ρ * uy, cz + ρ * uz)
        pts[j] = p
        exposed[j] = !SASA._occluded(p, cands, crds, rads, probe, i)
    end
    return (; pts, exposed)
end

# --------------------------------------------------------------------------
# headless report
# --------------------------------------------------------------------------

"""
    sasa_hydro_report(; probe = 1.4, n_exp = 4096) -> Nothing

Print the patch geometry `SASA.sasa` derives for hydration-shell dummy
placement on [`packed_cluster_scene`](@ref): how many atoms are buried vs
exposed, total SASA, and the per-exposed-atom area/centroid/patch-radius
table.
"""
function sasa_hydro_report(; probe::Float64 = 1.4, n_exp::Int = 4096)
    sc = packed_cluster_scene(; probe)
    area, centroid, patch_rg2, exposed = SASA.sasa(sc.mol; n_exp, probe)
    n_atoms = length(area)
    println(sc.title, "  (probe = ", probe, ", n_exp = ", n_exp, ")")
    @printf("  %d atoms: %d exposed (dummy sites), %d buried\n",
            n_atoms, count(exposed), count(!, exposed))
    @printf("  total SASA %.3f A^2\n\n", sum(area))
    for i in axes(centroid, 2)
        exposed[i] || continue
        @printf(
            "  atom %3d: area %8.3f  centroid (% .3f, % .3f, % .3f)  patch radius %.3f A\n",
            i, area[i], centroid[1, i], centroid[2, i], centroid[3, i],
            sqrt(max(patch_rg2[i], 0.0))
        )
    end
    return nothing
end

# --------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------

"""
    sasa_hydro_figure(; n_exp = 4096, n_pts = 400, probe = 1.4) -> Figure

Hydration-shell dummy placement on a real packed object
([`packed_cluster_scene`](@ref), ~80 atoms): every exposed atom gets a dummy
marker (its patch centroid from `SASA.sasa`, sized by its SASA area) so the
result reads as an actual shell coating the cluster's outside, with its
buried core visibly left uncoated.

Split out from [`vis_sasa_hydro`](@ref) so the plot can be assembled, saved
or inspected without a window.
"""
function sasa_hydro_figure(; n_exp::Int = 4096, n_pts::Int = 400, probe::Float64 = 1.4)
    sc = packed_cluster_scene(; probe)
    crds = Molecules.coords_cartesian(sc.mol)
    rads = Molecules.radii(sc.mol)
    natoms = size(crds, 2)
    area, centroid, patch_rg2, exposed = SASA.sasa(sc.mol; n_exp, probe)

    exp_idx = findall(exposed)
    bur_idx = findall(!, exposed)
    i_exposed = exp_idx[argmax(area[exp_idx])]                        # most-exposed atom
    highlight = isempty(bur_idx) ? (i_exposed,) : (i_exposed, first(bur_idx))

    fig = Figure(size = (1000, 900))
    ax = Axis3(
        fig[1, 1];
        title = @sprintf(
            "%s: hydration-shell dummy placement\n%d exposed (dummy sites), %d buried  (probe = %.2f A)",
            sc.title, length(exp_idx), length(bur_idx), probe),
        titlesize = 15, aspect = :data, azimuth = 1.1π,
        xlabel = "x", ylabel = "y", zlabel = "z")

    # atom spheres: small, uniform, low-alpha -- reads as the cluster's shape
    # rather than competing with the dummies for attention. The two
    # highlighted atoms get a stronger alpha so their sample points land on
    # something visible.
    for i in 1:natoms
        ρ = Float32(rads[i] + probe)
        c = Point3f(crds[1, i], crds[2, i], crds[3, i])
        mesh!(ax, Sphere(c, ρ); color = (ATOM_COLOR, i in highlight ? 0.35 : 0.10),
              transparency = true, shading = NoShading)
    end

    # raw sample points, only on the two highlighted atoms
    for i in highlight
        st = atom_sample_points(sc.mol, i, n_pts, probe)
        ex = [Point3f(p...) for (p, e) in zip(st.pts, st.exposed) if e]
        oc = [Point3f(p...) for (p, e) in zip(st.pts, st.exposed) if !e]
        isempty(oc) || scatter!(ax, oc; color = OCCLUDED_COLOR, markersize = 5)
        isempty(ex) || scatter!(ax, ex; color = EXPOSED_COLOR, markersize = 5)
    end

    # hydration-shell dummies: one per exposed atom, sized by its SASA area
    # so a bigger patch reads as a bigger dummy instead of every marker
    # looking identical.
    amax = maximum(area[exp_idx])
    dummy_pts = [Point3f(centroid[1, i], centroid[2, i], centroid[3, i]) for i in exp_idx]
    dummy_sizes = [6.0 + 14.0 * sqrt(area[i] / amax) for i in exp_idx]
    scatter!(ax, dummy_pts; color = DUMMY_COLOR, markersize = dummy_sizes,
             marker = :circle, strokewidth = 0.5, strokecolor = (:black, 0.3))

    els = [ MarkerElement(color = ATOM_COLOR, marker = :circle, markersize = 14),
            MarkerElement(color = DUMMY_COLOR, marker = :circle, markersize = 14),
            MarkerElement(color = EXPOSED_COLOR, marker = :circle, markersize = 10),
            MarkerElement(color = OCCLUDED_COLOR, marker = :circle, markersize = 10)]
    Legend(fig[2, 1], els,
           [   "atom (expanded radius)", "hydration-shell dummy (size ~ patch area)",
               "sample point: exposed (highlighted atoms only)",
               "sample point: occluded (highlighted atoms only)"];
           orientation = :horizontal, framevisible = false, nbanks = 2, labelsize = 12)
    return fig
end

"""
    vis_sasa_hydro(; n_exp = 4096, n_pts = 400, probe = 1.4) -> Nothing

Display [`sasa_hydro_figure`](@ref): every exposed atom's hydration-shell
dummy on a real ~80-atom packed cluster, plus the raw sample points on one
buried and one exposed atom to show the underlying mechanism. Blocks until
the window is closed.

# Keywords
- `n_exp`: sample points per atom for the `SASA.sasa` computation itself
  (all atoms).
- `n_pts`: sample points drawn for the two highlighted atoms only.
- `probe`: solvent probe radius.
"""
vis_sasa_hydro(; n_exp::Int = 4096, n_pts::Int = 400, probe::Float64 = 1.4) =
    wait(display(sasa_hydro_figure(; n_exp, n_pts, probe)))
