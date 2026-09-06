# Visual check for the hydration-shell dummy cloud `SASA.shell_points`
# builds; every solvent-accessible sample point becomes one dummy, so the
# shell is a resolved layer over the molecular surface rather than one marker
# per exposed atom. 
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
const CONCAVE_COLOR  = RGBf(0.20, 0.55, 0.90)
const CAVITY_COLOR   = RGBf(0.35, 0.80, 0.35)
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
    sasa_hydro_report(; probe = 1.4, n_target = nothing) -> Nothing

Print the hydration-shell dummy cloud on [`packed_cluster_scene`](@ref): buried
vs exposed atoms, and how the cloud behaves as the global budget varies. The
area column should stay flat Budgets shown bracket CRYSOL's `--fb` range (`F(10) = 55` to
`F(18) = 2584`, default `F(17) = 1597`).
"""
function sasa_hydro_report(; probe::Float64 = 1.4, n_target::Union{Nothing,Int} = nothing)
    sc = packed_cluster_scene(; probe)
    area, exposed = SASA.sasa(sc.mol; probe)
    println(sc.title, "  (probe = ", probe, ")")
    @printf("  %d atoms: %d exposed, %d buried\n",
            length(area), count(exposed), count(!, exposed))
    @printf("  total SASA %.3f A^2\n\n", sum(area))

    @printf("  %-10s %-9s %-14s %-8s %-8s %-8s\n",
            "n_target", "dummies", "shell area A^2", "convex", "concave", "cavity")
    for n in (55, 233, 610, 1597, 2584)
        pts, pa, cl = SASA.shell_points(sc.mol; probe, n_target = n)
        @printf("  %-10d %-9d %-14.3f %-8d %-8d %-8d\n", n, size(pts, 2), sum(pa),
                count(==(SASA.CONVEX), cl), count(==(SASA.CONCAVE), cl),
                count(==(SASA.CAVITY), cl))
    end

    pts, pa, _ = SASA.shell_points(sc.mol; probe, n_target)
    @printf("\n  auto budget: %d dummies, %.4f A^2 each\n",
            size(pts, 2), isempty(pa) ? 0.0 : first(pa))
    return nothing
end

# --------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------

"""
    sasa_hydro_figure(; n_target = nothing, n_show = 400, probe = 1.4) -> Figure

The hydration-shell dummy cloud on a real packed object
([`packed_cluster_scene`](@ref), ~80 atoms): every accessible sample point from
`SASA.shell_points` is drawn as one dummy, so the shell reads as a layer
coating the cluster's outside with its buried core visibly left uncoated.

Raw sample points for one exposed and one buried atom are overlaid in the
exposed/occluded colours, so the cloud can be seen coming *from* the occlusion
test rather than being asserted alongside it.

Split out from [`vis_sasa_hydro`](@ref) so the plot can be assembled, saved
or inspected without a window.
"""
function sasa_hydro_figure(; n_target::Union{Nothing,Int} = nothing, n_show::Int = 400, probe::Float64 = 1.4)
    sc = packed_cluster_scene(; probe)
    crds = Molecules.coords_cartesian(sc.mol)
    rads = Molecules.radii(sc.mol)
    natoms = size(crds, 2)
    area, exposed = SASA.sasa(sc.mol; probe)
    shell, shell_area, shell_cls = SASA.shell_points(sc.mol; probe, n_target)

    exp_idx = findall(exposed)
    bur_idx = findall(!, exposed)
    i_exposed = exp_idx[argmax(area[exp_idx])]                        # most-exposed atom
    highlight = isempty(bur_idx) ? (i_exposed,) : (i_exposed, first(bur_idx))

    fig = Figure(size = (1000, 900))
    ax = Axis3(
        fig[1, 1];
        title = @sprintf(
            "%s: hydration-shell dummy cloud\n%d dummies (%d convex, %d concave, %d cavity) over %d exposed atoms, %d buried",
            sc.title, size(shell, 2),
            count(==(SASA.CONVEX), shell_cls), count(==(SASA.CONCAVE), shell_cls),
            count(==(SASA.CAVITY), shell_cls), length(exp_idx), length(bur_idx)),
        titlesize = 15, aspect = :data, azimuth = 1.1π,
        xlabel = "x", ylabel = "y", zlabel = "z")

    # atom spheres: small, uniform, low-alpha.
    for i in 1:natoms
        ρ = Float32(rads[i] + probe)
        c = Point3f(crds[1, i], crds[2, i], crds[3, i])
        mesh!(ax, Sphere(c, ρ); color = (ATOM_COLOR, i in highlight ? 0.35 : 0.10),
                transparency = true, shading = NoShading)
    end

    # the shell itself: one marker per accessible point, i.e. one per dummy
    # `Scattering.hydration` will place. Every point of a given atom carries the
    # same area, so a uniform marker size is the honest rendering.
    bead_colour = Dict(SASA.CONVEX => DUMMY_COLOR, SASA.CONCAVE => CONCAVE_COLOR,
                        SASA.CAVITY => CAVITY_COLOR)
    scatter!(ax, [Point3f(shell[1, k], shell[2, k], shell[3, k]) for k in axes(shell, 2)];
            color = [bead_colour[c] for c in shell_cls], markersize = 5)

    # raw sample points, only on the two highlighted atoms
    for i in highlight
        st = atom_sample_points(sc.mol, i, n_show, probe)
        ex = [Point3f(p...) for (p, e) in zip(st.pts, st.exposed) if e]
        oc = [Point3f(p...) for (p, e) in zip(st.pts, st.exposed) if !e]
        isempty(oc) || scatter!(ax, oc; color = OCCLUDED_COLOR, markersize = 5)
        isempty(ex) || scatter!(ax, ex; color = EXPOSED_COLOR, markersize = 5)
    end

    els = [ MarkerElement(color = ATOM_COLOR, marker = :circle, markersize = 14),
            MarkerElement(color = DUMMY_COLOR, marker = :circle, markersize = 10),
            MarkerElement(color = CONCAVE_COLOR, marker = :circle, markersize = 10),
            MarkerElement(color = CAVITY_COLOR, marker = :circle, markersize = 10),
            MarkerElement(color = EXPOSED_COLOR, marker = :circle, markersize = 10),
            MarkerElement(color = OCCLUDED_COLOR, marker = :circle, markersize = 10)]
    Legend(fig[2, 1], els,
            [   "atom (expanded radius)",
                @sprintf("bead: convex (%.3f A^2 each)",
                        isempty(shell_area) ? 0.0 : first(shell_area)),
                "bead: concave", "bead: cavity",
                "sample point: exposed (highlighted atoms only)",
                "sample point: occluded (highlighted atoms only)"];
            orientation = :horizontal, framevisible = false, nbanks = 2, labelsize = 12)
    return fig
end

"""
    vis_sasa_hydro(; n_target = SASA.SHELL_POINTS, n_show = 400, probe = 1.4) -> Nothing

Display [`sasa_hydro_figure`](@ref): the whole hydration-shell dummy cloud on a
real ~80-atom packed cluster, plus the raw sample points on one buried and one
exposed atom to show the underlying mechanism. Blocks until the window is
closed.

# Keywords
-   `n_target`: global dummy budget; `nothing` derives it from accessible area.
                CRYSOL's `--fb` analogue, and what `Scattering.hydration` pays for.
-   `n_show`:   sample points drawn for the two highlighted atoms only.
-   `probe`:    solvent probe radius.
"""
vis_sasa_hydro(; n_target::Union{Nothing,Int} = nothing, n_show::Int = 400, probe::Float64 = 1.4) =
    wait(display(sasa_hydro_figure(; n_target, n_show, probe)))
