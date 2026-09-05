# Optional visual check for the Shrake-Rupley occlusion machinery in
# `src/Molecule/SASA.jl`. Kept out of the module so the geometry core stays
# free of plotting deps; GLMakie lives in visualize/Project.toml.
#
# Every scene uses `FixedRadii`, a test-only `RadiiSource` defined below, so the
# geometry is exact and reproducible instead of depending on the atomic-radii
# database.
#
# Run with:
#   julia --project=visualize -e 'include("test/tools/sasa_vis.jl"); vis_sasa_cases()'
#   julia --project=visualize -e 'include("test/tools/sasa_vis.jl"); vis_sasa_mesh_convergence()'
#   julia --project=visualize -e 'include("test/tools/sasa_vis.jl"); vis_sasa_molecule(cluster_scene().mol)'
#
# Numbers only, no window (works headless):
#   julia --project=visualize -e 'include("test/tools/sasa_vis.jl"); sasa_scene_report()'

using ScatterNet
using ScatterNet.Interfaces: Interfaces, RadiiSource
using ScatterNet.Molecule.Molecules: Molecules, Molecule
using ScatterNet.Molecule.SASA: SASA
using ScatterNet.Molecule.SASA.PlasticMap: plastic_points
using Printf: @printf, @sprintf
using GLMakie
using GLMakie.Makie: Tesselation   # not re-exported by GLMakie itself

# --------------------------------------------------------------------------
# radii source
# --------------------------------------------------------------------------

"""
Test-only [`RadiiSource`](@ref) mapping element labels to hand-picked radii, so
the demo scenes have the exact geometry the occlusion cases need rather than
whatever the atomic-radii database happens to hold.
"""
struct FixedRadii <: RadiiSource
    table::Dict{String,Float64}
end

"""
    lookup(src::FixedRadii, ions) -> Vector{Tuple{String,Union{Float64,Nothing}}}

Resolve each label against `src.table`, `nothing` for anything absent.
"""
function Interfaces.lookup(src::FixedRadii, ions)
    out = Vector{Tuple{String,Union{Float64,Nothing}}}(undef, length(ions))
    for (i, ion) in enumerate(ions)
        s = String(ion)
        out[i] = (s, get(src.table, s, nothing))
    end
    return out
end

# --------------------------------------------------------------------------
# pure geometry: point states, no plotting
# --------------------------------------------------------------------------

const EXPOSED_COLOR  = RGBf(1.00, 0.62, 0.13)   # warm/bright
const OCCLUDED_COLOR = RGBf(0.24, 0.26, 0.32)   # dark/desaturated

"""
    atom_point_states(mol, i, n, probe) -> NamedTuple

Shrake-Rupley state of every sample point on atom `i`: the `n` plastic-sequence
directions mapped onto the expanded sphere of radius `radii(mol)[i] + probe`,
each tested with `SASA._occluded` against every other atom.

The demo scenes are tiny, so the full index list is passed as the candidate set
instead of replicating `sasa_atoms`' KD-tree range query; `i` is passed as
`self` so the atom never occludes its own points.

This deliberately samples *every* atom, including the ones `sasa_atoms` would
never sample: `SASA._classify` resolves the fully-exposed and fully-engulfed
regimes exactly from the neighbour list, and those atoms reach no point loop at
all in the real code. Sampling them anyway is what makes the pictures worth
looking at -- you can see the state the exact predicate inferred. `status`
reports which path `sasa_atoms` actually took, so a panel can say so, and
`area` still agrees with `sasa_atoms(mol; n_occ = n, n_exp = n, probe = probe)[i]` in every case.

# Arguments
- `mol`: molecule to sample.
- `i`: atom index (column of `coords_cartesian(mol)`).
- `n`: number of sample points.
- `probe`: solvent probe radius.

# Returns
`(; pts, exposed, n_exposed, n_total, frac, area, rho, center, status, sampled)`
where `pts` is the `n` sample points, `exposed[j]` is `true` when `pts[j]` is
solvent accessible, `rho` is the expanded radius the points sit on, `status` is
a [`SASA.Coverage`](@ref) value, and
`sampled` is `false` when the exact pre-filter settles the atom on its own.
"""
function atom_point_states(mol::Molecule, i::Int, n::Int, probe::Float64)
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

    n_exposed = count(exposed)
    status = SASA._classify(i, cands, crds, rads, probe)
    return (;   pts, exposed, n_exposed, n_total = n, frac = n_exposed / n,
                area = 4π * ρ^2 * (n_exposed / n), rho = ρ,
                center = (cx, cy, cz), status,
                sampled = status == SASA.AMBIGUOUS)
end

"""
    two_sphere_exposed_area(ρ, d) -> Float64

Exact solvent-accessible area of one of two equal expanded spheres of radius `ρ`
whose centers are `d` apart, `0 < d < 2ρ`: the full sphere less the spherical
cap of height `ρ - d/2` that lies inside its partner,
`4πρ² - 2πρ(ρ - d/2)`.

Used as the ground truth the mesh-convergence panels are scored against.
"""
function two_sphere_exposed_area(ρ::Float64, d::Float64)
    0.0 < d < 2ρ || throw(DomainError(d, "need 0 < d < 2ρ for an overlapping pair"))
    return 4π * ρ^2 - 2π * ρ * (ρ - d / 2)
end

"""
    two_sphere_exposed_frac(ρ, d) -> Float64

[`two_sphere_exposed_area`](@ref) as a fraction of the full sphere `4πρ²`.
"""
two_sphere_exposed_frac(ρ::Float64, d::Float64) = two_sphere_exposed_area(ρ, d) / (4π * ρ^2)

# --------------------------------------------------------------------------
# the four scenes
# --------------------------------------------------------------------------

"""
    exposed_scene(; probe = 1.4) -> NamedTuple

Two equal atoms whose expanded spheres just fail to touch (`r = 1.5`, so
`ρ = 2.9` at the default probe, centers `6.2` apart): the neighbour is a
candidate but occludes nothing, so every sample point on both atoms is exposed.

`focus` lists the atoms whose points get drawn.
"""
function exposed_scene(; probe::Float64 = 1.4)
    src = FixedRadii(Dict("A" => 1.5))
    mol = Molecules.create("fully exposed", ["A", "A"],
                            [(0.0, 0.0, 0.0), (6.2, 0.0, 0.0)]; radii_source = src)
    return (; mol, probe, focus = [1, 2], title = "1. Fully exposed")
end

"""
    partial_scene(; probe = 1.4, d = 3.5) -> NamedTuple

Two equal overlapping atoms (`r = 1.5`, `ρ = 2.9`, centers `d` apart) — the
occluded points form one clean spherical cap facing the partner, everything else
is exposed. Carries `analytic_frac`, the exact cap-derived exposed fraction from
[`two_sphere_exposed_frac`](@ref).
"""
function partial_scene(; probe::Float64 = 1.4, d::Float64 = 3.5)
    src = FixedRadii(Dict("A" => 1.5))
    mol = Molecules.create("partially occluded", ["A", "A"],
                            [(0.0, 0.0, 0.0), (d, 0.0, 0.0)]; radii_source = src)
    ρ = 1.5 + probe
    return (;   mol, probe, focus = [1], title = "2. Partially occluded",
                d, rho = ρ, analytic_frac = two_sphere_exposed_frac(ρ, d),
                analytic_area = two_sphere_exposed_area(ρ, d),
                azimuth = 0.20π)  # look down the +x side so the cap faces the camera
end

"""
    buried_scene(; probe = 1.4) -> NamedTuple

A small atom (`r = 0.5`, `ρ = 1.9`) sitting `3.0` from the center of a big one
(`r = 6.0`, `ρ = 7.4`). `3.0 + 1.9 = 4.9 < 7.4` with `2.5` Å of margin, so the
small atom's whole expanded sphere is strictly inside the big one's and every
sample point is occluded.
"""
function buried_scene(; probe::Float64 = 1.4)
    src = FixedRadii(Dict("BIG" => 6.0, "SML" => 0.5))
    mol = Molecules.create("fully buried", ["BIG", "SML"],
                            [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0)]; radii_source = src)
    return (; mol, probe, focus = [2], title = "3. Fully buried")
end

"""
    cluster_scene(; probe = 1.4, spacing = 4.0) -> NamedTuple

Seven equal atoms (`r = 1.5`): one at the origin caged by six octahedral
neighbours at `±spacing` along each axis. The caged atom keeps only the small
patches that peek between neighbours while the shell atoms stay mostly exposed,
so a single scene shows the whole spread of per-atom states.
"""
function cluster_scene(; probe::Float64 = 1.4, spacing::Float64 = 4.0)
    src = FixedRadii(Dict("A" => 1.5))
    s = spacing
    coords = [  (0.0, 0.0, 0.0),
                ( s, 0.0, 0.0), (-s, 0.0, 0.0),
                (0.0,  s, 0.0), (0.0, -s, 0.0),
                (0.0, 0.0,  s), (0.0, 0.0, -s)]
    mol = Molecules.create("cluster", fill("A", 7), coords; radii_source = src)
    return (; mol, probe, focus = collect(1:7), title = "4. Small cluster")
end

"""
    sasa_scenes(; probe = 1.4) -> Vector

The four demo scenes in panel order: fully exposed, partially occluded, fully
buried, small cluster.
"""
sasa_scenes(; probe::Float64 = 1.4) =
    [   exposed_scene(; probe), partial_scene(; probe), buried_scene(; probe),
        cluster_scene(; probe)]

# --------------------------------------------------------------------------
# headless report
# --------------------------------------------------------------------------

"""
    classification_label(status) -> String

Human-readable name for a [`SASA._classify`](@ref) verdict, noting whether
`sasa_atoms` samples that atom or settles it exactly from the neighbour list.
"""
classification_label(status::SASA.Coverage) =
    status == SASA.ALL_EXPOSED ?    "exact: no neighbour reaches it" :
    status == SASA.ALL_BURIED  ?    "exact: engulfed by a neighbour" :
                                    "sampled: caps overlap"

"""
    sasa_scene_report(; n = 512, probe = 1.4, ns = (64, 256, 1024, 4096)) -> Nothing

Print every number the figures annotate — per-atom exposed counts, fractions and
areas for the four scenes, the analytic cross-check on the overlapping pair, and
the mesh-refinement sweep — without opening a window. Handy on a headless box,
and the thing to run first if a panel ever looks wrong.
"""
function sasa_scene_report(;n::Int = 512, probe::Float64 = 1.4,
                            ns = (64, 256, 1024, 4096))
    for sc in sasa_scenes(; probe)
        println(sc.title, "  (probe = ", probe, ", n = ", n, ")")
        areas = SASA.sasa_atoms(sc.mol; n_occ = n, n_exp = n, probe = probe)
        for i in sc.focus
            st = atom_point_states(sc.mol, i, n, probe)
            @printf("  atom %d: exposed %4d/%4d  frac %.4f  area %8.3f  (sasa_atoms %8.3f)  [%s]\n",
                    i, st.n_exposed, st.n_total, st.frac, st.area, areas[i],
                    classification_label(st.status))
        end
        if haskey(sc, :analytic_frac)
            @printf("  analytic: frac %.6f  area %8.3f\n", sc.analytic_frac, sc.analytic_area)
        end
        @printf("  total SASA %.3f A^2\n\n", sum(areas))
    end

    sc = partial_scene(; probe)
    @printf("Mesh refinement on scene 2 (analytic frac %.6f, area %.3f)\n",
            sc.analytic_frac, sc.analytic_area)
    for k in ns
        st = atom_point_states(sc.mol, 1, k, probe)
        @printf("  n = %5d  frac %.6f  err %+.6f  area %8.3f  err %+.4f\n",
                k, st.frac, st.frac - sc.analytic_frac, st.area,
                st.area - sc.analytic_area)
    end
    return nothing
end

# --------------------------------------------------------------------------
# plotting helpers
# --------------------------------------------------------------------------

"""
    _draw_atoms!(ax, mol, probe, focus) -> Nothing

Draw every atom of `mol` at its expanded radius `r + probe`, so the sample points
are visibly sitting *on* a surface rather than floating.

Atoms in `focus` (the ones being sampled) get a faint wireframe shell — a filled
surface in front of them would wash out the point colours, and the wireframe
still pins the points to a sphere. Everything else is a neutral low-alpha solid,
since for an occluder it is the enclosed *volume* that explains the verdict.
"""
function _draw_atoms!(ax, mol::Molecule, probe::Float64, focus)
    crds = Molecules.coords_cartesian(mol)
    rads = Molecules.radii(mol)
    wire_alpha = length(focus) > 2 ? 0.18 : 0.35   # busy scenes need a fainter cage
    for i in axes(crds, 2)
        ρ = Float32(rads[i] + probe)
        c = Point3f(crds[1, i], crds[2, i], crds[3, i])
        if i in focus
            wireframe!( ax, Tesselation(Sphere(c, ρ), 24);
                        color = (RGBf(0.30, 0.55, 0.85), wire_alpha), linewidth = 0.5,
                        transparency = true)
        else
            mesh!(  ax, Sphere(c, ρ); color = (RGBf(0.55, 0.58, 0.64), 0.16),
                    transparency = true, shading = NoShading)
        end
    end
    return nothing
end

"""
    _draw_points!(ax, st; markersize) -> Nothing

Scatter one atom's sample points at their true positions on the expanded sphere,
coloured by state: [`EXPOSED_COLOR`](@ref) when solvent accessible,
[`OCCLUDED_COLOR`](@ref) when swallowed by a neighbour.
"""
function _draw_points!(ax, st; markersize = 6)
    ex = [Point3f(p...) for (p, e) in zip(st.pts, st.exposed) if e]
    oc = [Point3f(p...) for (p, e) in zip(st.pts, st.exposed) if !e]
    isempty(oc) || scatter!(ax, oc; color = OCCLUDED_COLOR, markersize = markersize)
    isempty(ex) || scatter!(ax, ex; color = EXPOSED_COLOR, markersize = markersize)
    return nothing
end

"""
    _state_legend!(fig, cell) -> Nothing

Put the shared exposed/occluded point-colour key into `fig[cell...]`.
"""
function _state_legend!(fig, cell)
    els = [ MarkerElement(color = EXPOSED_COLOR, marker = :circle, markersize = 12),
            MarkerElement(color = OCCLUDED_COLOR, marker = :circle, markersize = 12)]
    Legend(fig[cell...], els, ["exposed point", "occluded point"];
           orientation = :horizontal, framevisible = false)
    return nothing
end

# --------------------------------------------------------------------------
# entry points
# --------------------------------------------------------------------------

"""
    sasa_cases_figure(; n = 512, probe = 1.4) -> Figure

Build (but do not display) the four-regime figure. Split out from
[`vis_sasa_cases`](@ref) so the plot can be assembled, saved or inspected without
a window.
"""
function sasa_cases_figure(; n::Int = 512, probe::Float64 = 1.4)
    scenes = sasa_scenes(; probe)
    fig = Figure(size = (1400, 950))
    Label(  fig[0, 1:2],
            "Shrake-Rupley occlusion states  (probe = $(probe) A, n = $(n) points/atom)";
            fontsize = 20, font = :bold)

    for (k, sc) in enumerate(scenes)
        row = (k - 1) ÷ 2 + 1
        col = (k - 1) % 2 + 1

        sts = [atom_point_states(sc.mol, i, n, probe) for i in sc.focus]
        n_exp = sum(s -> s.n_exposed, sts)
        n_tot = sum(s -> s.n_total, sts)
        area  = sum(s -> s.area, sts)

        sub =   length(sc.focus) == 1 ?
                @sprintf(   "%d/%d exposed (%.3f)  area %.2f A^2",
                            n_exp, n_tot, n_exp / n_tot, area) :
                @sprintf(   "%d/%d exposed (%.3f)  total area %.2f A^2 over %d atoms",
                            n_exp, n_tot, n_exp / n_tot, area, length(sc.focus))
        extra = haskey(sc, :analytic_frac) ?
                @sprintf(   "\nanalytic cap fraction %.4f  (error %+.4f)",
                            sc.analytic_frac, n_exp / n_tot - sc.analytic_frac) : ""
        # say which branch sasa_atoms takes: the points shown are illustrative
        # for the two regimes the exact pre-filter resolves without sampling.
        paths = unique(classification_label(s.status) for s in sts)
        extra *= "\n" * join(paths, " + ")

        ax = Axis3( fig[row, col]; title = sc.title * "\n" * sub * extra,
                    titlesize = 14, aspect = :data,
                    azimuth = get(sc, :azimuth, 1.275π),
                    xlabel = "x", ylabel = "y", zlabel = "z")
        _draw_atoms!(ax, sc.mol, probe, sc.focus)
        for st in sts
            _draw_points!(ax, st; markersize = length(sc.focus) > 2 ? 4 : 6)
        end
    end

    _state_legend!(fig, (3, 1:2))
    return fig
end

"""
    vis_sasa_cases(; n = 512, probe = 1.4) -> Nothing

Four-panel figure over the distinct occlusion regimes — fully exposed, partially
occluded, fully buried, and a small cluster — with every sample point coloured by
its `SASA._occluded` verdict and each title carrying the exposed count, total,
fraction and resulting area. Blocks until the window is closed.

# Keywords
- `n`: sample points per atom.
- `probe`: solvent probe radius.
"""
vis_sasa_cases(; n::Int = 512, probe::Float64 = 1.4) =
    wait(display(sasa_cases_figure(; n, probe)))

"""
    sasa_mesh_convergence_figure(; ns = (64, 256, 1024, 4096), probe = 1.4, d = 3.5) -> Figure

Build (but do not display) the mesh-refinement figure; see
[`vis_sasa_mesh_convergence`](@ref).
"""
function sasa_mesh_convergence_figure(; ns = (64, 256, 1024, 4096),
                                        probe::Float64 = 1.4, d::Float64 = 3.5)
    sc = partial_scene(; probe, d)
    ncol = length(ns) <= 2 ? length(ns) : cld(length(ns), 2)
    fig = Figure(size = (360 * ncol + 60, 900))
    Label(fig[0, 1:ncol],
        @sprintf(
            "Mesh refinement on one overlapping pair  (rho = %.2f A, d = %.2f A)\nanalytic exposed fraction %.6f, area %.3f A^2",
            sc.rho, sc.d, sc.analytic_frac, sc.analytic_area);
        fontsize = 18, font = :bold)

    for (k, n) in enumerate(ns)
        row = (k - 1) ÷ ncol + 1
        col = (k - 1) % ncol + 1
        st = atom_point_states(sc.mol, 1, Int(n), probe)
        title = @sprintf(   "n = %d\nfrac %.4f  (err %+.4f)\narea %.2f  (err %+.3f)",
                            n, st.frac, st.frac - sc.analytic_frac,
                            st.area, st.area - sc.analytic_area)
        ax = Axis3( fig[row, col]; title, titlesize = 13, aspect = :data,
                    azimuth = get(sc, :azimuth, 1.275π),
                    xlabel = "x", ylabel = "y", zlabel = "z")
        _draw_atoms!(ax, sc.mol, probe, sc.focus)
        _draw_points!(ax, st; markersize = n <= 256 ? 9 : (n <= 1024 ? 6 : 3))
    end

    _state_legend!(fig, (cld(length(ns), ncol) + 1, 1:ncol))
    return fig
end

"""
    vis_sasa_mesh_convergence(; ns = (64, 256, 1024, 4096), probe = 1.4, d = 3.5) -> Nothing

The *same* two-overlapping-atom scene rendered once per entry of `ns`, so the
effect of mesh density is visible directly: coarse meshes give a ragged cap edge
and a noisy fraction, fine meshes converge on the analytic value
[`two_sphere_exposed_area`](@ref). Each panel is annotated with `n`, the measured
exposed fraction and the signed error against that exact value.

This is the picture behind `sasa_atoms`' two-pass split: a handful of points is
already enough to answer "is this atom buried at all?" (the `n_occ` pass), while
a fraction accurate to a percent needs an order of magnitude more (the `n_exp`
pass). Blocks until the window is closed.

# Keywords
- `ns`: point counts, one panel each.
- `probe`: solvent probe radius.
- `d`: center separation of the overlapping pair.
"""
vis_sasa_mesh_convergence(; ns = (64, 256, 1024, 4096), probe::Float64 = 1.4,
                            d::Float64 = 3.5) =
    wait(display(sasa_mesh_convergence_figure(; ns, probe, d)))

"""
    sasa_molecule_figure(mol; n_occ = 128, n_exp = 1024, probe = 1.4) -> Figure

Build (but do not display) the per-atom SASA figure; see
[`vis_sasa_molecule`](@ref).
"""
function sasa_molecule_figure(
    mol::Molecule; 
    n_occ::Int = 128, 
    n_exp::Int = 1024,
    probe::Float64 = 1.4
)
    areas = SASA.sasa_atoms(mol; n_occ = n_occ, n_exp = n_exp, probe = probe)
    crds = Molecules.coords_cartesian(mol)
    rads = Molecules.radii(mol)
    lo, hi = extrema(areas)
    hi = hi > lo ? hi : lo + 1.0

    fig = Figure(size = (900, 800))
    ax = Axis3(
        fig[1, 1];
        title = @sprintf(
            "%s: per-atom SASA  (probe = %.2f A, n_occ = %d, n_exp = %d)\ntotal %.2f A^2 over %d atoms",
            Molecules.name(mol), probe, n_occ, n_exp,
            sum(areas), length(areas)),
            titlesize = 14, aspect = :data,
            xlabel = "x", ylabel = "y", zlabel = "z")

    cmap = to_colormap(:inferno)
    for i in axes(crds, 2)
        t = (areas[i] - lo) / (hi - lo)
        col = cmap[clamp(round(Int, t * (length(cmap) - 1)) + 1, 1, length(cmap))]
        mesh!(ax, Sphere(Point3f(crds[1, i], crds[2, i], crds[3, i]), Float32(rads[i] + probe)); color = (col, 0.85), transparency = true)
    end
    Colorbar(fig[1, 2]; colormap = :inferno, limits = (lo, hi), label = "per-atom SASA (A^2)")
    return fig
end

"""
    vis_sasa_molecule(mol; n_occ = 128, n_exp = 1024, probe = 1.4) -> Nothing

Colour a whole molecule's atoms by their per-atom area from
`SASA.sasa_atoms`: each atom is drawn as a sphere at its expanded radius
`r + probe`, tinted from buried (dark) to fully exposed (bright), with a
colourbar in absolute Å². Blocks until the window is closed.

# Arguments
- `mol`: molecule to score and draw.

# Keywords
- `n_occ`: points for the coarse buried/not-buried pass.
- `n_exp`: points for the finer exposed-fraction pass.
- `probe`: solvent probe radius.
"""
vis_sasa_molecule(mol::Molecule; n_occ::Int = 128, n_exp::Int = 1024, probe::Float64 = 1.4) =
    wait(display(sasa_molecule_figure(mol; n_occ, n_exp, probe)))
