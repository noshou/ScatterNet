# Accuracy + performance harness for `src/Molecule/SASA.jl`. Lives here rather
# than under test/ because the summary figures need GLMakie; nothing in this
# file runs under `Pkg.test()`.
#
# Every system is built from `BenchRadii`, a test-only `RadiiSource` holding
# Bondi/Alvarez van der Waals radii, so geometry is exact and reproducible
# instead of depending on the atomic-radii database.
#
# Numbers only, no window (works headless):
#   julia --project=visualize -e 'include("visualize/sasa_bench.jl"); bench_report()'
#
# Write the raw sweep to CSV for external analysis:
#   julia --project=visualize -e 'include("visualize/sasa_bench.jl"); accuracy_sweep("sweep.csv")'
#
# Summary figures, in a window or straight to PNG (the latter works headless):
#   julia --project=visualize -e 'include("visualize/sasa_bench.jl"); vis_bench()'
#   julia --project=visualize -e 'include("visualize/sasa_bench.jl"); save_bench_figures("bench")'

using ScatterNet
using ScatterNet.Interfaces: Interfaces, RadiiSource
using ScatterNet.Molecule.Molecules: Molecules, Molecule
using ScatterNet.Molecule.SASA: SASA
using Random, Printf, Statistics, DelimitedFiles
using GLMakie

# ---------------------------------------------------------------------------
# radii
# ---------------------------------------------------------------------------

"Test-only [`RadiiSource`](@ref) backed by a fixed element -> radius table."
struct BenchRadii <: RadiiSource
    table::Dict{String,Float64}
end
Interfaces.lookup(s::BenchRadii, ions::AbstractVector{<:AbstractString}) =
    Tuple{String,Union{Float64,Nothing}}[
        (String(i), get(s.table, String(i), nothing)) for i in ions
    ]

"Bondi/Alvarez van der Waals radii (Å); `*_i` entries are Shannon ionic radii."
const BENCH_RAD = Dict(
    "H"=>1.10, "C"=>1.70, "N"=>1.55, "O"=>1.52, "F"=>1.47, "P"=>1.80, "S"=>1.80,
    "Cl"=>1.75, "Br"=>1.85, "I"=>1.98, "Si"=>2.10, "B"=>1.92,
    "Na"=>2.27, "K"=>2.75, "Cs"=>3.43, "Mg"=>1.73, "Ca"=>2.31,
    "Fe"=>2.04, "Cu"=>1.96, "Zn"=>2.01, "Ni"=>1.97, "Pt"=>2.13, "Ru"=>2.10,
    "Mo"=>2.10, "Ti"=>2.11,
    "Na_i"=>1.16, "Cl_i"=>1.81, "Cs_i"=>1.81, "O_i"=>1.40)

const BENCH_SRC = BenchRadii(BENCH_RAD)

"Build a molecule from the shared radii table."
bench_mol(name, els, crds) = Molecules.create(name, els, crds; radii_source = BENCH_SRC)

const BENCH_PROBE = 1.4
"Working point count for the sweeps (fast, but past the steep part of the curve)."
const BENCH_NEXP  = 2048
"Converged reference point count."
const BENCH_NREF  = 32768
"`area_tol` values swept, including 0 (never skip)."
const BENCH_TOLS  = [0.0, 0.05, 0.2, 0.5, 0.8, 2.0, 5.0]

# ---------------------------------------------------------------------------
# geometry builders
# ---------------------------------------------------------------------------

_grid(n) = vec([(i, j, k) for i in 0:n-1, j in 0:n-1, k in 0:n-1])

"Simple-cubic lattice, `n³` sites, lattice constant `a`."
sc_lattice(n, a) = [(a*i, a*j, a*k) for (i, j, k) in _grid(n)]

"Body-centred cubic lattice."
function bcc_lattice(n, a)
    p = sc_lattice(n, a)
    append!(p, [(a*(i+.5), a*(j+.5), a*(k+.5)) for (i, j, k) in _grid(n)])
    p
end

"Face-centred cubic lattice."
function fcc_lattice(n, a)
    p = sc_lattice(n, a)
    for (dx, dy, dz) in ((.5,.5,0.), (.5,0.,.5), (0.,.5,.5))
        append!(p, [(a*(i+dx), a*(j+dy), a*(k+dz)) for (i, j, k) in _grid(n)])
    end
    p
end

"Diamond-cubic lattice (FCC plus its quarter-diagonal copy)."
function diamond_lattice(n, a)
    f = fcc_lattice(n, a)
    vcat(f, [(x+a/4, y+a/4, z+a/4) for (x, y, z) in f])
end

"Rock-salt (NaCl) binary lattice -> `(coords, elements)`."
function rocksalt(n, a, eA, eB)
    p = NTuple{3,Float64}[]; e = String[]
    for (i, j, k) in _grid(2n)
        push!(p, (a/2*i, a/2*j, a/2*k)); push!(e, iseven(i+j+k) ? eA : eB)
    end
    p, e
end

"Caesium-chloride binary lattice -> `(coords, elements)`."
function cscl(n, a, eA, eB)
    p = sc_lattice(n, a); e = fill(eA, length(p))
    append!(p, [(a*(i+.5), a*(j+.5), a*(k+.5)) for (i, j, k) in _grid(n)])
    append!(e, fill(eB, n^3))
    p, e
end

"Zig-zag `nc`-carbon alkane chain with its hydrogens -> `(coords, elements)`."
function alkane(nc)
    p = NTuple{3,Float64}[]; e = String[]
    for i in 0:nc-1
        x = 1.26i; y = iseven(i) ? 0.0 : 0.51
        push!(p, (x, y, 0.0));          push!(e, "C")
        push!(p, (x, y+0.63,  0.90));   push!(e, "H")
        push!(p, (x, y+0.63, -0.90));   push!(e, "H")
    end
    p, e
end

"Planar benzene ring with hydrogens."
function benzene()
    p = NTuple{3,Float64}[]; e = String[]
    for k in 0:5
        θ = 2π*k/6
        push!(p, (1.39cos(θ), 1.39sin(θ), 0.0)); push!(e, "C")
        push!(p, (2.48cos(θ), 2.48sin(θ), 0.0)); push!(e, "H")
    end
    p, e
end

"Chair cyclohexane with hydrogens."
function cyclohexane()
    p = NTuple{3,Float64}[]; e = String[]
    for k in 0:5
        θ = 2π*k/6; z = iseven(k) ? 0.25 : -0.25
        push!(p, (1.46cos(θ), 1.46sin(θ), z));                       push!(e, "C")
        push!(p, (2.55cos(θ), 2.55sin(θ), 1.1z));                    push!(e, "H")
        push!(p, (1.30cos(θ), 1.30sin(θ), z + (iseven(k) ? 1.1 : -1.1))); push!(e, "H")
    end
    p, e
end

"Adamantane C₁₀ cage (carbons only) -- a compact, heavily self-occluding case."
adamantane() = ([(0.,0.,0.), (1.54,1.54,0.), (1.54,0.,1.54), (0.,1.54,1.54),
                 (2.57,0.90,0.90), (0.90,2.57,0.90), (0.90,0.90,2.57),
                 (2.57,2.57,2.57), (3.47,1.80,1.80), (1.80,3.47,1.80)],
                fill("C", 10))

"Octahedral ML₆ complex at metal-ligand distance `d`."
function octahedral(metal, ligand, d)
    p = [(0., 0., 0.)]; e = [metal]
    for v in ((d,0,0), (-d,0,0), (0,d,0), (0,-d,0), (0,0,d), (0,0,-d))
        push!(p, Float64.(v)); push!(e, ligand)
    end
    p, e
end

"Ferrocene-like sandwich: one Fe between two eclipsed C₅ rings."
function ferrocene()
    p = [(0., 0., 0.)]; e = ["Fe"]
    for s in (1.66, -1.66), k in 0:4
        θ = 2π*k/5 + (s > 0 ? 0.0 : π/5)
        push!(p, (1.21cos(θ), 1.21sin(θ), s)); push!(e, "C")
    end
    p, e
end

"Roughly spherical FCC metal nanoparticle carved from an `n³` cell."
function metal_cluster(n, a, el)
    pts = fcc_lattice(n, a)
    ctr = ntuple(t -> sum(p[t] for p in pts)/length(pts), 3)
    keep = [p for p in pts if sqrt(sum((p[t]-ctr[t])^2 for t in 1:3)) <= a*n*0.45]
    keep, fill(el, length(keep))
end

"`n` atoms placed uniformly at random in a cube chosen for number density `dens`."
function random_pack(n, dens, el; seed = 1)
    box = (n/dens)^(1/3); rng = MersenneTwister(seed)
    [(box*rand(rng), box*rand(rng), box*rand(rng)) for _ in 1:n], fill(el, n)
end

# ---------------------------------------------------------------------------
# diagnostics
# ---------------------------------------------------------------------------

"""
    classify_census(m, probe) -> NamedTuple

How many atoms `SASA._classify` resolves exactly (`exposed`, `buried`) versus
hands to point sampling (`ambiguous`).
"""
function classify_census(m::Molecule, probe::Float64)
    rads = Molecules.radii(m); rmax = Molecules.r_max(m)
    crds = Molecules.coords_cartesian(m); tree = SASA.KDTree(crds)
    c = Dict(SASA.ALL_EXPOSED => 0, SASA.ALL_BURIED => 0, SASA.AMBIGUOUS => 0)
    for i in axes(crds, 2)
        ρ = rads[i] + probe
        cand = SASA.inrange(tree, @view(crds[:, i]), ρ + rmax + probe)
        c[SASA._classify(i, cand, crds, rads, probe)] += 1
    end
    (exposed = c[SASA.ALL_EXPOSED], buried = c[SASA.ALL_BURIED],
     ambiguous = c[SASA.AMBIGUOUS])
end

"""
    skip_diagnostics(m; probe, n_occ, n_exp, area_tol) -> NamedTuple

How much the `area_tol` early exit changed the answer for one molecule.

A skip is detected purely through the public API: an atom scoring `> 0` at
`area_tol = 0` but exactly `0.0` at the given tolerance is one the witness pass
wrote off. That keeps this analysis from drifting out of sync with the source.
"""
function skip_diagnostics(
    m::Molecule; probe = BENCH_PROBE, n_occ = 512, n_exp = BENCH_NEXP, area_tol = 0.8)
    ref = SASA.sasa_atoms(m; probe, n_occ, n_exp, area_tol = 0.0)
    got = SASA.sasa_atoms(m; probe, n_occ, n_exp, area_tol)
    fired = [i for i in eachindex(ref) if ref[i] > 0 && got[i] == 0.0]
    lost  = isempty(fired) ? 0.0 : sum(ref[i] for i in fired)
    (   n_atoms = length(ref), n_skipped = length(fired), area_lost = lost,
        total = sum(got), ref_total = sum(ref),
        rel_error = sum(ref) == 0 ? 0.0 : (sum(got) - sum(ref))/sum(ref),
        worst_atom = isempty(fired) ? 0.0 : maximum(ref[i] for i in fired))
end

# ---------------------------------------------------------------------------
# the system catalogue
# ---------------------------------------------------------------------------

"""
    bench_systems() -> Vector{Tuple{String,String,Molecule}}

Every `(tag, class, molecule)` the sweeps run over: monoatomic radius scans,
dimer separation scans, four elemental and two binary crystal lattices swept
across packing density, organic molecules, organometallics, and random packings.
"""
function bench_systems()
    out = Tuple{String,String,Molecule}[]

    # 1. monoatomic, 40 radii from H-like to Cs-like
    for r in range(0.30, 3.00; length = 40)
        src = BenchRadii(Dict("X" => r))
        m = Molecules.create("lone", ["X"], [(0.,0.,0.)]; radii_source = src)
        push!(out, (@sprintf("lone_r%.2f", r), "monoatomic", m))
    end

    # 2. dimers, 40 separations each, at three radius ratios
    for (lbl, ra, rb) in (("1to1", 1.70, 1.70), ("2to1", 2.00, 1.00), ("5to1", 2.50, 0.50))
        src = BenchRadii(Dict("A" => ra, "B" => rb))
        for d in range(0.0, 2*(max(ra, rb) + BENCH_PROBE)*1.15; length = 40)
            m = Molecules.create("dim", ["A","B"], [(0.,0.,0.), (d,0.,0.)]; radii_source = src)
            push!(out, (@sprintf("dimer_%s_d%.2f", lbl, d), "dimer", m))
        end
    end

    # 3. elemental crystals, 40 packing densities each
    for (nm, build, el, nn) in (("SC", sc_lattice, "C", 4), ("BCC", bcc_lattice, "Fe", 3),
                                ("FCC", fcc_lattice, "Cu", 3), ("diamond", diamond_lattice, "Si", 2))
        r = BENCH_RAD[el]
        for a in range(1.2r, 5.0r; length = 40)
            pts = build(nn, a)
            push!(out, (@sprintf("%s_a%.2f", nm, a), "crystal", bench_mol(nm, fill(el, length(pts)), pts)))
        end
    end

    # 4. binary ionic crystals, 40 densities each
    for (nm, f, eA, eB) in (("NaCl", rocksalt, "Na_i", "Cl_i"), ("CsCl", cscl, "Cs_i", "Cl_i"))
        for a in range(2.0, 9.0; length = 40)
            pts, els = f(2, a, eA, eB)
            push!(out, (@sprintf("%s_a%.2f", nm, a), "crystal_binary", bench_mol(nm, els, pts)))
        end
    end

    # 5. organic
    for (nm, (p, e)) in (
        ("alkane_C6", alkane(6)), ("alkane_C12", alkane(12)),
        ("benzene", benzene()), ("cyclohexane", cyclohexane()),
        ("adamantane", adamantane()))
        push!(out, (nm, "organic", bench_mol(nm, e, p)))
    end

    # 6. organometallic / inorganic coordination
    for (nm, (p, e)) in (   ("ferrocene", ferrocene()),
                            ("FeCl6", octahedral("Fe", "Cl", 2.4)),
                            ("PtCl6", octahedral("Pt", "Cl", 2.3)),
                            ("MoO6",  octahedral("Mo", "O",  1.95)),
                            ("Ru_cluster", metal_cluster(3, 3.8, "Ru")),
                            ("Cu_cluster", metal_cluster(3, 3.6, "Cu")))
        push!(out, (nm, "organometallic", bench_mol(nm, e, p)))
    end

    # 7. random packings, 40 number densities
    for (i, dens) in enumerate(range(0.002, 0.09; length = 40))
        pts, els = random_pack(120, dens, "C"; seed = 100 + i)
        push!(out, (@sprintf("rand_d%.4f", dens), "random", bench_mol("r", els, pts)))
    end
    out
end

# ---------------------------------------------------------------------------
# sweeps
# ---------------------------------------------------------------------------

"""
    accuracy_sweep(csv_path = nothing; n_exp, n_occ, probe, tols) -> (header, rows)

Every system in [`bench_systems`](@ref) crossed with every `area_tol`, scoring
each against two references: the same `n_exp` at `area_tol = 0` (isolating the
skip error) and `BENCH_NREF` points at `area_tol = 0` (isolating sampling
error). Writes CSV when given a path.
"""
function accuracy_sweep(csv_path = nothing; n_exp = BENCH_NEXP, n_occ = 512,
                        probe = BENCH_PROBE, tols = BENCH_TOLS)
    rows = Vector{Any}[]
    for (tag, class, m) in bench_systems()
        ref  = SASA.sasa_atoms(m; probe, n_occ, n_exp, area_tol = 0.0)
        conv = SASA.sasa_atoms(m; probe, n_occ, n_exp = BENCH_NREF, area_tol = 0.0)
        cen  = classify_census(m, probe)
        sref, sconv = sum(ref), sum(conv)
        for tol in tols
            a = SASA.sasa_atoms(m; probe, n_occ, n_exp, area_tol = tol)
            fired = [i for i in eachindex(ref) if ref[i] > 0 && a[i] == 0.0]
            push!(rows, Any[tag, class, length(ref), tol, sum(a), sref, sconv,
                            length(fired),
                            isempty(fired) ? 0.0 : sum(ref[i] for i in fired),
                            isempty(fired) ? 0.0 : maximum(ref[i] for i in fired),
                            sref == 0 ? 0.0 : (sum(a) - sref)/sref,
                            sconv == 0 ? 0.0 : (sref - sconv)/sconv,
                            cen.exposed, cen.buried, cen.ambiguous])
        end
    end
    hdr = [ "tag","class","natoms","area_tol","total","ref_total","conv_total",
            "n_skipped","area_lost","worst_skipped_atom","rel_vs_ref",
            "rel_ref_vs_conv","n_exposed","n_buried","n_ambiguous"]
    if csv_path !== nothing
        open(csv_path, "w") do io
            writedlm(io, permutedims(hdr), ',')
            writedlm(io, permutedims(hcat(rows...)), ',')
        end
    end
    hdr, rows
end

"""
    convergence_sweep(; ns, probe) -> (header, rows)

Relative error against a `BENCH_NREF`-point reference as `n_exp` grows, for one
representative system per class.
"""
function convergence_sweep(;ns = (64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384),
                            probe = BENCH_PROBE)
    reps = filter(  t -> t[1] in ("SC_a2.04", "FCC_a2.35", "NaCl_a4.59",
                                "adamantane", "ferrocene", "rand_d0.0450"),
                    bench_systems())
    isempty(reps) && (reps = bench_systems()[1:6])
    rows = Vector{Any}[]
    for (tag, class, m) in reps
        conv = sum(SASA.sasa_atoms(m; probe, n_occ = 512, n_exp = BENCH_NREF, area_tol = 0.0))
        for n in ns
            v = sum(SASA.sasa_atoms(m; probe, n_occ = min(512, n), n_exp = n, area_tol = 0.0))
            push!(rows, Any[tag, class, n, v, conv, conv == 0 ? 0.0 : (v - conv)/conv])
        end
    end
    ["tag","class","n_exp","total","conv_total","rel_err"], rows
end

"""
    analytic_sweep(; probe) -> (header, rows)

Ground-truth check with no sampling reference at all: two EQUAL spheres of
expanded radius ρ at separation `d` have exposed area `4πρ² - 2πρ(ρ - d/2)`
each. Swept over radii and separations.
"""
function analytic_sweep(; probe = BENCH_PROBE, n_exp = 8192)
    rows = Vector{Any}[]
    for r in (0.5, 1.0, 1.5, 1.7, 2.1, 2.5, 3.0)
        ρ = r + probe
        for frac in range(0.08, 0.98; length = 12)
            d = 2ρ*frac
            src = BenchRadii(Dict("A" => r))
            m = Molecules.create("pair", ["A","A"], [(0.,0.,0.), (d,0.,0.)]; radii_source = src)
            exact = 4π*ρ^2 - 2π*ρ*(ρ - d/2)
            got = SASA.sasa_atoms(m; probe, n_occ = 512, n_exp, area_tol = 0.0)[1]
            push!(rows, Any[r, ρ, d, d/(2ρ), got, exact, (got - exact)/exact])
        end
    end
    ["radius","rho","d","d_over_2rho","measured","analytic","rel_err"], rows
end

"""
    speed_sweep(; probe) -> (header, rows)

Wall time against atom count, `n_exp`, and `area_tol`, plus the derived
per-atom cost. Uses the minimum of several runs to suppress GC noise.
"""
function speed_sweep(; probe = BENCH_PROBE, reps = 3)
    rows = Vector{Any}[]
    best(f) = minimum(@elapsed(f()) for _ in 1:reps)
    for n in (2, 3, 4, 5, 6)
        pts = fcc_lattice(n, 2.6)
        m = bench_mol("cu", fill("Cu", length(pts)), pts)
        na = length(pts)
        for n_exp in (512, 2048, 8192)
            for tol in (0.0, 0.8, 5.0)
                SASA.sasa(m; probe, n_occ = 512, n_exp, area_tol = tol)   # warm
                t = best(() -> SASA.sasa(m; probe, n_occ = 512, n_exp, area_tol = tol))
                push!(rows, Any[na, n_exp, tol, t*1e3, t*1e6/na])
            end
        end
    end
    ["natoms","n_exp","area_tol","ms","us_per_atom"], rows
end

# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------

_col(rows, j) = [r[j] for r in rows]
_where(rows, j, v) = [r for r in rows if r[j] == v]
_p(v, q) = (s = sort(v); isempty(s) ? NaN : s[clamp(ceil(Int, q*length(s)), 1, length(s))])

"""
    bench_report(; csv_path = nothing) -> Nothing

Print the full accuracy / skip-firing / speed breakdown without opening a
window. Pass `csv_path` to also dump the raw accuracy sweep as CSV.
"""
function bench_report(; csv_path = nothing)
    hdr, rows = accuracy_sweep(csv_path)
    classes = unique(_col(rows, 2))

    println("="^100)
    println("SASA sweep -- ", length(unique(_col(rows,1))), " systems x ",
            length(BENCH_TOLS), " area_tol = ", length(rows), " runs")
    println("  probe ", BENCH_PROBE, " A | n_exp ", BENCH_NEXP,
            " | n_occ 512 | reference n_exp ", BENCH_NREF)
    println("="^100)

    println("\n--- 1. SKIP FIRING + ERROR vs area_tol (all systems pooled) ---")
    @printf("%9s %9s %10s %13s %12s %12s %13s\n",
            "area_tol", "runs", "w/ skips", "atoms skipped", "mean err%",
            "worst err%", "max lost A^2")
    for tol in BENCH_TOLS
        rs = _where(rows, 4, tol)
        nsk = sum(_col(rs, 8)); nat = sum(_col(rs, 3))
        errs = 100 .* _col(rs, 11)
        @printf("%9.2f %9d %10d %12.3f%% %11.4f%% %11.4f%% %13.3f\n",
                tol, length(rs), count(>(0), _col(rs, 8)), 100*nsk/nat,
                mean(errs), minimum(errs), maximum(_col(rs, 9)))
    end

    println("\n--- 2. WORST relative error by class and area_tol (%) ---")
    @printf("%-16s %6s", "class", "n_sys")
    for tol in BENCH_TOLS; @printf("%10.2f", tol); end; println()
    for c in classes
        @printf("%-16s %6d", c, length(unique(_col(_where(rows,2,c), 1))))
        for tol in BENCH_TOLS
            rs = [r for r in rows if r[2] == c && r[4] == tol]
            @printf("%10.4f", isempty(rs) ? 0.0 : 100*minimum(_col(rs, 11)))
        end
        println()
    end

    println("\n--- 3. _classify census by class (% of atoms) ---")
    @printf("%-16s %9s %12s %12s %12s\n", "class", "atoms", "exact-exp", "exact-bur", "sampled")
    for c in classes
        rs = [r for r in rows if r[2] == c && r[4] == 0.0]
        tot = sum(_col(rs, 3))
        @printf("%-16s %9d %11.2f%% %11.2f%% %11.2f%%\n", c, tot,
                100*sum(_col(rs,13))/tot, 100*sum(_col(rs,14))/tot,
                100*sum(_col(rs,15))/tot)
    end

    println("\n--- 4. Sampling error at n_exp=", BENCH_NEXP, " vs ", BENCH_NREF,
            "-pt reference (area_tol=0) ---")
    @printf("%-16s %13s %13s %13s\n", "class", "mean |err|", "p95 |err|", "max |err|")
    for c in classes
        rs = [r for r in rows if r[2] == c && r[4] == 0.0]
        e = 100 .* abs.(_col(rs, 12))
        isempty(e) && continue
        @printf("%-16s %12.4f%% %12.4f%% %12.4f%%\n", c, mean(e), _p(e,0.95), maximum(e))
    end

    ch, crows = convergence_sweep()
    println("\n--- 5. Convergence in n_exp (rel. error vs ", BENCH_NREF, "-pt reference) ---")
    tags = unique(_col(crows, 1))
    @printf("%9s", "n_exp"); for t in tags; @printf("%15s", first(t, 14)); end; println()
    for n in unique(_col(crows, 3))
        @printf("%9d", n)
        for t in tags
            k = findfirst(x -> x[1] == t && x[3] == n, crows)
            @printf("%14.4f%%", k === nothing ? NaN : 100*crows[k][6])
        end
        println()
    end

    ah, arows = analytic_sweep()
    println("\n--- 6. Analytic ground truth: 2 equal spheres, 4pi.rho^2 - 2pi.rho(rho-d/2) ---")
    ae = 100 .* abs.(_col(arows, 7))
    @printf("  %d configs | radii %.1f-%.1f A | d/2rho %.2f-%.2f\n",
            length(arows), minimum(_col(arows,1)), maximum(_col(arows,1)),
            minimum(_col(arows,4)), maximum(_col(arows,4)))
    @printf("  |rel err|  mean %.4f%%  median %.4f%%  p95 %.4f%%  max %.4f%%\n",
            mean(ae), median(ae), _p(ae,0.95), maximum(ae))
    println("  worst 5:")
    for r in sort(arows, by = x -> -abs(x[7]))[1:5]
        @printf("    r=%.1f d/2rho=%.2f  measured %9.3f  analytic %9.3f  err %+.4f%%\n",
                r[1], r[4], r[5], r[6], 100*r[7])
    end

    sh, srows = speed_sweep()
    println("\n--- 7. Speed ---")
    @printf("%8s %9s %10s %11s %12s\n", "atoms", "n_exp", "area_tol", "ms", "us/atom")
    for r in srows
        @printf("%8d %9d %10.1f %11.3f %12.2f\n", r[1], r[2], r[3], r[4], r[5])
    end
    println("\n  scaling in atom count (n_exp 2048, area_tol 0):")
    for r in [x for x in srows if x[2] == 2048 && x[3] == 0.0]
        @printf("    %5d atoms -> %8.3f ms  (%6.2f us/atom)\n", r[1], r[4], r[5])
    end
    big = maximum(_col(srows, 1))
    t0 = only([r[4] for r in srows if r[1] == big && r[2] == 2048 && r[3] == 0.0])
    println("\n  area_tol speedup at ", big, " atoms, n_exp 2048:")
    for tol in (0.0, 0.8, 5.0)
        t = only([r[4] for r in srows if r[1] == big && r[2] == 2048 && r[3] == tol])
        @printf("    area_tol %4.1f -> %8.3f ms  (%.2fx)\n", tol, t, t0/t)
    end
    println("\n", "="^100)
    return nothing
end

# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

const _EXPOSED_C = RGBf(1.00, 0.62, 0.13)
const _OCC_C     = RGBf(0.24, 0.26, 0.32)

"Distinct colours for an arbitrary number of series."
_series_colors(n) = [get(GLMakie.Makie.ColorSchemes.viridis, i / max(n, 1)) for i in 1:n]

"""
    bench_figure(; n_exp, n_occ, probe) -> Figure

Six-panel summary of the accuracy sweep:

1. area lost to the `area_tol` early exit, and how many atoms it skipped
2. worst relative error per system class against `area_tol`
3. what fraction of atoms `_classify` settles exactly vs hands to sampling
4. sampling error against a converged reference, by class
5. convergence in `n_exp` for one system per class
6. error against the analytic two-sphere cap, versus overlap

Every number here is produced by `test/sasa_bench.jl`; this only draws it.
"""
function bench_figure(; n_exp = BENCH_NEXP, n_occ = 512, probe = BENCH_PROBE)
    _, rows = accuracy_sweep(; n_exp, n_occ, probe)
    classes = unique(_col(rows, 2))
    cols = _series_colors(length(classes))

    fig = Figure(size = (1500, 950))
    Label(fig[0, 1:3],
          "SASA accuracy and cost  ($(length(unique(_col(rows,1)))) systems, " *
          "probe $(probe) A, n_exp $(n_exp))"; fontsize = 20, font = :bold)

    # --- 1. skip firing + area lost vs area_tol -----------------------------
    ax1 = Axis(fig[1, 1]; title = "area_tol early exit",
               xlabel = "area_tol (A^2)", ylabel = "atoms skipped (%)")
    ax1b = Axis(fig[1, 1]; ylabel = "total area lost (A^2)",
                yaxisposition = :right, yticklabelcolor = _EXPOSED_C)
    hidespines!(ax1b); hidexdecorations!(ax1b)
    pct  = Float64[]; lost = Float64[]
    for tol in BENCH_TOLS
        rs = _where(rows, 4, tol)
        push!(pct, 100 * sum(_col(rs, 8)) / sum(_col(rs, 3)))
        push!(lost, sum(_col(rs, 9)))
    end
    scatterlines!(ax1, BENCH_TOLS, pct; color = _OCC_C, markersize = 10)
    scatterlines!(ax1b, BENCH_TOLS, lost; color = _EXPOSED_C, markersize = 10)

    # --- 2. worst relative error per class ---------------------------------
    ax2 = Axis(fig[1, 2]; title = "worst relative error by class",
               xlabel = "area_tol (A^2)", ylabel = "worst error (%)")
    for (k, c) in enumerate(classes)
        ys = [(rs = [r for r in rows if r[2] == c && r[4] == tol];
               isempty(rs) ? 0.0 : 100 * minimum(_col(rs, 11))) for tol in BENCH_TOLS]
        scatterlines!(ax2, BENCH_TOLS, ys; color = cols[k], label = c, markersize = 8)
    end
    axislegend(ax2; position = :lb, labelsize = 9, framevisible = false)

    # --- 3. _classify census ------------------------------------------------
    ax3 = Axis(fig[1, 3]; title = "how atoms are resolved",
               ylabel = "% of atoms", xticks = (1:length(classes), classes),
               xticklabelrotation = pi/5, xticklabelsize = 9)
    exact_e = Float64[]; exact_b = Float64[]; samp = Float64[]
    for c in classes
        rs = [r for r in rows if r[2] == c && r[4] == 0.0]
        t = sum(_col(rs, 3))
        push!(exact_e, 100*sum(_col(rs,13))/t)
        push!(exact_b, 100*sum(_col(rs,14))/t)
        push!(samp,    100*sum(_col(rs,15))/t)
    end
    n = length(classes)
    barplot!(ax3, repeat(1:n, 3), vcat(exact_e, exact_b, samp);
             stack = repeat(1:3, inner = n),
             color = repeat([_EXPOSED_C, _OCC_C, RGBf(0.35,0.55,0.75)], inner = n))
    Legend(fig[2, 3], [PolyElement(color = c) for c in
                       (_EXPOSED_C, _OCC_C, RGBf(0.35,0.55,0.75))],
           ["exact: exposed", "exact: buried", "sampled"];
           orientation = :horizontal, framevisible = false, labelsize = 9)

    # --- 4. sampling error by class ----------------------------------------
    ax4 = Axis(fig[3, 1]; title = "sampling error vs $(BENCH_NREF)-pt reference",
               ylabel = "|relative error| (%)", yscale = log10,
               xticks = (1:length(classes), classes),
               xticklabelrotation = pi/5, xticklabelsize = 9)
    for (k, c) in enumerate(classes)
        rs = [r for r in rows if r[2] == c && r[4] == 0.0]
        e = filter(>(0), 100 .* abs.(_col(rs, 12)))
        isempty(e) && continue
        scatter!(ax4, fill(k, length(e)), e; color = (cols[k], 0.45), markersize = 6)
        scatter!(ax4, [k], [mean(e)]; color = :black, marker = :hline, markersize = 22)
    end

    # --- 5. convergence in n_exp -------------------------------------------
    _, crows = convergence_sweep(; probe)
    ax5 = Axis(fig[3, 2]; title = "convergence in n_exp",
               xlabel = "n_exp", ylabel = "|relative error| (%)",
               xscale = log2, yscale = log10)
    for (k, t) in enumerate(unique(_col(crows, 1)))
        rs = [r for r in crows if r[1] == t]
        ys = max.(1e-4, 100 .* abs.(_col(rs, 6)))
        scatterlines!(ax5, Float64.(_col(rs, 3)), ys;
                      color = cols[mod1(k, length(cols))], label = t, markersize = 8)
    end
    axislegend(ax5; position = :lb, labelsize = 9, framevisible = false)

    # --- 6. analytic ground truth ------------------------------------------
    _, arows = analytic_sweep(; probe)
    # every radius traces the same curve: the exposed fraction depends only on
    # d/2rho, so the point set meets the cap boundary identically at any scale.
    ax6 = Axis(fig[3, 3]; title = "error vs analytic cap (all radii coincide)",
               xlabel = "d / 2rho  (0 = coincident, 1 = tangent)",
               ylabel = "relative error (%)")
    hlines!(ax6, [0.0]; color = (:black, 0.3))
    rr = unique(_col(arows, 1))
    for (k, r) in enumerate(rr)
        rs = [x for x in arows if x[1] == r]
        scatterlines!(ax6, Float64.(_col(rs, 4)), 100 .* _col(rs, 7);
                      color = _series_colors(length(rr))[k],
                      label = "r = $(r) A", markersize = 7)
    end
    axislegend(ax6; position = :lb, labelsize = 8, framevisible = false, nbanks = 2)

    fig
end

"""
    speed_figure() -> Figure

Wall time against atom count and `n_exp`, plus the per-atom cost, showing
where `area_tol` actually buys anything.
"""
function speed_figure()
    _, srows = speed_sweep()
    fig = Figure(size = (1150, 430))
    Label(fig[0, 1:2], "SASA cost"; fontsize = 20, font = :bold)

    ax1 = Axis(fig[1, 1]; title = "wall time vs atom count",
               xlabel = "atoms", ylabel = "ms", xscale = log10, yscale = log10)
    nexps = unique(_col(srows, 2))
    cs = _series_colors(length(nexps))
    for (k, ne) in enumerate(nexps)
        rs = [r for r in srows if r[2] == ne && r[3] == 0.0]
        scatterlines!(ax1, Float64.(_col(rs, 1)), _col(rs, 4);
                      color = cs[k], label = "n_exp $(ne)", markersize = 9)
    end
    axislegend(ax1; position = :lt, labelsize = 9, framevisible = false)

    ax2 = Axis(fig[1, 2]; title = "per-atom cost, by area_tol",
               xlabel = "atoms", ylabel = "us / atom")
    tols = unique(_col(srows, 3))
    ct = _series_colors(length(tols))
    for (k, tol) in enumerate(tols)
        rs = [r for r in srows if r[3] == tol && r[2] == 2048]
        scatterlines!(ax2, Float64.(_col(rs, 1)), _col(rs, 5);
                      color = ct[k], label = "area_tol $(tol)", markersize = 9)
    end
    axislegend(ax2; position = :lt, labelsize = 9, framevisible = false)
    fig
end

"""
    vis_bench(; kwargs...) -> Nothing

Build both figures and display them, blocking until the windows are closed.
"""
function vis_bench(; kwargs...)
    f1 = bench_figure(; kwargs...)
    f2 = speed_figure()
    display(f2)
    wait(display(f1))
    return nothing
end

"""
    save_bench_figures(prefix = "bench"; kwargs...) -> Vector{String}

Render both figures to `<prefix>_accuracy.png` and `<prefix>_speed.png`
without needing a window. Returns the paths written.
"""
function save_bench_figures(prefix::AbstractString = "bench"; kwargs...)
    paths = String[]
    for (suffix, f) in (("accuracy", bench_figure(; kwargs...)), ("speed", speed_figure()))
        p = "$(prefix)_$(suffix).png"
        save(p, f); push!(paths, p)
    end
    paths
end
