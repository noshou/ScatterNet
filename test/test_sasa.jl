using ScatterNet: Interfaces
using ScatterNet.Molecule.SASA: SASA, sasa, sasa_atoms, _occluded
using .Molecules: create

# ---------------------------------------------------------------------------
# Injected radii source.
#
# Every molecule below is built with a fixed lookup table instead of the
# bundled SQLite radii DB, so the geometry under test is exactly what the test
# says it is and nothing here depends on `data/atomic_radii.sqlite3`.
#
# NOTE: `include` evaluates at module top level regardless of the enclosing
# `@testset`, so this `struct` is legal here.
# ---------------------------------------------------------------------------
struct SasaTestRadii <: Interfaces.RadiiSource
    table::Dict{String,Float64}
end

function Interfaces.lookup(s::SasaTestRadii, ions::AbstractVector{<:AbstractString})
    out = Vector{Tuple{String,Union{Float64,Nothing}}}(undef, length(ions))
    for i in eachindex(ions)
        k = String(ions[i])
        out[i] = (k, get(s.table, k, nothing))
    end
    return out
end

# element letter -> radius (Å). Deliberately not real elements: these are shapes.
const SASA_SRC = SasaTestRadii(Dict(
    "q" => 1.5,   # workhorse
    "a" => 1.0,
    "b" => 2.0,
    "c" => 0.5,   # tiny
    "d" => 5.0,   # huge
))

sasa_mol(elms, crds) = create("sasa-test", elms, crds; radii_source = SASA_SRC)

"Analytic area of a lone expanded sphere: 4π(r + probe)²."
sasa_full(r, probe) = 4π * (r + probe)^2

"""
Exposed area of one of two EQUAL spheres of expanded radius ρ whose centers are
`d` apart (0 < d < 2ρ): the occluded part is a spherical cap of area
`2πρ(ρ - d/2)`, so the exposed part is `4πρ² - 2πρ(ρ - d/2)`.
"""
sasa_cap_exposed(ρ, d) = 4π * ρ^2 - 2π * ρ * (ρ - d / 2)

"Cubic lattice of `n`×`n`×`n` points with spacing `a`, in (i, j, k) row-major order."
function sasa_cube(n::Int, a::Float64)
    out = NTuple{3,Float64}[]
    for i in 0:n-1, j in 0:n-1, k in 0:n-1
        push!(out, (a * i, a * j, a * k))
    end
    return out
end

@testset "SASA" begin

    @testset "argument contract" begin
        m = sasa_mol(["q", "q"], [(0.0, 0.0, 0.0), (4.0, 0.0, 0.0)])

        @test_throws DomainError sasa_atoms(m; n_occ = 0, n_exp = 100, probe = 1.4)
        @test_throws DomainError sasa_atoms(m; n_occ = -5, n_exp = 100, probe = 1.4)
        @test_throws DomainError sasa_atoms(m; n_occ = 10, n_exp = 0, probe = 1.4)
        @test_throws DomainError sasa_atoms(m; n_occ = 10, n_exp = -1, probe = 1.4)
        @test_throws DomainError sasa_atoms(m; n_occ = 100, n_exp = 10, probe = 1.4)      # n_exp < n_occ
        @test_throws DomainError sasa_atoms(m; n_occ = 10, n_exp = 100, probe = -1e-9)    # negative probe

        # same checks reach through the `sasa` wrapper
        @test_throws DomainError sasa(m; n_occ = 0, n_exp = 100, probe = 1.4)
        @test_throws DomainError sasa(m; n_occ = 100, n_exp = 10, probe = 1.4)
        @test_throws DomainError sasa(m; n_occ = 10, n_exp = 100, probe = -1.0)

        # boundary cases that MUST be legal
        @test sasa(m; n_occ = 100, n_exp = 100, probe = 1.4) > 0.0        # n_occ == n_exp
        @test sasa(m; n_occ = 1, n_exp = 1, probe = 1.4) > 0.0            # smallest legal counts
        @test sasa(m; n_occ = 10, n_exp = 100, probe = 0.0) > 0.0         # probe == 0 is legal
    end

    @testset "single isolated atom is analytically exact" begin
        # A lone atom has no candidate other than itself, and `_occluded` skips
        # self, so EVERY sample point is exposed: p_exp/n_exp == 1 exactly and
        # the area is 4π(r+probe)² to the last bit, for any point counts.
        for (el, r) in (("c", 0.5), ("a", 1.0), ("q", 1.5), ("d", 5.0))
            m = sasa_mol([el], [(3.0, -2.0, 7.0)])   # centering puts it at the origin
            for probe in (0.0, 1.4, 2.5)
                for (n_occ, n_exp) in ((1, 1), (1, 100), (10, 100), (77, 1000), (250, 250))
                    @test sasa(m; n_occ = n_occ, n_exp = n_exp, probe = probe) == sasa_full(r, probe)
                end
            end
        end
    end

    @testset "two well-separated atoms: exact sum of two spheres" begin
        probe = 1.4
        m = sasa_mol(["a", "b"], [(0.0, 0.0, 0.0), (50.0, 0.0, 0.0)])
        areas = sasa_atoms(m; n_occ = 20, n_exp = 1000, probe = probe)
        @test areas[1] == sasa_full(1.0, probe)
        @test areas[2] == sasa_full(2.0, probe)
        @test sasa(m; n_occ = 20, n_exp = 1000, probe = probe) == sasa_full(1.0, probe) + sasa_full(2.0, probe)

        # exactly touching-but-not-overlapping expanded spheres:
        # ρ₁ + ρ₂ = 2.4 + 3.4 = 5.8, centers 6.0 apart -> still fully exposed.
        m2 = sasa_mol(["a", "b"], [(0.0, 0.0, 0.0), (6.0, 0.0, 0.0)])
        @test sasa(m2; n_occ = 20, n_exp = 1000, probe = probe) == sasa_full(1.0, probe) + sasa_full(2.0, probe)
    end

    @testset "complete engulfment: inner atom is exactly zero" begin
        # r_small = 0.5, r_large = 5.0, probe = 1.0, centers 1.0 apart:
        #   d + r_small + probe = 1.0 + 1.5 = 2.5  <<  r_large + probe = 6.0
        # so the small atom's whole expanded sphere is strictly inside the
        # large one (margin 3.5 Å), and no point of the large sphere can be
        # reached by the small one (6.0 - 2.5 = 3.5 Å of clearance).
        probe = 1.0
        m = sasa_mol(["d", "c"], [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
        for (n_occ, n_exp) in ((1, 100), (20, 1000), (50, 5000))
            areas = sasa_atoms(m; n_occ = n_occ, n_exp = n_exp, probe = probe)
            @test areas[2] == 0.0                       # engulfed: EXACTLY zero
            @test areas[1] == sasa_full(5.0, probe)     # host: full analytic area
            @test sasa(m; n_occ = n_occ, n_exp = n_exp, probe = probe) == sasa_full(5.0, probe)
        end

        # concentric variant (d = 0) with strictly different radii is the same
        # clean case: ρ_small = 2.4 < ρ_large = 3.4.
        mc = sasa_mol(["a", "b"], [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)])
        ac = sasa_atoms(mc; n_occ = 50, n_exp = 1000, probe = 1.4)
        @test ac[1] == 0.0
        @test ac[2] == sasa_full(2.0, 1.4)
    end

    @testset "coincident identical atoms are exactly zero" begin
        # Two atoms with the SAME center and the SAME radius. Every sample point
        # of atom i sits at distance exactly ρ from atom j's center, and the
        # sampled test `dst² <= ρ_c²` is exactly on its boundary there -- so
        # when this went through point sampling it leaked 3-15% of the points
        # either way depending on rounding in `ρ*u` and the sum of squares.
        #
        # `_classify` now decides it before any point is generated: d = 0 and
        # ρᵢ = ρⱼ satisfies `d + ρᵢ <= ρⱼ`, so each atom is engulfed by the
        # other and both are exactly 0. The knife edge is gone, not hidden, so
        # this asserts the exact value at every point count.
        probe = 1.4
        for n_exp in (1, 100, 500, 2000, 10000)
            areas = sasa_atoms(sasa_mol(
                ["q", "q"], [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]); n_occ = min(50, n_exp), n_exp = n_exp, probe = probe)
            @test areas == [0.0, 0.0]
        end
    end

    @testset "partial occlusion: monotone in separation" begin
        probe = 1.4
        ρ = 1.5 + probe
        iso = 2 * sasa_full(1.5, probe)
        ds = (10.0, 8.0, 2 * ρ, 5.0, 4.0, 3.0, 2.0, 1.0, 0.5)

        prev_tot = Inf
        prev = (Inf, Inf)
        for d in ds
            m = sasa_mol(["q", "q"], [(0.0, 0.0, 0.0), (d, 0.0, 0.0)])
            a = sasa_atoms(m; n_occ = 50, n_exp = 4000, probe = probe)
            tot = sum(a)

            @test tot <= prev_tot + 1e-9            # total never increases
            @test a[1] <= prev[1] + 1e-9            # nor does either atom
            @test a[2] <= prev[2] + 1e-9

            if d >= 2 * ρ
                @test tot == iso                    # no overlap -> exact
            else
                @test tot < iso                     # overlap -> strictly less
                @test a[1] < sasa_full(1.5, probe)
                @test a[2] < sasa_full(1.5, probe)
            end
            prev_tot = tot
            prev = (a[1], a[2])
        end
    end

    @testset "analytic spherical cap (strongest correctness check)" begin
        # Two EQUAL spheres, expanded radius ρ = 2.9, centers d apart with
        # 0 < d < 2ρ. Each atom's exposed area is 4πρ² - 2πρ(ρ - d/2) exactly.
        #
        # Tolerance: relative, empirically tuned. With n_exp = 40_000 the
        # measured worst relative error over these eight separations was
        # 3.12e-4 (typical ~1e-4). 1% therefore leaves ~30x headroom, which is
        # what a quasi-random (non-i.i.d., so not variance-bounded) point set
        # warrants -- tight enough that a real bias of even a percent fails.
        probe = 1.4
        ρ = 1.5 + probe
        rtol = 0.01
        for d in (0.5, 1.0, 1.7, 2.5, 3.3, 4.2, 5.0, 5.7)
            @test 0.0 < d < 2ρ
            exact = sasa_cap_exposed(ρ, d)
            m = sasa_mol(["q", "q"], [(0.0, 0.0, 0.0), (d, 0.0, 0.0)])
            a = sasa_atoms(m; n_occ = 64, n_exp = 40_000, probe = probe)
            @test abs(a[1] - exact) < rtol * exact
            @test abs(a[2] - exact) < rtol * exact
            @test abs(sum(a) - 2 * exact) < rtol * 2 * exact
        end
    end

    @testset "convergence: more points -> closer to the analytic cap" begin
        # Quasi-random error is not monotone point-by-point (n_exp = 500 can
        # beat n_exp = 1000), so this compares a genuinely coarse count against
        # a genuinely fine one, averaged over several separations. Measured
        # mean relative error over these seven separations: 2.96e-2 at
        # n_exp = 32, 2.13e-2 at 64, 1.16e-2 at 128, 1.39e-3 at 1000,
        # 2.38e-4 at 16_000. The coarse/fine gap is ~90x, so the 10x margin
        # asserted below is not knife-edge.
        probe = 1.4
        ρ = 1.5 + probe
        ds = (0.6, 1.4, 2.2, 3.0, 3.8, 4.6, 5.4)

        mean_rel_err(n_exp) = sum(ds) do d
            exact = sasa_cap_exposed(ρ, d)
            m = sasa_mol(["q", "q"], [(0.0, 0.0, 0.0), (d, 0.0, 0.0)])
            a = sasa_atoms(m; n_occ = min(32, n_exp), n_exp = n_exp, probe = probe)
            (abs(a[1] - exact) + abs(a[2] - exact)) / (2 * exact)
        end / length(ds)

        coarse = mean_rel_err(64)
        fine   = mean_rel_err(16_000)
        @test fine < coarse / 10
        @test fine < 1e-3
    end

    @testset "n_occ invariance (coarse pass is a prefix of the fine pass)" begin
        # `n_occ` only gates the buried early-exit; its points are the first
        # `n_occ` terms of the SAME plastic sequence the fine pass reuses. So
        # for a molecule with no buried atom, varying n_occ at fixed n_exp must
        # give the IDENTICAL answer, bit for bit.
        probe = 1.4

        m3 = sasa_mol(["q", "q", "q"], [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (1.5, 2.5, 0.0)])
        ref3 = sasa_atoms(m3; n_occ = 50, n_exp = 2000, probe = probe)
        for n_occ in (2, 5, 20, 100, 500, 2000)
            @test sasa_atoms(m3; n_occ = n_occ, n_exp = 2000, probe = probe) == ref3
        end

        m8 = sasa_mol(fill("q", 8), sasa_cube(2, 4.0))
        ref8 = sasa_atoms(m8; n_occ = 64, n_exp = 1000, probe = probe)
        for n_occ in (8, 32, 64, 256, 1000)
            @test sasa_atoms(m8; n_occ = n_occ, n_exp = 1000, probe = probe) == ref8
        end
        @test all(>(0.0), ref8)   # nothing is buried in this configuration
    end

    @testset "regression: a small n_occ can no longer bury an exposed atom" begin
        # A finite sample can PROVE exposure (one unoccluded point is a witness)
        # but can never prove burial. This trimer's third atom is genuinely
        # exposed, and at n_occ = 1 its single coarse point happens to land
        # occluded -- which used to zero the atom outright. The rule-of-three
        # demotion means an unwitnessed atom is now confirmed against the full
        # n_exp set instead, so the answer no longer depends on n_occ at all.
        probe = 1.4
        m = sasa_mol(["q", "q", "q"], [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (1.5, 2.5, 0.0)])
        ref = sasa_atoms(m; n_occ = 64, n_exp = 2000, probe = probe)
        @test ref[3] > 0.0
        for n_occ in (1, 2, 3, 8, 64)
            @test sasa_atoms(m; n_occ = n_occ, n_exp = 2000, probe = probe) == ref
        end

        # raising area_tol re-arms the early exit: with an absurd tolerance the
        # coarse pass is allowed to call an unwitnessed atom buried again.
        lax = sasa_atoms(m; n_occ = 1, n_exp = 2000, probe = probe, area_tol = 1e9)
        @test lax[3] == 0.0
        @test lax[1] == ref[1] && lax[2] == ref[2]   # witnessed atoms unaffected
        @test_throws DomainError sasa_atoms(m; n_occ = 64, n_exp = 2000, probe = probe, area_tol = -1.0)
    end

    @testset "defaults" begin
        m = sasa_mol(["q", "q"], [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0)])

        # a bare call == spelling every default out explicitly
        ref = sasa_atoms(m; probe = 1.4, n_occ = 512, n_exp = 4096, area_tol = 0.8)
        @test sasa_atoms(m) == ref
        @test sasa(m) == sum(ref)
        @test sasa(m) isa Float64
        @test sasa_atoms(m) isa Vector{Float64}
        @test length(sasa_atoms(m)) == 2

        # each keyword is independently overridable
        @test sasa_atoms(m; probe = 0.0) != ref
        @test sasa_atoms(m; n_exp = 512) != ref

        # a lone atom is still exact through the default path
        @test sasa(sasa_mol(["q"], [(0.0, 0.0, 0.0)])) == sasa_full(1.5, 1.4)
    end

    @testset "the default n_exp is accurate enough to not be the limiting error" begin
        # The default must land the sampled answer within ~0.1% of the analytic
        # cap -- an order of magnitude under the ~1.3% a 0.05 A radius shift
        # already costs, so refining the mesh further chases noise the radius
        # table cannot justify.
        probe = 1.4
        ρ = 1.5 + probe
        worst = 0.0
        for d in (1.0, 2.0, 3.0, 3.5, 4.0, 5.0)
            m = sasa_mol(["q", "q"], [(0.0, 0.0, 0.0), (d, 0.0, 0.0)])
            got = sasa_atoms(m)[1]
            worst = max(worst, abs(got - sasa_cap_exposed(ρ, d)) / sasa_cap_exposed(ρ, d))
        end
        @test worst < 0.005          # measured ~1e-3; generous headroom
    end

    @testset "exact pre-filter: the two decidable regimes need no sampling" begin
        # `_classify` settles these from the neighbour list alone, so the answer
        # is exact at ANY point count -- n_exp = 1 is enough. Sampling could not
        # produce these numbers reliably.
        probe = 1.4

        # no neighbour reaches the surface -> exactly 4πρ², fully exposed
        far = sasa_mol(["q", "q"], [(0.0, 0.0, 0.0), (10.0, 0.0, 0.0)])
        ρ = Molecules.radii(far)[1] + probe
        @test sasa_atoms(far; n_occ = 1, n_exp = 1, probe = probe) == [4π * ρ^2, 4π * ρ^2]

        # exact tangency (d == ρᵢ + ρⱼ) cuts a measure-zero cap: still full
        tan_ = sasa_mol(["q", "q"], [(0.0, 0.0, 0.0), (2ρ, 0.0, 0.0)])
        @test sasa_atoms(tan_; n_occ = 1, n_exp = 1, probe = probe) == [4π * ρ^2, 4π * ρ^2]

        # one neighbour engulfs the other -> inner is exactly 0.0, outer exactly full
        eng = sasa_mol(["c", "d"], [(3.0, 0.0, 0.0), (0.0, 0.0, 0.0)])
        ρs = Molecules.radii(eng)[1] + probe; ρb = Molecules.radii(eng)[2] + probe
        @test 3.0 + ρs <= ρb                       # geometry really is engulfment
        for n_exp in (1, 16, 4096)
            a = sasa_atoms(eng; n_occ = 1, n_exp = n_exp, probe = probe)
            @test a[1] == 0.0
            @test a[2] == 4π * ρb^2
        end

    end

    @testset "regression: sasa must not be identically zero (self-occlusion)" begin
        # Every sample point of atom i lies at exactly rads[i] + probe from atom
        # i's own center, so an occlusion test that did not skip `self` would
        # report `dst² <= ρ_self²` for every point of every atom and `sasa`
        # would return 0.0 for EVERY molecule. `_occluded`'s `self` argument is
        # the guard; these assertions fail loudly if it is ever dropped.
        probe = 1.4
        for m in (
            sasa_mol(["q"], [(0.0, 0.0, 0.0)]),
            sasa_mol(["a", "b"], [(0.0, 0.0, 0.0), (2.5, 0.0, 0.0)]),
            sasa_mol(["q", "q", "q"], [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (1.0, 1.7, 0.0)]),
            sasa_mol(fill("q", 27), sasa_cube(3, 2.0)),
        )
            @test sasa(m; n_occ = 50, n_exp = 1000, probe = probe) > 0.0
            @test any(>(0.0), sasa_atoms(m; n_occ = 50, n_exp = 1000, probe = probe))
        end
    end

    @testset "sasa == sum(sasa_atoms)" begin
        probe = 1.4
        configs = (
            sasa_mol(["a"], [(0.0, 0.0, 0.0)]),
            sasa_mol(["a", "b"], [(0.0, 0.0, 0.0), (30.0, 0.0, 0.0)]),
            sasa_mol(["q", "q"], [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)]),
            sasa_mol(["d", "c"], [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]),
            sasa_mol(["q", "a", "b", "c"], [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 3.0, 0.0), (1.0, 1.0, 1.0)]),
            sasa_mol(fill("q", 27), sasa_cube(3, 2.0)),
        )
        for m in configs, (n_occ, n_exp) in ((1, 200), (32, 1500), (100, 100))
            @test sasa(m; n_occ = n_occ, n_exp = n_exp, probe = probe) === sum(sasa_atoms(m; n_occ = n_occ, n_exp = n_exp, probe = probe))
        end
    end

    @testset "sasa_atoms: shape, type, and per-atom bounds" begin
        probe = 1.4
        elms = ["q", "a", "b", "c", "d"]
        crds = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 3.5, 0.0),
                (1.0, 1.0, 1.0), (-4.0, 2.0, 1.0)]
        rs   = [1.5, 1.0, 2.0, 0.5, 5.0]
        m = sasa_mol(elms, crds)
        areas = sasa_atoms(m; n_occ = 40, n_exp = 1200, probe = probe)

        @test areas isa Vector{Float64}
        @test length(areas) == length(elms)
        @test length(areas) == size(Molecules.coords_cartesian(m), 2)
        for (i, r) in enumerate(rs)
            @test areas[i] >= 0.0
            @test areas[i] <= sasa_full(r, probe)   # can never exceed a full sphere
        end
        @test sasa(m; n_occ = 40, n_exp = 1200, probe = probe) >= 0.0
    end

    @testset "determinism (quasi-random, not random)" begin
        probe = 1.4
        m = sasa_mol(["q", "a", "b", "c"], [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 3.0, 0.0), (1.0, 1.0, 1.0)])
        a = sasa_atoms(m; n_occ = 32, n_exp = 1500, probe = probe)
        for _ in 1:3
            @test sasa_atoms(m; n_occ = 32, n_exp = 1500, probe = probe) == a       # bit-identical
        end
        @test sasa(m; n_occ = 32, n_exp = 1500, probe = probe) === sasa(m; n_occ = 32, n_exp = 1500, probe = probe)

        # a freshly-built but geometrically identical molecule agrees too
        m2 = sasa_mol(["q", "a", "b", "c"], [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 3.0, 0.0), (1.0, 1.0, 1.0)])
        @test sasa_atoms(m2; n_occ = 32, n_exp = 1500, probe = probe) == a
    end

    @testset "probe scaling on a lone atom is exactly (r + probe)²" begin
        m = sasa_mol(["a"], [(1.0, 2.0, 3.0)])         # r = 1.0
        v0 = sasa(m; n_occ = 10, n_exp = 100, probe = 0.0)
        v1 = sasa(m; n_occ = 10, n_exp = 100, probe = 1.5)
        v2 = sasa(m; n_occ = 10, n_exp = 100, probe = 4.0)
        @test v0 == sasa_full(1.0, 0.0)
        @test v1 / v0 == (2.5 / 1.0)^2
        @test v2 / v0 == (5.0 / 1.0)^2
        @test check_float(v2 / v1, (5.0 / 2.5)^2)
    end

    @testset "_occluded unit tests" begin
        # atom 1 at origin (r = 1.0), atom 2 at x = 3.0 (r = 1.0), probe = 0.5
        # -> expanded radii 1.5 each.
        crds = [0.0 3.0; 0.0 0.0; 0.0 0.0]
        rads = [1.0, 1.0]
        probe = 0.5

        # strictly inside candidate 2's expanded sphere
        @test _occluded((3.2, 0.0, 0.0), [1, 2], crds, rads, probe, 1)
        @test _occluded((3.0, 0.4, -0.3), [2], crds, rads, probe, 1)

        # strictly outside every candidate's expanded sphere
        @test !_occluded((10.0, 0.0, 0.0), [1, 2], crds, rads, probe, 1)
        @test !_occluded((1.5 + 1e-9, 0.0, 0.0), [1, 2], crds, rads, probe, 2)

        # `self` is skipped even for a point ON self's own expanded sphere ...
        @test !_occluded((1.5, 0.0, 0.0), [1], crds, rads, probe, 1)
        @test !_occluded((0.0, 0.0, 1.5), [1], crds, rads, probe, 1)
        # ... and even for a point deep INSIDE self (self's own center)
        @test !_occluded((0.0, 0.0, 0.0), [1], crds, rads, probe, 1)
        @test !_occluded((0.0, 0.0, 0.0), [1, 2], crds, rads, probe, 1)
        # the same point IS occluded once atom 1 is not `self`
        @test _occluded((0.0, 0.0, 0.0), [1, 2], crds, rads, probe, 2)

        # degenerate candidate lists
        @test !_occluded((0.0, 0.0, 0.0), Int[], crds, rads, probe, 1)
        @test !_occluded((3.0, 0.0, 0.0), Int[], crds, rads, probe, 1)   # inside 2, but not a candidate
        @test !_occluded((3.0, 0.0, 0.0), [2], crds, rads, probe, 2)     # only self

        # boundary: `dst² <= ρ_c²` is inclusive, so exactly on the surface counts
        @test _occluded((1.5, 0.0, 0.0), [1], crds, rads, probe, 2)

        # probe widens the occluding sphere
        @test !_occluded((2.0, 0.0, 0.0), [1], crds, rads, 0.5, 2)   # ρ = 1.5 < 2.0
        @test _occluded((2.0, 0.0, 0.0), [1], crds, rads, 1.5, 2)    # ρ = 2.5 > 2.0

        @test _occluded((3.2, 0.0, 0.0), [1, 2], crds, rads, probe, 1) isa Bool
    end

    @testset "larger molecule: 3x3x3 cubic lattice" begin
        # spacing a = 2.0, r = 1.5, probe = 1.4 -> ρ = 2.9. A point of the
        # centre atom in direction u is occluded by the axial neighbour ê iff
        # u·ê >= a/(2ρ) = 0.345; over the six axial neighbours the worst-case
        # direction (1,1,1)/√3 still achieves 0.577 > 0.345, so the centre atom
        # is provably occluded in EVERY direction -> exactly 0 for any n_exp.
        probe = 1.4
        r = 1.5
        crds = sasa_cube(3, 2.0)
        m = sasa_mol(fill("q", 27), crds)
        areas = sasa_atoms(m; n_occ = 50, n_exp = 1000, probe = probe)
        total = sum(areas)
        iso = 27 * sasa_full(r, probe)

        @test length(areas) == 27
        @test all(>=(0.0), areas)
        @test all(a -> a <= sasa_full(r, probe), areas)
        @test 0.0 < total < iso
        @test total == sasa(m; n_occ = 50, n_exp = 1000, probe = probe)

        # index of (i, j, k) in sasa_cube(3, ·): k varies fastest
        idx(i, j, k) = 9i + 3j + k + 1
        centre  = idx(1, 1, 1)
        corners = [idx(i, j, k) for i in (0, 2), j in (0, 2), k in (0, 2)]
        faces   = [idx(1, 1, 0), idx(1, 1, 2), idx(1, 0, 1), idx(1, 2, 1), idx(0, 1, 1), idx(2, 1, 1)]

        # the fully engulfed interior atom is EXACTLY zero
        @test areas[centre] == 0.0
        # this must be robust to the point count
        for n_exp in (200, 1000, 4000)
            @test sasa_atoms(m; n_occ = min(50, n_exp), n_exp = n_exp, probe = probe)[centre] == 0.0
        end

        # corners are the most exposed atoms in the lattice
        @test all(>(0.0), areas[corners])
        for c in corners, f in faces
            @test areas[f] < areas[c]
        end
        @test maximum(areas) == maximum(areas[corners])

        # every non-corner atom is strictly less exposed than every corner
        interior = setdiff(1:27, corners)
        @test maximum(areas[interior]) < minimum(areas[corners])

        # a corner atom (3 neighbours) keeps a good fraction of its sphere;
        # the lattice as a whole loses most of its isolated area
        @test 0.05 < total / iso < 0.5

        # loosening the lattice raises the total; tightening it lowers it
        loose = sasa(sasa_mol(fill("q", 27), sasa_cube(3, 3.0)); n_occ = 50, n_exp = 1000, probe = probe)
        tight = sasa(sasa_mol(fill("q", 27), sasa_cube(3, 1.5)); n_occ = 50, n_exp = 1000, probe = probe)
        @test tight < total < loose < iso
    end

    @testset "PlasticMap is reachable through SASA" begin
        # SASA's sampling is the plastic sequence; the submodule is re-exported
        # from here and its own behaviour is covered in test_plasticmap.jl.
        @test SASA.PlasticMap === PlasticMap
        @test length(SASA.PlasticMap.plastic_points(8)) == 8
    end
end
