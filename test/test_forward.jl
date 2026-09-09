# Exercises src/Scattering/Forward.jl -- the module-level forward model that
# composes `vacuo`/`excluded`/`hydration` -> `gram` -> `intensity` ->
# `intensity_calc` into `I_calc(q)`, and the `Scattering`-level configuration
# constants those defaults come from.
#
# The reference side re-does the composition by hand from the per-species
# primitives, so a mismatch is a wiring bug in Forward.jl, not a physics bug
# (the physics is checked in test_partialwave.jl / test_scatterers.jl /
# test_intensity.jl).

using .Scattering: forward, gram_matrix, species_multipoles,
                   forward_cache, ForwardCache, mean_atomic_radius,
                   excluded_volume_factor, contrast_matrix,
                   gram, intensity, intensity_calc, contrast_vector,
                   partial_wave_weights, vacuo, excluded, hydration,
                   SHELL_THICKNESS, PROBE_RADIUS, SHELL_N_TARGET, SHELL_CLASSES,
                   DRO_UNIT, FORM_FACTOR_SOURCE, B_LM_CHUNK
using .Molecules: create, elms
using ScatterNet.Molecule: SASA
using .Molecules: radii
using LinearAlgebra: issymmetric, eigvals

fwd_mol() = create("gly", ["n", "c", "c", "o", "o", "h", "h", "h"],
    [(-1.9, 0.2, 0.1), (-0.5, -0.3, 0.0), (0.6, 0.7, -0.1),
     ( 1.8, 0.2, 0.0), (0.4, 1.9, -0.2), (-2.6, -0.5, 0.0),
     (-0.4, -1.0, 0.8), (0.7, 1.3, 0.8)])
fwd_q      = [0.0, 0.03, 0.07, 0.15, 0.31]
fwd_E      = 9000.0
fwd_lmax   = 4
fwd_chunk  = UInt64(3)

@testset "Forward" begin

    @testset "module-level config constants" begin
        @test SHELL_THICKNESS === 3.0
        @test PROBE_RADIUS    === 1.4
        @test SHELL_N_TARGET  === nothing
        @test SHELL_CLASSES   === (SASA.CONVEX, SASA.CONCAVE, SASA.CAVITY)
        @test DRO_UNIT        === 0.03
        @test B_LM_CHUNK isa Unsigned
        @test FORM_FACTOR_SOURCE isa ScatterNet.Interfaces.FormFactorSource
        # the primitives really do read these as their defaults
        m = fwd_mol()
        @test hydration(m, fwd_q, 2, fwd_chunk) ==
              hydration(m, fwd_q, 2, fwd_chunk; thickness = SHELL_THICKNESS,
                        probe = PROBE_RADIUS, n_target = SHELL_N_TARGET,
                        classes = SHELL_CLASSES)
    end

    @testset "species_multipoles: 5 species in canonical order" begin
        m  = fwd_mol()
        Bs = species_multipoles(m, fwd_q, fwd_lmax, fwd_E; chunk = fwd_chunk)
        @test length(Bs) == 5
        K = (fwd_lmax + 1) * (fwd_lmax + 2) ÷ 2
        @test all(B -> size(B, 2) == K && size(B, 3) == length(fwd_q), Bs)
        @test size(Bs[1], 1) == 2                     # vac: anomalous channels at 9 keV
        @test all(B -> size(B, 1) == 1, Bs[2:5])      # dummies: one real channel

        # element-by-element identical to calling the primitives directly
        ref_vac = vacuo(m, fwd_q, fwd_lmax, elms(m), fwd_E, fwd_chunk)
        ref_ex  = excluded(m, fwd_q, fwd_lmax, fwd_chunk)
        ref_sh  = hydration(m, fwd_q, fwd_lmax, fwd_chunk)
        @test Bs[1] == ref_vac
        @test Bs[2] == ref_ex
        @test Bs[3] == ref_sh.convex
        @test Bs[4] == ref_sh.concave
        @test Bs[5] == ref_sh.cavity
    end

    @testset "gram_matrix: (5,5,Q), symmetric, PSD, == gram(species, weights)" begin
        m = fwd_mol()
        G = gram_matrix(m, fwd_q, fwd_lmax, fwd_E; chunk = fwd_chunk)
        @test size(G) == (5, 5, length(fwd_q))
        for k in axes(G, 3)
            @test issymmetric(G[:, :, k])
            @test minimum(eigvals(G[:, :, k])) > -1e-9
        end
        G_ref = gram(collect(species_multipoles(m, fwd_q, fwd_lmax, fwd_E; chunk = fwd_chunk)),
                     partial_wave_weights(fwd_lmax))
        @test G == G_ref
    end

    @testset "forward(G, m, c, dns, ρ) == hand-assembled I_calc" begin
        m   = fwd_mol()
        G   = gram_matrix(m, fwd_q, fwd_lmax, fwd_E; chunk = fwd_chunk)
        pars = (0.9, -0.4, 0.331, (1.2, 0.8, -0.3))
        v    = contrast_vector(pars[3], pars[4])
        ref  = intensity_calc(intensity(G, v), pars[1], pars[2])
        @test forward(G, pars...) == ref
    end

    @testset "forward(mol, …) convenience == gram_matrix + forward(G, …)" begin
        m = fwd_mol()
        G = gram_matrix(m, fwd_q, fwd_lmax, fwd_E; chunk = fwd_chunk)
        got = forward(m, fwd_q, fwd_lmax, fwd_E;
                      m = 1.7, c = 2.5, dns = 0.334, ρ = (1.0, 1.0, 0.0),
                      chunk = fwd_chunk)
        @test got ≈ forward(G, 1.7, 2.5, 0.334, (1.0, 1.0, 0.0))
    end

    @testset "detector map: m scales, c offsets" begin
        m = fwd_mol()
        G = gram_matrix(m, fwd_q, fwd_lmax, fwd_E; chunk = fwd_chunk)
        base = forward(G, 1.0, 0.0, 0.334, (1.0, 1.0, 0.0))
        @test forward(G, 3.0, 0.0, 0.334, (1.0, 1.0, 0.0)) ≈ 3 .* base
        @test forward(G, 1.0, 7.0, 0.334, (1.0, 1.0, 0.0)) ≈ base .+ 7
        @test all(>=(0.0), base)                      # PSD ⇒ vᵀGv ≥ 0
    end

    @testset "a class outside `classes` is inert (its ρ_k does nothing)" begin
        m  = fwd_mol()
        G  = gram_matrix(m, fwd_q, fwd_lmax, fwd_E;
                         chunk = fwd_chunk, classes = (SASA.CONVEX, SASA.CONCAVE))
        a  = forward(G, 1.0, 0.0, 0.334, (1.0, 1.0, 0.0))
        b  = forward(G, 1.0, 0.0, 0.334, (1.0, 1.0, 42.0))   # cavity unbuilt ⇒ B=0
        @test a ≈ b
    end

    # ---------------------------------------------------------------------
    # CRYSOL's fitted excluded-volume radius r0
    # ---------------------------------------------------------------------

    @testset "mean_atomic_radius is the plain mean of the per-atom radii" begin
        mo = fwd_mol()
        r  = radii(mo)
        @test mean_atomic_radius(mo) ≈ sum(r) / length(r)
    end

    @testset "excluded_volume_factor: r0 == r_m is exactly the identity" begin
        rm = 1.62
        @test excluded_volume_factor(fwd_q, rm, rm) == ones(length(fwd_q))
    end

    @testset "excluded_volume_factor: q=0 scales the excluded volume by c1^3" begin
        rm, r0 = 1.62, 1.78
        g = excluded_volume_factor([0.0], rm, r0)
        @test g[1] ≈ (r0 / rm)^3
    end

    @testset "excluded_volume_factor: exact for a dummy at the mean radius" begin
        # The single-envelope approximation is exact for an atom of radius r_m:
        # G(q)*f(V_m, q) must equal the directly-expanded dummy f(c1^3*V_m, q).
        rm, r0 = 1.62, 1.71
        c1 = r0 / rm
        vm = (4π / 3) * rm^3
        gauss(v, q) = v * exp(-q^2 * v^(2 / 3) / (4π))
        g = excluded_volume_factor(fwd_q, rm, r0)
        for (k, q) in enumerate(fwd_q)
            @test g[k] * gauss(vm, q) ≈ gauss(c1^3 * vm, q)
        end
    end

    @testset "excluded_volume_factor: r0 > r_m damps with q, r0 < r_m lifts" begin
        rm = 1.62
        up   = excluded_volume_factor(fwd_q, rm, 1.8)
        down = excluded_volume_factor(fwd_q, rm, 1.4)
        # both start at c1^3 and move monotonically in q away from it
        @test issorted(up ./ up[1];   rev = true)
        @test issorted(down ./ down[1])
    end

    @testset "excluded_volume_factor: domain errors" begin
        @test_throws DomainError excluded_volume_factor(fwd_q, 0.0, 1.6)
        @test_throws DomainError excluded_volume_factor(fwd_q, 1.6, -1.0)
    end

    @testset "contrast_matrix scales only the ex species" begin
        g = excluded_volume_factor(fwd_q, 1.62, 1.75)
        v = contrast_vector(0.334, (1.2, 0.8, -0.3))
        V = contrast_matrix(0.334, (1.2, 0.8, -0.3), g)
        @test size(V) == (5, length(fwd_q))
        for k in eachindex(g)
            @test V[2, k] ≈ v[2] * g[k]
            for a in (1, 3, 4, 5)
                @test V[a, k] == v[a]          # every other species untouched
            end
        end
    end

    @testset "intensity(G, V) reduces to intensity(G, v) on a constant V" begin
        mo = fwd_mol()
        G  = gram_matrix(mo, fwd_q, fwd_lmax, fwd_E; chunk = fwd_chunk)
        v  = contrast_vector(0.334, (1.0, 1.0, 0.0))
        V  = repeat(v, 1, length(fwd_q))
        @test intensity(G, V) ≈ intensity(G, v)
    end

    @testset "forward_cache carries G, the q grid and r_m" begin
        mo = fwd_mol()
        fc = forward_cache(mo, fwd_q, fwd_lmax, fwd_E; chunk = fwd_chunk)
        @test fc isa ForwardCache
        @test fc.G == gram_matrix(mo, fwd_q, fwd_lmax, fwd_E; chunk = fwd_chunk)
        @test fc.qvals == collect(Float64, fwd_q)
        @test fc.r_m == mean_atomic_radius(mo)
    end

    @testset "forward(cache, …): r0 = nothing / r_m is the uncorrected model" begin
        mo = fwd_mol()
        fc = forward_cache(mo, fwd_q, fwd_lmax, fwd_E; chunk = fwd_chunk)
        ref = forward(fc.G, 1.7, 2.5, 0.334, (1.0, 1.0, 0.0))
        @test forward(fc, 1.7, 2.5, 0.334, (1.0, 1.0, 0.0)) == ref
        @test forward(fc, 1.7, 2.5, 0.334, (1.0, 1.0, 0.0); r0 = nothing) == ref
        @test forward(fc, 1.7, 2.5, 0.334, (1.0, 1.0, 0.0); r0 = fc.r_m) == ref
    end

    @testset "forward(cache, …; r0) == hand-assembled q-dependent contraction" begin
        mo   = fwd_mol()
        fc   = forward_cache(mo, fwd_q, fwd_lmax, fwd_E; chunk = fwd_chunk)
        r0   = fc.r_m * 1.08
        pars = (0.9, -0.4, 0.331, (1.2, 0.8, -0.3))
        g    = excluded_volume_factor(fc.qvals, fc.r_m, r0)
        V    = contrast_matrix(pars[3], pars[4], g)
        ref  = intensity_calc(intensity(fc.G, V), pars[1], pars[2])
        @test forward(fc, pars...; r0 = r0) == ref
    end

    @testset "r0 actually changes the curve, and stays non-negative" begin
        mo = fwd_mol()
        fc = forward_cache(mo, fwd_q, fwd_lmax, fwd_E; chunk = fwd_chunk)
        base = forward(fc, 1.0, 0.0, 0.334, (1.0, 1.0, 0.0))
        big  = forward(fc, 1.0, 0.0, 0.334, (1.0, 1.0, 0.0); r0 = fc.r_m * 1.15)
        @test !(big ≈ base)
        @test all(>=(0.0), big)          # PSD ⇒ v(q)ᵀ G v(q) ≥ 0 at every r0
    end

    @testset "forward(mol, …; r0) convenience == cache + forward(cache, …; r0)" begin
        mo = fwd_mol()
        fc = forward_cache(mo, fwd_q, fwd_lmax, fwd_E; chunk = fwd_chunk)
        r0 = fc.r_m * 0.93
        @test forward(mo, fwd_q, fwd_lmax, fwd_E;
                      m = 1.7, c = 2.5, dns = 0.334, ρ = (1.0, 1.0, 0.0),
                      r0 = r0, chunk = fwd_chunk) ≈
              forward(fc, 1.7, 2.5, 0.334, (1.0, 1.0, 0.0); r0 = r0)
    end

    @testset "the whole fit-parameter path is AD-differentiable" begin
        # r0/dns/ρ are HMC parameters in stage 1, so a Dual must survive the
        # contrast -> contraction -> detector-map chain without a Float64 cast.
        mo = fwd_mol()
        fc = forward_cache(mo, fwd_q, fwd_lmax, fwd_E; chunk = fwd_chunk)
        f(p) = sum(forward(fc, p[1], p[2], p[3], (p[4], p[5], p[6]); r0 = p[7]))
        p0 = [1.7, 2.5, 0.334, 1.0, 1.0, 0.0, fc.r_m * 1.05]
        g  = ForwardDiff.gradient(f, p0)
        @test all(isfinite, g)
        # finite-difference check on r0, the new parameter
        h  = 1e-6
        pp = copy(p0); pp[7] += h
        pm = copy(p0); pm[7] -= h
        @test g[7] ≈ (f(pp) - f(pm)) / (2h) rtol = 1e-5
    end
end
