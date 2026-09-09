# Exercises src/Scattering/Intensity.jl: the species Gram matrix `gram`, the
# contrast contraction `intensity` = v' G v, the detector map `intensity_calc`,
# and `contrast_vector`.
#
# The Gram / quadratic-form algebra is checked against hand-written 6-term
# (3-species) and 15-term (5-species) expansions assembled directly from
# `self_scatter` / `cross_scatter`, i.e. the term tables in SCATTERNET_MATH.md
# "Multi-species amplitude algebra and the fitted parameters", subsections 2
# and 7. The `B_lm` inputs are deterministic pseudo-random arrays: the assembly
# is linear in them and does not depend on where they came from, so there is no
# pipeline input on the reference side.
using .Scattering: gram, intensity, intensity_calc, contrast_vector, DRO_UNIT,
                   partial_wave_weights, self_scatter, cross_scatter
using LinearAlgebra: Symmetric, eigvals, issymmetric

# A packed-(l,m) B_lm array (C, K, Q) filled from a deterministic LCG so every
# run and every reader sees the same numbers.
function iy_B(C::Int, lMax::Int, Q::Int; seed::Int = 0)
    K = (lMax + 1) * (lMax + 2) ÷ 2
    B = Array{ComplexF64,3}(undef, C, K, Q)
    s = UInt64(seed) + 0x9e3779b97f4a7c15
    nextf() = (s = 6364136223846793005 * s + 1442695040888963407; Float64(s >> 11) / 2.0^53 - 0.5)
    @inbounds for i in eachindex(B)
        B[i] = complex(nextf(), nextf())
    end
    return B
end

iy_lMax = 3
iy_Q    = 5
iy_w    = partial_wave_weights(iy_lMax)

@testset "Intensity" begin

    @testset "DRO_UNIT is CRYSOL's --dro contrast unit" begin
        @test DRO_UNIT == 0.03
    end

    @testset "gram: shape, symmetry, and the self / cross diagonal identity" begin
        Bs = [iy_B(2, iy_lMax, iy_Q; seed = 1),
              iy_B(1, iy_lMax, iy_Q; seed = 2),
              iy_B(1, iy_lMax, iy_Q; seed = 3)]
        G = gram(Bs, iy_w)
        @test G isa Array{Float64,3}
        @test size(G) == (3, 3, iy_Q)
        for k in 1:iy_Q
            @test issymmetric(G[:, :, k])
        end
        # the diagonal is exactly self_scatter, and self_scatter == cross_scatter(B, B)
        for a in 1:3
            s = self_scatter(Bs[a], iy_w)
            @test G[a, a, :] == s
            @test cross_scatter(Bs[a], Bs[a], iy_w) == s
        end
        # the off-diagonal is exactly cross_scatter
        @test G[1, 2, :] == cross_scatter(Bs[1], Bs[2], iy_w)
        @test G[1, 3, :] == cross_scatter(Bs[1], Bs[3], iy_w)
        @test G[2, 3, :] == cross_scatter(Bs[2], Bs[3], iy_w)
    end

    @testset "gram: every G(:,:,q) is positive semidefinite" begin
        Bs = [iy_B(2, iy_lMax, iy_Q; seed = 7),
              iy_B(2, iy_lMax, iy_Q; seed = 8),
              iy_B(1, iy_lMax, iy_Q; seed = 9),
              iy_B(1, iy_lMax, iy_Q; seed = 10),
              iy_B(1, iy_lMax, iy_Q; seed = 11)]
        G = gram(Bs, iy_w)
        for k in 1:iy_Q
            λ = eigvals(Symmetric(G[:, :, k]))
            @test minimum(λ) ≥ -1e-9
        end
    end

    @testset "gram: single-species and all-zero edge cases" begin
        B0 = zeros(ComplexF64, 1, length(iy_w), iy_Q)
        @test all(iszero, gram([B0], iy_w))
        B1 = iy_B(1, iy_lMax, iy_Q; seed = 15)
        G1 = gram([B1], iy_w)
        @test size(G1) == (1, 1, iy_Q)
        @test G1[1, 1, :] == self_scatter(B1, iy_w)
        # a zero species contributes a zero row/column, nothing else
        G = gram([B1, B0], iy_w)
        @test all(iszero, G[1, 2, :]) && all(iszero, G[2, 2, :])
        @test G[1, 1, :] == self_scatter(B1, iy_w)
    end

    @testset "gram accepts a tuple identically to a vector" begin
        B1 = iy_B(1, iy_lMax, iy_Q; seed = 21)
        B2 = iy_B(1, iy_lMax, iy_Q; seed = 22)
        @test gram((B1, B2), iy_w) == gram([B1, B2], iy_w)
    end

    @testset "gram: shape guards" begin
        good = iy_B(1, iy_lMax, iy_Q; seed = 25)
        @test_throws ArgumentError gram(typeof(good)[], iy_w)             # no species
        @test_throws ArgumentError gram([good], partial_wave_weights(iy_lMax + 1))  # K mismatch
        @test_throws ArgumentError gram([good, iy_B(1, iy_lMax, iy_Q + 1; seed = 26)], iy_w)  # Q mismatch
    end

    @testset "intensity: 3-species collapse equals the 6-term expansion" begin
        Bv = iy_B(2, iy_lMax, iy_Q; seed = 31)
        Be = iy_B(1, iy_lMax, iy_Q; seed = 32)
        Bs = iy_B(1, iy_lMax, iy_Q; seed = 33)
        dns, dro = 1.07, 0.031

        Svv = self_scatter(Bv, iy_w)
        See = self_scatter(Be, iy_w)
        Sss = self_scatter(Bs, iy_w)
        Sve = cross_scatter(Bv, Be, iy_w)
        Svs = cross_scatter(Bv, Bs, iy_w)
        Ses = cross_scatter(Be, Bs, iy_w)

        ref = Svv .- 2dns .* Sve .+ 2dro .* Svs .+
              dns^2 .* See .+ dro^2 .* Sss .- (2 * dns * dro) .* Ses

        got = intensity(gram([Bv, Be, Bs], iy_w), [1.0, -dns, dro])
        @test all(check_float.(got, ref))
    end

    @testset "intensity: 5-species equals the 15-term expansion" begin
        B = [iy_B(2, iy_lMax, iy_Q; seed = 40 + i) for i in 1:5]
        dns = 0.98
        ρ   = (1.0, 0.7, -0.3)
        d   = (DRO_UNIT * ρ[1], DRO_UNIT * ρ[2], DRO_UNIT * ρ[3])

        S(a, b) = a == b ? self_scatter(B[a], iy_w) : cross_scatter(B[a], B[b], iy_w)

        ref =  S(1, 1) .+ dns^2 .* S(2, 2) .+                                   # vac,vac  ex,ex
               d[1]^2 .* S(3, 3) .+ d[2]^2 .* S(4, 4) .+ d[3]^2 .* S(5, 5) .+    # shk,shk
               (-2dns) .* S(1, 2) .+                                            # vac,ex
               2d[1] .* S(1, 3) .+ 2d[2] .* S(1, 4) .+ 2d[3] .* S(1, 5) .+       # vac,shk
               (-2dns) .* (d[1] .* S(2, 3) .+ d[2] .* S(2, 4) .+ d[3] .* S(2, 5)) .+  # ex,shk
               2 .* (d[1] * d[2] .* S(3, 4) .+ d[1] * d[3] .* S(3, 5) .+ d[2] * d[3] .* S(4, 5))  # shj,shk

        v = [1.0, -dns, d[1], d[2], d[3]]
        @test v == contrast_vector(dns, ρ)
        got = intensity(gram(B, iy_w), v)
        @test all(check_float.(got, ref))
    end

    @testset "intensity: quadratic form is ≥ 0 for real v (PSD)" begin
        B = [iy_B(2, iy_lMax, iy_Q; seed = 60 + i) for i in 1:4]
        G = gram(B, iy_w)
        for v in ([1.0, -1.0, 0.03, 0.03],
                  [1.0, 0.5, -2.0, 3.0],
                  [0.3, -1.1, 2.2, -0.7],
                  zeros(4))
            @test all(≥(-1e-9), intensity(G, v))
        end
        @test all(iszero, intensity(G, zeros(4)))
    end

    @testset "intensity: bilinear in v" begin
        B = [iy_B(1, iy_lMax, iy_Q; seed = 81), iy_B(1, iy_lMax, iy_Q; seed = 82)]
        G = gram(B, iy_w)
        v = [1.0, -0.9]
        # I(t v) = t² I(v)
        @test all(check_float.(intensity(G, 2.0 .* v), 4.0 .* intensity(G, v)))
        @test_throws DimensionMismatch intensity(G, [1.0, 2.0, 3.0])
    end

    @testset "intensity_calc is affine in (m, c)" begin
        B = [iy_B(1, iy_lMax, iy_Q; seed = 71), iy_B(1, iy_lMax, iy_Q; seed = 72)]
        I = intensity(gram(B, iy_w), [1.0, -1.0])
        @test intensity_calc(I, 1.0, 0.0) == I
        @test all(check_float.(intensity_calc(I, 2.5, 0.0), 2.5 .* I))
        @test all(check_float.(intensity_calc(I, 1.0, 4.0), I .+ 4.0))
        @test all(check_float.(intensity_calc(I, 3.0, -1.5), 3.0 .* I .- 1.5))
        @test intensity_calc(I, 2, 0) isa Vector{Float64}   # integer m, c accepted
    end

    @testset "contrast_vector: values, length, species order" begin
        @test contrast_vector(0.334, (1.0, 1.0, 0.0)) == [1.0, -0.334, 0.03, 0.03, 0.0]
        @test contrast_vector(0.334, 1.0)             == [1.0, -0.334, 0.03]
        @test length(contrast_vector(0.3, (1.0, 1.0, 1.0))) == 5
        @test length(contrast_vector(0.3, 1.0))             == 3
        @test contrast_vector(0.3, (1.0, 0.0, 0.0))[4:5] == [0.0, 0.0]   # ρ = 0 zeroes a shell
        @test contrast_vector(0.334, (1, 1, 0)) isa Vector{Float64}       # integer ρ accepted
    end
end
