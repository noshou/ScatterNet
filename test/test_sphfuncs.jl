# Exercises src/Scattering/SphFuncs.jl against closed forms and exact identities
# (Unsold's theorem, the Y = P̄ e^{imφ} definition, the Bessel recurrence),
# not just tabulated numbers.
using .SphFuncs: sphHarm, sphBess, legendre_sphPlm, SphHarmError, SphBessError

y00 = 1.0 / (2.0 * sqrt(π))
y10(θ) = sqrt(3.0 / (4.0 * π)) * cos(θ)
y11(θ, φ) = (-(sqrt(3.0 / (8.0 * π))) * sin(θ)) * cis(φ)
# l = 2, Condon-Shortley phase included
y20(θ) = sqrt(5.0 / (16.0 * π)) * (3.0 * cos(θ)^2 - 1.0)
y21(θ, φ) = -sqrt(15.0 / (8.0 * π)) * sin(θ) * cos(θ) * cis(φ)
y22(θ, φ) = sqrt(15.0 / (32.0 * π)) * sin(θ)^2 * cis(2φ)
idx(l, m) = l * (l + 1) ÷ 2 + m + 1
nrows(lMax) = (lMax + 1) * (lMax + 2) ÷ 2

# closed forms for the first few spherical Bessel functions
j0(x) = x == 0.0 ? 1.0 : sin(x) / x
j1(x) = x == 0.0 ? 0.0 : sin(x) / (x * x) - cos(x) / x
j2(x) = x == 0.0 ? 0.0 : (3.0 / x^3 - 1.0 / x) * sin(x) - 3.0 * cos(x) / x^2

@testset "SphFuncs" begin

    @testset "sphHarm known values" begin
        for (θv, φv) in ((π / 2, 0.0), (0.0, 0.7))
            y = sphHarm(1, [θv], [φv])
            @test check_complex(complex(y00), y[idx(0, 0), 1])
            @test check_complex(complex(y10(θv)), y[idx(1, 0), 1])
            @test check_complex(y11(θv, φv), y[idx(1, 1), 1])
        end
    end

    @testset "sphHarm known values through l = 2" begin
        for (θv, φv) in ((0.37, 0.9), (2.1, -1.3), (π / 2, π), (1.0, 0.0))
            y = sphHarm(2, [θv], [φv])
            @test check_complex(complex(y00),      y[idx(0, 0), 1])
            @test check_complex(complex(y10(θv)),  y[idx(1, 0), 1])
            @test check_complex(y11(θv, φv),       y[idx(1, 1), 1])
            @test check_complex(complex(y20(θv)),  y[idx(2, 0), 1])
            @test check_complex(y21(θv, φv),       y[idx(2, 1), 1])
            @test check_complex(y22(θv, φv),       y[idx(2, 2), 1])
        end
    end

    @testset "sphHarm packing: shape and row index" begin
        θ = [0.3, 1.1, 2.7]; φ = [0.2, -1.0, 3.0]
        for lMax in 0:5
            y = sphHarm(lMax, θ, φ)
            @test y isa Matrix{ComplexF64}
            @test size(y) == (nrows(lMax), 3)
        end
        # the packing is dense and collision-free over the (l, m) triangle
        seen = [idx(l, m) for l in 0:5 for m in 0:l]
        @test sort(seen) == collect(1:nrows(5))
        # a bigger lMax only appends rows; the shared prefix is bit-identical
        @test sphHarm(2, θ, φ) == sphHarm(5, θ, φ)[1:nrows(2), :]
    end

    @testset "sphHarm columns are independent points" begin
        θ = [0.3, 1.1, 2.7, 0.0, π]; φ = [0.2, -1.0, 3.0, 0.5, 1.5]
        y = sphHarm(4, θ, φ)
        for i in eachindex(θ)
            @test y[:, i] == sphHarm(4, [θ[i]], [φ[i]])[:, 1]
        end
    end

    @testset "sphHarm: Unsold's theorem, sum_m |Y_lm|^2 = (2l+1)/4pi" begin
        # Y_{l,-m} = (-1)^m conj(Y_{lm}), so the stored m >= 0 half determines
        # the full sum; this pins normalization AND phase convention exactly.
        θ = [0.05, 0.7, 1.5707, 2.4, 3.09]; φ = [0.0, 1.2, -2.5, 3.0, 0.4]
        lMax = 6
        y = sphHarm(lMax, θ, φ)
        for i in eachindex(θ), l in 0:lMax
            s = abs2(y[idx(l, 0), i]) + 2 * sum(abs2(y[idx(l, m), i]) for m in 1:l; init = 0.0)
            @test check_float(s, (2l + 1) / (4π))
        end
    end

    @testset "sphHarm: Y_lm = legendre_sphPlm(l, m, cos θ) * exp(i m φ)" begin
        θ = [0.3, 1.1, 2.7]; φ = [0.2, -1.0, 3.0]
        y = sphHarm(4, θ, φ)
        for i in eachindex(θ), l in 0:4, m in 0:l
            @test check_complex(y[idx(l, m), i], legendre_sphPlm(l, m, cos(θ[i])) * cis(m * φ[i]))
        end
    end

    @testset "sphHarm: m = 0 harmonics are real and phi-independent" begin
        for l in 0:4
            a = sphHarm(4, [0.8], [0.0])[idx(l, 0), 1]
            b = sphHarm(4, [0.8], [2.5])[idx(l, 0), 1]
            @test check_float(imag(a), 0.0)
            @test check_complex(a, b)
        end
    end

    @testset "sphHarm: phi enters only as the phase exp(i m phi)" begin
        θ = [0.8, 2.2]
        y0 = sphHarm(4, θ, [0.0, 0.0])
        for φv in (0.3, 1.9, -2.2, 2π)
            y = sphHarm(4, θ, [φv, φv])
            for i in 1:2, l in 0:4, m in 0:l
                @test check_complex(y[idx(l, m), i], y0[idx(l, m), i] * cis(m * φv))
                @test check_float(abs(y[idx(l, m), i]), abs(y0[idx(l, m), i]))
            end
        end
    end

    @testset "sphHarm at the poles: only m = 0 survives" begin
        for θv in (0.0, π)
            y = sphHarm(4, [θv], [0.9])
            for l in 0:4
                @test check_float(abs(y[idx(l, 0), 1]), sqrt((2l + 1) / (4π)))
                for m in 1:l
                    @test check_float(abs(y[idx(l, m), 1]), 0.0)
                end
            end
        end
    end

    @testset "sphHarm accepts non-Float64 and non-Vector 1-D inputs" begin
        @test sphHarm(2, [0, 1], [0, 1]) ≈ sphHarm(2, [0.0, 1.0], [0.0, 1.0])
        @test sphHarm(2, 0.0:1.0:1.0, 0.0:1.0:1.0) ≈ sphHarm(2, [0.0, 1.0], [0.0, 1.0])
        @test all(isfinite, sphHarm(6, [0.0, π, 1e-12, π - 1e-12], fill(0.3, 4)))
    end

    @testset "sphHarm exception contract" begin
        @test_throws SphHarmError sphHarm(-1, [1.0], [1.0])
        @test_throws SphHarmError sphHarm(-5, [1.0], [1.0])
        @test_throws SphHarmError sphHarm(2, [1.0], [1.0, 2.0])
        @test_throws SphHarmError sphHarm(2, [1.0, 2.0], [1.0])
        @test_throws SphHarmError sphHarm(2, Float64[], Float64[])
        @test_throws SphHarmError sphHarm(2, Float64[], [1.0])
        @test_throws SphHarmError sphHarm(2, [1.0], Float64[])
        @test_throws SphHarmError sphHarm(2, zeros(2, 2), [1.0])
        @test_throws SphHarmError sphHarm(2, [1.0], zeros(2, 2))
        # lMax = 0 is the boundary of the valid range, not an error
        @test size(sphHarm(0, [1.0], [1.0])) == (1, 1)
    end

    @testset "sphBess known values" begin
        for x in (1.0, 0.0)
            j = sphBess([x], [1.0], 1)
            @test check_float(j0(x), j[1, 1, 1])
            @test check_float(j1(x), j[2, 1, 1])
        end
    end

    @testset "sphBess closed forms for l = 0, 1, 2" begin
        rs = [0.0, 0.5, 1.0, 3.7]; qs = [0.0, 0.25, 2.0]
        j = sphBess(rs, qs, 2)
        for (qi, q) in enumerate(qs), (ri, r) in enumerate(rs)
            x = q * r
            @test check_float(j[1, qi, ri], j0(x))
            @test check_float(j[2, qi, ri], j1(x))
            @test check_float(j[3, qi, ri], j2(x))
        end
    end

    @testset "sphBess shape is (lMax+1, |q|, |r|)" begin
        j = sphBess([1.0, 2.0, 3.0], [0.1, 0.2], 3)
        @test j isa Array{Float64,3}
        @test size(j) == (4, 2, 3)
        @test size(sphBess([1.0], [1.0], 0)) == (1, 1, 1)
        # a bigger lMax only appends the leading dimension
        @test sphBess([1.0, 2.0], [0.5], 2) == sphBess([1.0, 2.0], [0.5], 5)[1:3, :, :]
    end

    @testset "sphBess depends only on the product q*r" begin
        @test sphBess([2.0], [3.0], 4)[:, 1, 1] == sphBess([3.0], [2.0], 4)[:, 1, 1]
        @test sphBess([1.0], [6.0], 4)[:, 1, 1] == sphBess([6.0], [1.0], 4)[:, 1, 1]
        # and the (q, r) grid really is the outer product
        j = sphBess([1.0, 2.0], [3.0, 5.0], 3)
        for (qi, q) in enumerate([3.0, 5.0]), (ri, r) in enumerate([1.0, 2.0])
            @test j[:, qi, ri] == sphBess([q * r], [1.0], 3)[:, 1, 1]
        end
    end

    @testset "sphBess: j_l(0) = delta_{l0}" begin
        for (rs, qs) in (([0.0], [1.0]), ([1.0], [0.0]), ([0.0], [0.0]))
            j = sphBess(rs, qs, 5)
            @test check_float(j[1, 1, 1], 1.0)
            @test all(l -> check_float(j[l, 1, 1], 0.0), 2:6)
        end
    end

    @testset "sphBess satisfies the three-term recurrence" begin
        # j_{l-1}(x) + j_{l+1}(x) = (2l+1)/x * j_l(x)
        lMax = 8
        for x in (0.1, 1.0, 2.5, 7.3, 30.0)
            j = sphBess([x], [1.0], lMax)
            for l in 1:(lMax - 1)
                @test check_float(j[l, 1, 1] + j[l + 2, 1, 1], (2l + 1) / x * j[l + 1, 1, 1])
            end
        end
    end

    @testset "sphBess: bounds and small-x asymptotics" begin
        j = sphBess([1.0], collect(0.05:0.37:10.0), 6)
        @test all(isfinite, j)
        @test all(v -> abs(v) <= 1.0 + 1e-12, j)      # |j_l(x)| <= 1 for real x >= 0
        # j_l(x) -> x^l / (2l+1)!! as x -> 0
        x = 1.0e-3
        js = sphBess([x], [1.0], 4)
        dfact = 1.0
        for l in 0:4
            l > 0 && (dfact *= (2l + 1))
            @test isapprox(js[l + 1, 1, 1], x^l / dfact; rtol = 1e-5)
        end
    end

    @testset "sphBess accepts non-Float64 inputs" begin
        @test sphBess([1, 2], [1], 2) ≈ sphBess([1.0, 2.0], [1.0], 2)
        @test sphBess(0.0:1.0:2.0, 1.0:1.0:1.0, 2) ≈ sphBess([0.0, 1.0, 2.0], [1.0], 2)
    end

    @testset "sphBess exception contract" begin
        @test_throws SphBessError sphBess(Float64[], [1.0], 1)
        @test_throws SphBessError sphBess([1.0], Float64[], 1)
        @test_throws SphBessError sphBess(Float64[], Float64[], 1)
        @test_throws SphBessError sphBess([-1.0], [1.0], 1)
        @test_throws SphBessError sphBess([1.0, -1e-12], [1.0], 1)   # any negative entry
        @test_throws SphBessError sphBess([1.0], [-1.0], 1)
        @test_throws SphBessError sphBess([1.0], [1.0, -1.0], 1)
        @test_throws SphBessError sphBess([1.0], [1.0], -1)
        @test_throws SphBessError sphBess([1.0], [1.0], -7)
        # 0 is on the allowed side of every boundary
        @test size(sphBess([0.0], [0.0], 0)) == (1, 1, 1)
    end

    @testset "legendre_sphPlm closed forms" begin
        for x in (-1.0, -0.6, 0.0, 0.25, 1.0)
            @test check_float(legendre_sphPlm(0, 0, x), 1.0 / (2.0 * sqrt(π)))
            @test check_float(legendre_sphPlm(1, 0, x), sqrt(3.0 / (4.0 * π)) * x)
            # Condon-Shortley phase: P̄_1^1 is negative for x in (-1, 1)
            @test check_float(legendre_sphPlm(1, 1, x),
                              -sqrt(3.0 / (8.0 * π)) * sqrt(max(0.0, 1.0 - x^2)))
            @test check_float(legendre_sphPlm(2, 0, x),
                              sqrt(5.0 / (16.0 * π)) * (3x^2 - 1.0))
            @test check_float(legendre_sphPlm(2, 2, x),
                              sqrt(15.0 / (32.0 * π)) * (1.0 - x^2))
        end
        @test legendre_sphPlm(3, 2, 0.5) isa Float64
        @test legendre_sphPlm(2, 0, 1) isa Float64        # Integer x is accepted
    end
end
