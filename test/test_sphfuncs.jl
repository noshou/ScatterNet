# 1:1 port of scattering/test/test_sphfuncs.ml
using .SphFuncs: sphHarm, sphBess, SphHarmError, SphBessError

y00 = 1.0 / (2.0 * sqrt(pi))
y10(θ) = sqrt(3.0 / (4.0 * pi)) * cos(θ)
y11(θ, φ) = (-(sqrt(3.0 / (8.0 * pi))) * sin(θ)) * cis(φ)
idx(l, m) = l * (l + 1) ÷ 2 + m + 1

@testset "SphFuncs" begin
    @testset "sphHarm known values" begin
        for (θv, φv) in ((pi / 2, 0.0), (0.0, 0.7))
            y = sphHarm(1, [θv], [φv])
            @test check_complex(complex(y00), y[idx(0, 0), 1])
            @test check_complex(complex(y10(θv)), y[idx(1, 0), 1])
            @test check_complex(y11(θv, φv), y[idx(1, 1), 1])
        end
    end

    @testset "sphHarm exception contract" begin
        @test_throws SphHarmError sphHarm(-1, [1.0], [1.0])
        @test_throws SphHarmError sphHarm(2, [1.0], [1.0, 2.0])
        @test_throws SphHarmError sphHarm(2, Float64[], Float64[])
        @test_throws SphHarmError sphHarm(2, zeros(2, 2), [1.0])
    end

    @testset "sphBess known values" begin
        j0(x) = x == 0.0 ? 1.0 : sin(x) / x
        j1(x) = x == 0.0 ? 0.0 : sin(x) / (x * x) - cos(x) / x
        for x in (1.0, 0.0)
            j = sphBess([x], [1.0], 1)
            @test check_float(j0(x), j[1, 1, 1])
            @test check_float(j1(x), j[2, 1, 1])
        end
    end

    @testset "sphBess exception contract" begin
        @test_throws SphBessError sphBess(Float64[], [1.0], 1)
        @test_throws SphBessError sphBess([1.0], Float64[], 1)
        @test_throws SphBessError sphBess([-1.0], [1.0], 1)
        @test_throws SphBessError sphBess([1.0], [-1.0], 1)
        @test_throws SphBessError sphBess([1.0], [1.0], -1)
    end
end
