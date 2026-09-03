using .Molecules: create, r, theta, phi, radii, vols, coords, MoleculeError

@testset "Molecules" begin
    @testset "two atoms on x axis" begin
        m = create("test", ["h", "h"], [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)])
        @test length(r(m)) == 2
        @test check_float(r(m)[1], 1.0) && check_float(r(m)[2], 1.0)
        @test length(theta(m)) == 2 && length(phi(m)) == 2
        @test check_float(theta(m)[1], pi / 2) && check_float(phi(m)[1], 0.0)
        @test check_float(theta(m)[2], pi / 2) && check_float(phi(m)[2], pi)
    end

    @testset "centering shifts to centroid" begin
        m = create("test", ["h", "h"], [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)])
        @test check_float(r(m)[1], 1.0) && check_float(r(m)[2], 1.0)
    end

    @testset "single atom r = 0, theta not NaN" begin
        m = create("test", ["h"], [(5.0, 5.0, 5.0)])
        @test check_float(r(m)[1], 0.0)
        @test !isnan(theta(m)[1])
    end

    @testset "empty coords raises" begin
        @test_throws MoleculeError create("empty", String[], NTuple{3,Float64}[])
    end

    @testset "vols for known elements" begin
        m = create(
            "test", ["fe", "o", "rn"],
            [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
        )
        @test length(vols(m)) == 3
        ev(x) = (4.0 / 3.0) * pi * x^3
        @test check_float(vols(m)[1], ev(1.274))
        @test check_float(vols(m)[3], ev(2.24))
    end

    @testset "unknown element raises only when vols is forced" begin
        m = create("test", ["zzzz"], [(0.0, 0.0, 0.0)])   # must not raise
        @test_throws MoleculeError vols(m)
    end

    @testset "create computes r/theta/phi/vols; repeat access stable" begin
        m = create(
            "test", ["o", "h", "h"],
            [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
        )
        @test length(r(m)) == 3
        @test check_float(r(m)[1], 0.4714045208)
        @test length(theta(m)) == 3 && length(phi(m)) == 3 && length(vols(m)) == 3
        @test radii(m) === radii(m)          # lazy accessor: memoized, same object on repeat
    end

    @testset "length mismatch raises" begin
        @test_throws MoleculeError create("bad", ["o", "h"], [(0.0, 0.0, 0.0)])
    end
end
