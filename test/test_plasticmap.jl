using .PlasticMap: Vec3, plastic_points, plastic_points!, PLASTIC_RATIO

nrm(p) = sqrt(p[1]^2 + p[2]^2 + p[3]^2)
dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]

@testset "PlasticMap" begin

    @testset "PLASTIC_RATIO is the real root of x³ = x + 1" begin
        ρ = PLASTIC_RATIO
        @test 1.0 < ρ < 2.0
        @test check_float(ρ^3, ρ + 1.0)
    end

    @testset "plastic_points: shape, type, edges" begin
        p = plastic_points(50)
        @test p isa Vector{Vec3}
        @test eltype(p) === NTuple{3,Float64}
        @test length(p) == 50
        @test isempty(plastic_points(0))
        @test plastic_points(1) == [plastic_points(3)[1]]
    end

    @testset "deterministic" begin
        @test plastic_points(200) == plastic_points(200)
    end

    @testset "every point is a unit vector" begin
        for p in plastic_points(5_000)
            @test abs(nrm(p) - 1.0) < 1e-12
        end
    end

    @testset "prefix stability: coarse set is an exact prefix of any finer set" begin
        fine = plastic_points(1_000)
        for k in (0, 1, 2, 7, 64, 257, 999)
            @test plastic_points(k) == fine[1:k]
        end
    end

    @testset "plastic_points! grows in place and matches a fresh allocation" begin
        buf = plastic_points(64)
        old = copy(buf)
        ret = plastic_points!(buf, 64, 256)
        @test ret === buf                     # grows the same array, no realloc of the handle
        @test length(buf) == 256
        @test buf[1:64] == old                # existing entries untouched
        @test buf == plastic_points(256)      # tail is the real continuation
    end

    @testset "repeated grows equal one shot" begin
        buf = Vector{Vec3}()
        for want in (10, 10, 11, 50, 200, 200, 512)
            have = length(buf)
            plastic_points!(buf, have, want)
            @test length(buf) == max(have, want)
        end
        @test length(buf) == 512
        @test buf == plastic_points(512)
    end

    @testset "no-op when want <= have" begin
        buf = plastic_points(128)
        snap = copy(buf)
        @test plastic_points!(buf, 128, 128) === buf
        @test buf == snap
        @test plastic_points!(buf, 128, 30) === buf
        @test buf == snap                     # length unchanged, nothing dropped
    end

    @testset "error contract" begin
        buf = plastic_points(10)
        @test_throws DomainError   plastic_points!(Vector{Vec3}(), -1, 5)
        @test_throws DomainError   plastic_points!(buf, 10, -1)
        @test_throws ArgumentError plastic_points!(buf, 9, 20)    # have too small for buf
        @test_throws ArgumentError plastic_points!(buf, 11, 20)   # have too large for buf
        @test_throws ArgumentError plastic_points!(Vector{Vec3}(), 3, 9)
    end

    @testset "even coverage (low-discrepancy sanity)" begin
        n = 4_096
        pts = plastic_points(n)

        # 1. centroid of a uniform sphere sample sits near the origin
        mx = sum(p -> p[1], pts) / n
        my = sum(p -> p[2], pts) / n
        mz = sum(p -> p[3], pts) / n
        @test sqrt(mx^2 + my^2 + mz^2) < 0.02

        # 2. octants are balanced
        oct = zeros(Int, 8)
        for p in pts
            oct[1 + (p[1] > 0) + 2 * (p[2] > 0) + 4 * (p[3] > 0)] += 1
        end
        @test all(c -> abs(c - n / 8) < 0.12 * n / 8, oct)

        # 3. height z is uniform on [-1, 1) (the projection is equal-area by construction)
        zb = zeros(Int, 10)
        for p in pts
            zb[clamp(floor(Int, (p[3] + 1.0) / 2.0 * 10) + 1, 1, 10)] += 1
        end
        @test all(c -> abs(c - n / 10) < 0.15 * n / 10, zb)

        # 4. spherical caps hold their area fraction, for a few axes and sizes
        for u in (  Vec3((1.0, 0.0, 0.0)), 
                    Vec3((0.0, 0.0, 1.0)), 
                    Vec3((1, 1, 1) ./ sqrt(3))
        )
            for freq in (0.1, 0.25, 0.5)
                cosθ = 1.0 - 2.0 * freq                 # cap of area fraction `freq`
                hit = count(p -> dot3(p, u) >= cosθ, pts) / n
                @test abs(hit - freq) < 0.03
            end
        end
    end

    @testset "type stability" begin
        buf = Vector{Vec3}()
        @inferred plastic_points(16)
        @inferred plastic_points!(buf, 0, 16)
    end
end
