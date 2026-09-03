using .AtomicRadii: lookup

lookup_one(ion) = lookup([ion])[1][2]

@testset "AtomicRadii" begin
    @testset "exact ion match, both token orderings" begin
        @test check_float(0.49, lookup_one("fe3+"))
        @test check_float(0.49, lookup_one("fe+3"))
        @test check_float(2.2, lookup_one("au1-"))
    end

    @testset "bare element falls through to atomic_radii" begin
        @test check_float(1.274, lookup_one("fe"))
        @test check_float(2.24, lookup_one("rn"))
    end

    @testset "nearest-charge fallback" begin           # fe has 2/3/4/6+, not 5+
        v5, v4 = lookup_one("fe5+"), lookup_one("fe4+")
        @test v5 !== nothing && v4 !== nothing && check_float(v5, v4)
    end

    @testset "total miss" begin
        @test lookup_one("zzzz9+") === nothing
    end

    @testset "batch: order/count preserved, repeats deduped" begin
        input = ["fe3+", "zzzz9+", "fe3+", "rn"]
        res = lookup(input)
        @test [k for (k, _) in res] == input
        r1, m, r2, r3 = (v for (_, v) in res)
        @test r1 !== nothing && m === nothing && r3 !== nothing && check_float(r1, r2)
    end
end
