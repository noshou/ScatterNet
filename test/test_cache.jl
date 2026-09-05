# Exercises src/Molecule/Cache.jl (`Lazy` / `make` / `force`), which is
# `include`d straight into `Molecules` rather than living in its own module.
using .Molecules: Lazy, make, force

# `done` is the private "has the thunk run?" flag; reading it is the only way to
# assert laziness without also forcing the value we are trying to prove unforced.
_forced(c) = getfield(c, :done)

@testset "Cache" begin
    @testset "make builds an unforced Lazy{T}" begin
        c = make(Int, () -> 1)
        @test c isa Lazy{Int}
        @test !_forced(c)
        @test force(c) == 1
        @test _forced(c)
    end

    @testset "force memoizes: thunk runs exactly once" begin
        calls = Ref(0)
        c = make(Int, () -> (calls[] += 1; 42))
        a = force(c); b = force(c)
        @test a == 42 && b == 42 && calls[] == 1
        # ...and stays at one however many more times it is forced
        for _ in 1:10; force(c); end
        @test calls[] == 1
    end

    @testset "memoized value is the identical object on repeat" begin
        calls = Ref(0)
        c = make(Vector{Float64}, () -> (calls[] += 1; [1.0, 2.0]))
        v = force(c)
        @test v === force(c) === force(c)     # not merely ==: the same object
        push!(v, 3.0)                          # mutating the memo is visible
        @test force(c) == [1.0, 2.0, 3.0]
        @test calls[] == 1
    end

    @testset "make does not run the thunk until force" begin
        calls = Ref(0)
        c = make(Int, () -> (calls[] += 1; 1))
        @test calls[] == 0
        @test !_forced(c)
        force(c)
        @test calls[] == 1
    end

    @testset "the thunk's value is type-asserted against T" begin
        # `c.value = c.f()::T` is an assertion, not a convert: an Int-valued
        # thunk in a Lazy{Float64} is a TypeError, it is not promoted to 1.0.
        @test_throws TypeError force(make(Float64, () -> 1))
        # an abstract T admits any subtype, unconverted
        c = make(Real, () -> 1)
        @test force(c) === 1
        @test force(make(Any, () -> "x")) == "x"
        # the failed assertion must not have marked the cache as done
        bad = make(Float64, () -> 1)
        @test_throws TypeError force(bad)
        @test !_forced(bad)
    end

    @testset "a throwing thunk propagates and leaves the cache unforced" begin
        calls = Ref(0)
        c = make(Int, () -> (calls[] += 1; calls[] == 1 ? error("boom") : 99))
        @test_throws ErrorException force(c)
        @test calls[] == 1
        @test !_forced(c)                      # lock released, nothing memoized
        # documented consequence of not memoizing failures: the next force retries
        @test force(c) == 99
        @test calls[] == 2
        @test _forced(c)
        @test force(c) == 99 && calls[] == 2    # and the retry's value is memoized
    end

    @testset "concurrent force from multiple tasks runs the thunk once" begin
        calls = Threads.Atomic{Int}(0)
        c = make(Int, () -> (Threads.atomic_add!(calls, 1); 7))
        tasks = [Threads.@spawn force(c) for _ in 1:8]   # spawn all, then join all
        rs = fetch.(tasks)
        @test all(==(7), rs) && calls[] == 1
    end

    @testset "concurrent force under real contention still runs the thunk once" begin
        # a slow thunk guarantees the other tasks arrive while the lock is held,
        # which is the case the ReentrantLock in `Lazy` actually exists for
        calls = Threads.Atomic{Int}(0)
        c = make(Vector{Int}, () -> (Threads.atomic_add!(calls, 1); sleep(0.05); [1, 2, 3]))
        results = Vector{Vector{Int}}(undef, 16)
        @sync for i in 1:16
            Threads.@spawn results[i] = force(c)
        end
        @test calls[] == 1
        @test all(r -> r === results[1], results)   # every task saw the same object
        @test results[1] == [1, 2, 3]
    end

    @testset "distinct caches are independent" begin
        @test force(make(Int, () -> 1)) == 1
        @test force(make(Int, () -> 2)) == 2
        a = make(Vector{Int}, () -> Int[]); b = make(Vector{Int}, () -> Int[])
        @test force(a) !== force(b)          # equal values, separate objects
        # forcing one cache leaves an unrelated one untouched
        c = make(Vector{Int}, () -> Int[])
        force(a)
        @test !_forced(c)
    end

    @testset "force is type stable in T" begin
        @test @inferred(force(make(Int, () -> 1))) === 1
        @test @inferred(force(make(Float64, () -> 1.5))) === 1.5
        @test @inferred(force(make(Vector{Float64}, () -> [1.0]))) == [1.0]
    end
end
