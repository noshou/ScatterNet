using .Cache: make, force

@testset "Cache" begin
    @testset "force memoizes: thunk runs exactly once" begin
        calls = Ref(0)
        c = make(Int, () -> (calls[] += 1; 42))
        a = force(c); b = force(c)
        @test a == 42 && b == 42 && calls[] == 1
    end

    @testset "make does not run the thunk until force" begin
        calls = Ref(0)
        _ = make(Int, () -> (calls[] += 1; 1))
        @test calls[] == 0
    end

    @testset "concurrent force from multiple tasks runs the thunk once" begin
        calls = Threads.Atomic{Int}(0)
        c = make(Int, () -> (Threads.atomic_add!(calls, 1); 7))
        tasks = [Threads.@spawn force(c) for _ in 1:8]   # spawn all, then join all
        rs = fetch.(tasks)
        @test all(==(7), rs) && calls[] == 1
    end

    @testset "distinct caches are independent" begin
        @test force(make(Int, () -> 1)) == 1
        @test force(make(Int, () -> 2)) == 2
    end
end
