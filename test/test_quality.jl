# Not an OCaml port: package hygiene + type-stability guards.
using Aqua, JET
using .SphFuncs: sphHarm, sphBess, legendre_sphPlm
using .Molecules: create, r, radii, vols
using .AtomicRadii: resolve_one, lookup

@testset "Aqua" begin
    Aqua.test_all(ScatterNet; ambiguities = false)
    Aqua.test_ambiguities(ScatterNet)
end

@testset "type stability (@inferred)" begin
    θ = collect(range(0.1, pi - 0.1; length = 8)); φ = collect(range(0.0, 2pi; length = 8))
    @inferred sphHarm(4, θ, φ)
    @inferred sphBess([1.0, 2.0], [0.1, 0.5, 1.0], 4)
    @inferred legendre_sphPlm(3, 2, 0.5)
    @inferred Union{Float64,Nothing} resolve_one("fe3+")
    @inferred lookup(["fe3+", "o2-"])
    m = @inferred create("t", ["o", "h", "h"],
                         [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)])
    @inferred r(m); @inferred radii(m); @inferred vols(m)
end

@testset "JET (focused type-stability analysis)" begin
    @test_opt target_modules = (SphFuncs,) sphBess([1.0, 2.0], [0.1, 0.5], 3)
    @test_opt target_modules = (SphFuncs,) sphHarm(3, [0.4, 1.2], [0.1, 2.0])
    @test_opt target_modules = (Molecules,) create("t", ["o", "h"],
        [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
    @test_opt target_modules = (AtomicRadii,) resolve_one("fe3+")
end
