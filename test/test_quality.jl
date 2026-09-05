# package hygiene + type-stability guards.
using Aqua, JET
using .SphFuncs: sphHarm, sphBess, legendre_sphPlm
using .Molecules: create, coords_cartesian, coords_spherical, radii, vols, r_max,
                elms, name, Molecule
using .Molecules: resolve_one, _resolve_all, tryparse_ion, ion_key, nearest_ion
using ScatterNet.Molecule: SASA

@testset "Aqua" begin
    Aqua.test_all(ScatterNet; ambiguities = false)
    Aqua.test_ambiguities(ScatterNet)
end

@testset "type stability (@inferred)" begin
    θ = collect(range(0.1, π - 0.1; length = 8)); φ = collect(range(0.0, 2pi; length = 8))
    @inferred sphHarm(4, θ, φ)
    @inferred sphBess([1.0, 2.0], [0.1, 0.5, 1.0], 4)
    @inferred legendre_sphPlm(3, 2, 0.5)
    @inferred Union{Float64,Nothing} resolve_one("fe3+")
    @inferred _resolve_all(["fe3+", "o2-"])
    @inferred Union{Molecules.Ion,Nothing} tryparse_ion("fe3+")
    @inferred ion_key(Molecules.Ion("fe", 3))
    @inferred Union{String,Nothing} nearest_ion("fe", 5)

    m = @inferred create(
        "t", ["o", "h", "h"],
        [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
    )
    @test m isa Molecule
    @inferred coords_cartesian(m)
    @inferred coords_spherical(m)
    @inferred radii(m)
    @inferred vols(m)
    @inferred r_max(m)
    @inferred elms(m)
    @inferred name(m)
end

@testset "type stability of the SASA entry points (@inferred)" begin
    m = create("t", ["o", "h", "h"], [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)])
    @test (@inferred SASA.sasa_atoms(m; n_occ = 16, n_exp = 64, probe = 1.4)) isa Vector{Float64}
    @test (@inferred SASA.sasa(m; n_occ = 16, n_exp = 64, probe = 1.4)) isa Float64
end

@testset "JET (focused type-stability analysis)" begin
    @test_opt target_modules = (SphFuncs,) sphBess([1.0, 2.0], [0.1, 0.5], 3)
    @test_opt target_modules = (SphFuncs,) sphHarm(3, [0.4, 1.2], [0.1, 2.0])
    @test_opt target_modules = (SphFuncs,) legendre_sphPlm(3, 2, 0.5)
    @test_opt target_modules = (Molecules,) create("t", ["o", "h"],
        [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
    @test_opt target_modules = (Molecules,) resolve_one("fe3+")
    @test_opt target_modules = (Molecules,) _resolve_all(["fe3+", "o2-"])
    @test_opt target_modules = (Molecules,) tryparse_ion("fe3+")

    let m = create("t", ["o", "h", "h"], [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)])
        @test_opt target_modules = (Molecules,) radii(m)
        @test_opt target_modules = (Molecules,) vols(m)
        @test_opt target_modules = (Molecules,) r_max(m)
    end
end

@testset "JET: SASA's per-atom loop is free of runtime dispatch" begin
    # `KDTree(crds)` cannot infer to a concrete type -- the point dimension is a
    # runtime property of the Matrix. `sasa_atoms` therefore crosses a FUNCTION BARRIER
    # into `_sasa_loop!`, which Julia specialises on the concrete tree type.
    #
    # The barrier call is itself one dynamic dispatch, but exactly one per
    # `sasa_atoms` call rather than one `inrange` dispatch per atom (which is
    # what this used to be, and was worth ~21% of runtime).
    _reports(f, types) =
        JET.get_reports(JET.report_opt(f, types; target_modules = (SASA,)))

    m = create("t", ["o", "h", "h"],
               [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)])
    crds = Molecules.coords_cartesian(m)
    TT   = typeof(SASA.KDTree(crds))

    # the loop that runs once per atom must be completely clean
    @test isempty(_reports(SASA._sasa_loop!,
        (Vector{Float64}, TT, Matrix{Float64}, Vector{Float64}, Float64,
         Vector{SASA.PlasticMap.Vec3}, Float64, Int, Int, Float64)))

    # the entry points carry exactly the one barrier dispatch, and it is the
    # barrier -- not `inrange`, and not anything inside the per-atom loop.
    for f in (SASA.sasa_atoms, SASA.sasa)
        rs = _reports(f, (Molecule,))
        @test length(rs) <= 1
        @test all(r -> occursin("_sasa_loop!", sprint(show, r)), rs)
    end
end

