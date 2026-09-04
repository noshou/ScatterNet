using Test
using ScatterNet
using ScatterNet: Interfaces
using ScatterNet.Molecule: Cache, AtomicRadii, Molecules
using ScatterNet.Scattering: SphFuncs, FormFactorXrayDB

check_float(a, b) = abs(a - b) < 1e-9
check_complex(a, b) = abs(a - b) < 1e-9

@testset "ScatterNet" begin
    include("test_cache.jl")
    include("test_atomicradii.jl")
    include("test_molecules.jl")
    include("test_sphfuncs.jl")
    include("test_formfactor.jl")
    include("test_quality.jl")
end
