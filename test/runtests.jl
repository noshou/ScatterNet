using Test
using PythonCall     # loads FormFactorXrayDBExt, the xraydb backend
using ScatterNet
using ScatterNet: Interfaces
using ScatterNet.Constants: DEFAULT_ATOL
using ScatterNet.Molecule: Molecules
using ScatterNet.Molecule.SASA: PlasticMap
using ScatterNet.Scattering: SphFuncs, FormFactorXrayDB

check_float(a, b; atol = DEFAULT_ATOL) = abs(a - b) < atol
check_complex(a, b; atol = DEFAULT_ATOL) = abs(a - b) < atol

@testset "ScatterNet" begin
    include("test_cache.jl")
    include("test_atomicradii.jl")
    include("test_molecules.jl")
    include("test_sphfuncs.jl")
    include("test_formfactor.jl")
    include("test_plasticmap.jl")
    include("test_sasa.jl")
    include("test_quality.jl")
end
