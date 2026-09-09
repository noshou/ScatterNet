# to run the full test suite: 
#   julia --project=test test/runtests.jl

using Test
using ForwardDiff
using ScatterNet
using ScatterNet: Interfaces
using ScatterNet.ABSOLUTE_TOLERANCE: DEFAULT_ATOL
using ScatterNet.Interfaces: AtomicRadii
using ScatterNet.Molecule: Molecules
using ScatterNet.Molecule.SASA: PlasticMap
using ScatterNet: Scattering
using ScatterNet.Scattering: SphFuncs
using ScatterNet.Interfaces: FormFactor

check_float(a, b; atol = DEFAULT_ATOL) = abs(a - b) < atol
check_complex(a, b; atol = DEFAULT_ATOL) = abs(a - b) < atol

@testset "ScatterNet" begin
    include("test_cache.jl")
    include("test_atomicradii.jl")
    include("test_molecules.jl")
    include("test_sphfuncs.jl")
    include("test_partialwave.jl")
    include("test_scatterers.jl")
    include("test_intensity.jl")
    include("test_forward.jl")
    include("test_formfactor.jl")
    include("test_plasticmap.jl")
    include("test_sasa.jl")
    include("test_quality.jl")
end
