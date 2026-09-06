using Test
# Python airlock: nothing at the top level of this suite loads PythonCall, so a
# default `Pkg.test()` never provisions the CondaPkg Python env. The xraydb
# form-factor tests in `test_formfactor.jl` are opt-in -- run them with
# `SCATTERNET_TEST_XRAYDB=1`, which is also the only path that loads PythonCall
# and lets CondaPkg build the numpy + xraydb env.
using ScatterNet
using ScatterNet: Interfaces
using ScatterNet.ABSOLUTE_TOLERANCE: DEFAULT_ATOL
using ScatterNet.Interfaces: AtomicRadii
using ScatterNet.Molecule: Molecules
using ScatterNet.Molecule.SASA: PlasticMap
using ScatterNet.Scattering: SphFuncs
using ScatterNet.Interfaces: FormFactorXrayDB

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
