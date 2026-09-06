using Test
# The xraydb form-factor backend is exercised on every run: PythonCall and
# CondaPkg are test-environment dependencies, so `Pkg.test()` provisions the
# numpy + xraydb conda env and `test_formfactor.jl` / the `vacuo` half of
# `test_scatterers.jl` run unconditionally.
#
# The package itself stays Python-free: PythonCall/CondaPkg are `[weakdeps]` of
# the top-level Project.toml, so `using ScatterNet` alone loads neither.
import PythonCall
using ScatterNet
using ScatterNet: Interfaces
using ScatterNet.ABSOLUTE_TOLERANCE: DEFAULT_ATOL
using ScatterNet.Interfaces: AtomicRadii
using ScatterNet.Molecule: Molecules
using ScatterNet.Molecule.SASA: PlasticMap
using ScatterNet: Scattering
using ScatterNet.Scattering: SphFuncs
using ScatterNet.Interfaces: FormFactorXrayDB

check_float(a, b; atol = DEFAULT_ATOL) = abs(a - b) < atol
check_complex(a, b; atol = DEFAULT_ATOL) = abs(a - b) < atol

@testset "ScatterNet" begin
    include("test_cache.jl")
    include("test_atomicradii.jl")
    include("test_molecules.jl")
    include("test_sphfuncs.jl")
    include("test_partialwave.jl")
    include("test_scatterers.jl")
    include("test_formfactor.jl")
    include("test_plasticmap.jl")
    include("test_sasa.jl")
    include("test_quality.jl")
end
