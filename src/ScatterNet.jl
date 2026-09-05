"""
Thin aggregator. Two grouped submodules — `Molecule` (`Molecules`,
`SASA.PlasticMap`) and `Scattering` (`SphFuncs`, `FormFactorXrayDB`) — over the
shared `Interfaces` markers and `Constants` (numeric tolerances).
"""
module ScatterNet

include("Constants.jl")
include("Interfaces.jl")
include("Molecule/Molecule.jl")
include("Scattering/Scattering.jl")

using .Constants: Constants
using .Interfaces: Interfaces
using .Molecule: Molecule
using .Scattering: Scattering

end # module
