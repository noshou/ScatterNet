"""
Thin aggregator. Two grouped submodules — `Molecule` (`AtomicRadii`,
`Molecules`) and `Scattering` (`SphFuncs`, `FormFactorXrayDB`) — over the shared
`Interfaces` markers.
"""
module ScatterNet

import CondaPkg  # keep as a direct dep for the root CondaPkg.toml

include("Interfaces.jl")
include("Molecule/Molecule.jl")
include("Scattering/Scattering.jl")

using .Interfaces: Interfaces
using .Molecule: Molecule
using .Scattering: Scattering

end # module
