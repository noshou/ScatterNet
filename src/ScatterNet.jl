"""
Thin aggregator. Two grouped submodules — `Structure` (`AtomicRadii`,
`Molecules`) and `Scattering` (`SphFuncs`, `FormFactorXrayDB`) — over the shared
`Interfaces` markers.
"""
module ScatterNet

import CondaPkg  # keep as a direct dep for the root CondaPkg.toml

include("Interfaces.jl")
include("Structure/Structure.jl")
include("Scattering/Scattering.jl")

using .Interfaces: Interfaces
using .Structure: Structure
using .Scattering: Scattering

end # module
