"""
Thin aggregator. Three grouped submodules — `Molecule` (`AtomicRadii`,
`Molecules`), `Scattering` (`SphFuncs`, `FormFactorXrayDB`) and `SASA`
(`PlasticMap`) — over the shared `Interfaces` markers.
"""
module ScatterNet

import CondaPkg  # keep as a direct dep for the root CondaPkg.toml

include("Interfaces.jl")
include("Molecule/Molecule.jl")
include("Scattering/Scattering.jl")
include("SASA/SASA.jl")

using .Interfaces: Interfaces
using .Molecule: Molecule
using .Scattering: Scattering
using .SASA: SASA

end # module
