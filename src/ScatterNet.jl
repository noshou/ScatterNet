"""
Thin aggregator. `Interfaces` is the swappable-backend facade: it declares the
`RadiiSource`/`FormFactorSource` markers plus their generic functions and
encapsulates the bundled backends as its own submodules.
"""
module ScatterNet

include("ABSOLUTE_TOLERANCE.jl")
include("Interfaces/Interfaces.jl")
include("Molecule/Molecule.jl")
include("Scattering/Scattering.jl")

using .ABSOLUTE_TOLERANCE: ABSOLUTE_TOLERANCE
using .Interfaces: Interfaces
using .Molecule: Molecule
using .Scattering: Scattering

end # module
