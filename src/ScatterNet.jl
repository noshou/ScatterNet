"""
Thin aggregator. `AtomicRadii` (a `RadiiSource` implementation) sits beside
two grouped submodules — `Molecule` (`Molecules`, `SASA.PlasticMap`) and
`Scattering` (`SphFuncs`, `FormFactorXrayDB`) — over the shared `Interfaces`
markers and `ABSOLUTE_TOLERANCE` (`DEFAULT_ATOL`, the one numeric constant
used outside its own module -- everything else this package needs, like the
plastic ratio, lives right next to the one thing that uses it). `AtomicRadii`
is independent of `Molecule`: it only implements `Interfaces.RadiiSource`,
the abstract backend `Molecule.Molecules.create` consults.
"""
module ScatterNet

include("ABSOLUTE_TOLERANCE.jl")
include("Interfaces.jl")
include("AtomicRadii/AtomicRadii.jl")
include("Molecule/Molecule.jl")
include("Scattering/Scattering.jl")

using .ABSOLUTE_TOLERANCE: ABSOLUTE_TOLERANCE
using .Interfaces: Interfaces
using .AtomicRadii: AtomicRadii
using .Molecule: Molecule
using .Scattering: Scattering

end # module
