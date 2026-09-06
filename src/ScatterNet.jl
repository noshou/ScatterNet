"""
Thin aggregator. `Interfaces` is the swappable-backend facade: it declares the
`RadiiSource`/`FormFactorSource` markers plus their generic functions and
encapsulates the two bundled backends as its own submodules
(`Interfaces.AtomicRadii`, `Interfaces.FormFactorXrayDB`) — every non-test
consumer reaches them only as `Interfaces.<feature>`. Beside it sit two grouped
submodules — `Molecule` (`Molecules`, `SASA.PlasticMap`) and `Scattering`
(`SphFuncs`, `PartialWave`) — over the shared `Interfaces` facade and
`ABSOLUTE_TOLERANCE` (`DEFAULT_ATOL`, the one numeric constant used outside its
own module -- everything else this package needs, like the plastic ratio, lives
right next to the one thing that uses it). The `AtomicRadii` backend is
independent of `Molecule`: it only implements `Interfaces.RadiiSource`, the
abstract backend `Molecule.Molecules.create` consults.
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
