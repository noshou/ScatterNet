"""
Real-space structure: the `Molecule` type and SASA point sampling. Depends on
the top-level `Interfaces` markers.

Lazy memoization and the default radii backend are implementation details
folded directly into `Molecules` (no `Cache`/`AtomicRadii` module of their
own -- see `Molecules.jl`'s docstring), since neither is used anywhere else.
A caller who wants a different radii source implements `Interfaces.RadiiSource`
and passes it to `Molecules.create`'s `radii_source` keyword; they never touch
that internal code. `Molecules` and `SASA` are the two supported submodules.
"""
module Molecule

include("Molecules.jl")
include("SASA.jl")

using .Molecules: Molecules
using .SASA: SASA

end # module
