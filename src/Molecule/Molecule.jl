module Molecule

include("Molecules.jl")
include("SASA.jl")

using .Molecules: Molecules
using .SASA: SASA

end # module
