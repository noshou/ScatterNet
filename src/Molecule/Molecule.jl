"""
Real-space structure: atomic/ionic radii and the molecule geometry built on
them. Depends on the top-level `Interfaces` markers; `Cache` is local here.
"""
module Molecule

include("Cache.jl")
include("AtomicRadii.jl")
include("Molecules.jl")

using .Cache: Cache
using .AtomicRadii: AtomicRadii
using .Molecules: Molecules

end # module
