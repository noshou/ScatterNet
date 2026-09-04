"""
Scattering primitives: spherical harmonics / Bessel functions and X-ray form
factors. Depends on the top-level `Interfaces` markers.
"""
module Scattering

include("SphFuncs.jl")
include("FormFactorXrayDB.jl")

using .SphFuncs: SphFuncs
using .FormFactorXrayDB: FormFactorXrayDB

end # module
