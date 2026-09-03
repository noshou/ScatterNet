"""
ScatterNet — forward model for X-ray/neutron scattering shape reconstruction.
Submodules: `Cache`, `Interfaces`, `AtomicRadii`, `SphFuncs`, `Molecules`, `FormFactorXrayDB`.
"""
module ScatterNet

import CondaPkg  # keep as a direct dep for the root CondaPkg.toml

include("Cache.jl")
include("Interfaces.jl")
include("AtomicRadii.jl")
include("SphFuncs.jl")
include("Molecules.jl")
include("FormFactorXrayDB.jl")

using .Cache: Cache
using .Interfaces: Interfaces, RadiiSource, FormFactorSource
using .AtomicRadii: AtomicRadii, Ion, tryparse_ion, resolve_one, AtomicRadiiSource
using .SphFuncs: SphFuncs, sphHarm, sphBess, legendre_sphPlm, SphHarmError, SphBessError
using .Molecules: Molecules, Molecule, MoleculeError
using .FormFactorXrayDB: FormFactorXrayDB, FF, FormFactorError,
                         compute_form_factors, FormFactorSourceXrayDB

export Molecule, MoleculeError, Ion, tryparse_ion, resolve_one, AtomicRadiiSource,
       sphHarm, sphBess, legendre_sphPlm, SphHarmError, SphBessError,
       compute_form_factors, FormFactorSourceXrayDB, FF, FormFactorError

end # module
