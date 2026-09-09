"""
Swappable-backend interfaces. Declares the backend markers (`RadiiSource`,
`FormFactorSource`) and their generic functions then encapsulates the two 
bundled backends as child submodules:

    -   `AtomicRadii` — atomic/ionic radii from a bundled SQLite file
        (`AtomicRadiiSource <: RadiiSource`, extends [`lookup`](@ref)).

    -   `FormFactor` — X-ray form factors from Waasmaier-Kirfel and Chantler tables
        (`FormFactorSourceTables <: FormFactorSource`, extends
        [`form_factor_table`](@ref) / [`form_factors`](@ref) /
        [`form_factor_log`](@ref)).

Non-test code elsewhere in the package reaches both backends ONLY through
this module (`Interfaces.<name>`); the submodule internals are not
re-exported.
"""
module Interfaces

export  RadiiSource, FormFactorSource, AtomicRadiiSource, FormFactorSourceTables,
        FormFactorError, lookup, form_factor_table, form_factors, form_factor_log

"A source of atomic/ionic radii. Implement [`lookup`](@ref) for a concrete subtype."
abstract type RadiiSource end

"""
    lookup(src::RadiiSource, ions) -> Vector{Tuple{String,Union{Float64,Nothing}}}

Resolve each ion/element string to a radius in Å, or `nothing` if unknown. One
entry per input, in input order.

# Arguments
- `src`: the radii backend to query.
- `ions`: ion/element strings, e.g. `["fe3+", "o2-", "fe"]`.
"""
function lookup end

"A form-factor backend. See `FormFactor` for the reference implementation."
abstract type FormFactorSource end

"""
    form_factor_table([src::FormFactorSource,] energy::Real, ions, qvals) -> FF

Build a form-factor container for `ions` at one `energy` over the `qvals` grid.

# Arguments
- `src`: the form-factor backend to query (optional).
- `energy`: photon energy in eV.
- `ions`: vector of ion strings.
- `qvals`: vector of q values in Å⁻¹.
"""
function form_factor_table end

"""
    form_factors(t, ions, qvals) -> Matrix{ComplexF64}

Per-ion form-factor rows from a container `t` previously built by
[`form_factor_table`](@ref), as a `(length(ions), length(qvals))` matrix: row
`i` is the row for `ions[i]`, columns aligned to `qvals` in input order. Pass
the per-atom ion vector and the result is exactly the `f_atoms` matrix
[`compute_B_lm`](@ref) takes — no mapping step in between.

The queried `qvals` must be grid points of `t`, and every `ions[i]` must be
present in `t`; either violation throws [`FormFactorError`](@ref).

# Arguments
    - `t`: form-factor container from [`form_factor_table`](@ref).
    - `ions`: ion strings to fetch, one per output row.
    - `qvals`: vector of q values; each must match a grid point of `t` exactly.
"""
function form_factors end

"""
    form_factor_log(t) -> Vector{String}

Construction-time diagnostics for the container `t` built by
[`form_factor_table`](@ref) — one line per ion the backend could not resolve
in full, in the order they were encountered.

# Arguments
- `t`: form-factor container from [`form_factor_table`](@ref).
"""
function form_factor_log end

include("AtomicRadii/AtomicRadii.jl")
include("FormFactor/FormFactor.jl")

using .AtomicRadii: AtomicRadii, AtomicRadiiSource
using .FormFactor: FormFactor, FormFactorSourceTables, FormFactorError

end # module
