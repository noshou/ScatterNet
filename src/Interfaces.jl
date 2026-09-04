"Swappable-backend markers: `RadiiSource` and `FormFactorSource`."
module Interfaces

export RadiiSource, lookup, FormFactorSource

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

"A form-factor backend. See `FormFactorXrayDB` for the reference implementation's `create`/`log`/`lookup`."
abstract type FormFactorSource end

end # module
