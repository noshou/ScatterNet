"Swappable-backend markers: `RadiiSource` and `FormFactorSource`."
module Interfaces

export RadiiSource, lookup, FormFactorSource

"A source of atomic/ionic radii. Implement `lookup(::T, ions) -> Vector{Tuple{String,Union{Float64,Nothing}}}`."
abstract type RadiiSource end
function lookup end

"A form-factor backend. See `FormFactorXrayDB` for the reference implementation's `create`/`log`/`lookup`."
abstract type FormFactorSource end

end # module
