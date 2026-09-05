"The absolute tolerance used for floating-point equality checks across the forward model and its tests."
module ABSOLUTE_TOLERANCE

export DEFAULT_ATOL

"""
Default absolute tolerance for floating-point equality checks (`abs(a - b) < DEFAULT_ATOL`, 
or `isapprox(a, b; atol = DEFAULT_ATOL)`): a few orders of magnitude above `Float64` roundoff 
accumulated over a short chain of elementary operations.
"""
const DEFAULT_ATOL = 1.0e-9

end # module
