"Shared numeric constants used across the forward model and its tests."
module Constants

using Roots: find_zero

export DEFAULT_ATOL, PLASTIC_RATIO

"""
Default absolute tolerance for floating-point equality checks
(`abs(a - b) < DEFAULT_ATOL`, or `isapprox(a, b; atol = DEFAULT_ATOL)`):
a few orders of magnitude above `Float64` roundoff accumulated over a short
chain of elementary operations.
"""
const DEFAULT_ATOL = 1.0e-9

"Cubic `x³ − x − 1`; its real root in `(1, 2)` is the plastic ratio."
_plastic_poly(x) = x^3 - x - 1

"Plastic ratio ρ ≈ 1.324718, the real root of `x³ = x + 1`."
const PLASTIC_RATIO     = find_zero(_plastic_poly, (1.0, 2.0))
const PLASTIC_RATIO_SQR = PLASTIC_RATIO^2

"Fractional part of `x`, i.e. `x - floor(x)`, in `[0, 1)`."
_frac(x) = x - floor(x)

end # module
