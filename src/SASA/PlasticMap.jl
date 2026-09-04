"""
Even, incrementally-extensible point sets on the unit sphere, drawn from the 2-D
plastic (R₂) low-discrepancy sequence.

Sequence term `i` is a fixed function of `i` alone, so the first `k` points are a
byte-exact prefix of the first `k′ > k`. Grow a sampling in place with
[`plastic_points!`](@ref); any per-point work already cached over the shared
prefix stays valid, so a molecule-wide buffer can be extended for the few atoms
that need more resolution without disturbing the rest.

Points are stored as [`Vec3`](@ref) unit vectors: the trig runs once at
generation, and mapping onto an atom is `center .+ (r + probe) .* p`.

Not synchronized. Mutating calls (`plastic_points!`, and `plastic_points` via
it) must happen before any parallel fan-out; treat the buffer as read-only for
the duration of a threaded phase.
"""
module PlasticMap

using Roots: find_zero

export Vec3, plastic_points, plastic_points!, PLASTIC_RATIO

"A unit vector on the sphere, `(x, y, z)`. `isbits`, so `Vector{Vec3}` is contiguous and cheap to grow."
const Vec3 = NTuple{3,Float64}

"Cubic `x³ − x − 1`; its real root in `(1, 2)` is the plastic ratio."
_plastic_poly(x) = x^3 - x - 1

"Plastic ratio ρ ≈ 1.324718, the real root of `x³ = x + 1`."
const PLASTIC_RATIO     = find_zero(_plastic_poly, (1.0, 2.0))
const PLASTIC_RATIO_SQR = PLASTIC_RATIO^2

"Fractional part of `x`, i.e. `x - floor(x)`, in `[0, 1)`."
_frac(x) = x - floor(x)

"""
    _plastic_point(i::Int) -> Vec3

Unit-sphere point for 1-based plastic-sequence term `i`. The 2-D term
`(frac(i/ρ), frac(i/ρ²))` is read as (azimuth, height) and lifted to the sphere
through the equal-area cylindrical projection, so the points are uniform in area
rather than clustered at the poles.

Division by `ρ`/`ρ²` (rather than multiplication by reciprocals) keeps the
fractional part accurate; it degrades only once `i` nears the mantissa limit
(~1e15), far above any SASA point count.
"""
@inline function _plastic_point(i::Int)::Vec3
    φ = 2.0 * pi * _frac(i / PLASTIC_RATIO)
    z = 2.0 * _frac(i / PLASTIC_RATIO_SQR) - 1.0
    r = sqrt(max(0.0, 1.0 - z * z))
    return (r * cos(φ), r * sin(φ), z)
end

"""
    plastic_points!(buf::Vector{Vec3}, have::Int, want::Int) -> buf

Extend `buf` in place to hold the first `want` points of the plastic sequence,
given that it already holds the first `have`. Appends `want - have` points and
leaves the existing entries untouched; a no-op when `want <= have`.

# Arguments
- `buf`: point buffer to grow; its length must equal `have` on entry.
- `have`: count of valid leading entries already in `buf`; `have >= 0`.
- `want`: desired total point count; `want >= 0`.
"""
function plastic_points!(buf::Vector{Vec3}, have::Int, want::Int)
    have < 0 && throw(DomainError(have, "have must be >= 0"))
    want < 0 && throw(DomainError(want, "want must be >= 0"))
    length(buf) == have ||
        throw(ArgumentError("buf has $(length(buf)) entries but have = $have"))
    want <= have && return buf

    sizehint!(buf, want)
    @inbounds for i in (have + 1):want
        push!(buf, _plastic_point(i))
    end
    return buf
end

"""
    plastic_points(n::Int) -> Vector{Vec3}

The first `n` plastic-sequence points, freshly allocated. Convenience wrapper
over [`plastic_points!`](@ref) for callers that don't reuse a buffer.
"""
plastic_points(n::Int) = plastic_points!(Vector{Vec3}(), 0, n)

end # module
