"Complex spherical harmonics and spherical Bessel functions over point/grid inputs."
module SphFuncs

using SphericalHarmonics: SphericalHarmonics
using LegendrePolynomials: Plm
using Bessels: sphericalbesselj

export legendre_sphPlm, sphHarm, sphBess, SphHarmError, SphBessError

"Raised by [`sphHarm`](@ref) on invalid input."
struct SphHarmError <: Exception; msg::String end
"Raised by [`sphBess`](@ref) on invalid input."
struct SphBessError <: Exception; msg::String end
Base.showerror(io::IO, e::SphHarmError) = print(io, "SphHarmError: ", e.msg)
Base.showerror(io::IO, e::SphBessError) = print(io, "SphBessError: ", e.msg)

const _INV_SQRT_2PI = 1.0 / sqrt(2.0 * pi)

"""
    legendre_sphPlm(l::Integer, m::Integer, x::Real) -> Float64

Normalized associated Legendre P̄_l^m(x) with Condon–Shortley phase (GSL `legendre_sphPlm`).

# Arguments
- `l`: degree, `l >= 0`.
- `m`: order, `0 <= m <= l`.
- `x`: argument, typically `cos θ` in `[-1, 1]`.
"""
@inline legendre_sphPlm(l::Integer, m::Integer, x::Real)::Float64 =
    Plm(float(x), l, m; norm = Val(:normalized), csphase = true) * _INV_SQRT_2PI

"""
    sphHarm(lMax::Int, θ::AbstractArray{<:Real}, φ::AbstractArray{<:Real}) -> Matrix{ComplexF64}

Complex spherical harmonics Y_l^m for `l = 0..lMax`, `m = 0..l`. Rows are packed as
`l*(l+1)÷2 + m + 1`; columns are input points.

# Arguments
- `lMax`: maximum degree, `lMax >= 0`.
- `θ`: 1-D vector of polar angles; same length as `φ`.
- `φ`: 1-D vector of azimuthal angles; same length as `θ`.
"""
function sphHarm(lMax::Int, θ::AbstractArray{<:Real}, φ::AbstractArray{<:Real})::Matrix{ComplexF64}
    lMax < 0 && throw(SphHarmError("lMax must be >= 0"))
    (ndims(θ) == 1 && ndims(φ) == 1) || throw(SphHarmError("theta/phi must be 1-D"))
    (isempty(θ) || isempty(φ)) && throw(SphHarmError("theta/phi must be non-empty"))
    length(θ) == length(φ) || throw(SphHarmError("theta/phi length mismatch"))

    npts = length(θ)
    y = Matrix{ComplexF64}(undef, (lMax + 1) * (lMax + 2) ÷ 2, npts)
    S = SphericalHarmonics.cache(Int(lMax))
    θv, φv = vec(θ), vec(φ)
    @inbounds for i in 1:npts
        θi = Float64(θv[i])
        SphericalHarmonics.computePlmcostheta!(S, θi, lMax)
        Yi = SphericalHarmonics.computeYlm!(S, θi, Float64(φv[i]), lMax)
        for l in 0:lMax, m in 0:l
            y[l * (l + 1) ÷ 2 + m + 1, i] = Yi[(l, m)]
        end
    end
    return y
end

"""
    sphBess(r::AbstractArray{<:Real}, q::AbstractArray{<:Real}, lMax::Int) -> Array{Float64,3}

Spherical Bessel functions j_l for `l = 0..lMax` over the outer product `q ⊗ r`,
shape `(lMax+1, |q|, |r|)`.

# Arguments
- `r`: non-empty vector of radii, all `>= 0`.
- `q`: non-empty vector of q values, all `>= 0`.
- `lMax`: maximum order, `lMax >= 0`.
"""
function sphBess(r::AbstractArray{<:Real}, q::AbstractArray{<:Real}, lMax::Int)::Array{Float64,3}
    isempty(r) && throw(SphBessError("radii must be non-empty"))
    isempty(q) && throw(SphBessError("q grid must be non-empty"))
    any(<(0), r) && throw(SphBessError("radii must be >= 0"))
    any(<(0), q) && throw(SphBessError("q must be >= 0"))
    lMax < 0 && throw(SphBessError("lMax must be >= 0"))

    rv, qv = vec(r), vec(q)
    j = Array{Float64,3}(undef, lMax + 1, length(qv), length(rv))
    @inbounds for qi in eachindex(qv), ri in eachindex(rv)
        x = Float64(qv[qi]) * Float64(rv[ri])
        for l in 0:lMax
            j[l + 1, qi, ri] = sphericalbesselj(l, x)
        end
    end
    return j
end

end # module
