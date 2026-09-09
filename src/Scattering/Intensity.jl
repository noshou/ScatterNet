# Folds the per-species multipoles `B_lm` # (from `vacuo` / `excluded` / `hydration`) 
# into the orientationally-averaged `I(q)` via the species Gram matrix.
#
# SPECIES ORDERING:
#
#     1  vac         real atoms in vacuo
#     2  ex          excluded-volume dummies      (enters A_total with -dns)
#     3  sh_convex   hydration shell, convex beads   (enters with +dro_1)
#     4  sh_concave  hydration shell, concave beads  (enters with +dro_2)
#     5  sh_cavity   hydration shell, cavity beads   (enters with +dro_3)

"""
    gram(Bs, weights) -> Array{Float64,3}

Species Gram matrix `G` of shape `(n, n, Q)` for the `n` multipole arrays in
`Bs` (each `(C, K, Q)` as returned by [`compute_B_lm`](@ref); channel counts
`C` may differ between species).

`G[a, b, :]` is the paired partial-wave sum

    S_ab(q) = 4π Σ_c Σ_lm w_lm Re(B_a[c,lm](q) * conj(B_b[c,lm](q)))

i.e. [`cross_scatter`](@ref) off the diagonal and [`self_scatter`](@ref) on it
(`S_aa = S_ab|_{b=a}` identically, since `|z|² = Re(z * conj z)`). Every
`G[:, :, k]` is real, symmetric and positive semidefinite.

# Arguments
    - `Bs`: an `AbstractVector` (or `Tuple`) of `n ≥ 1` arrays
    `<:AbstractArray{<:Complex,3}`, each `(C, K, Q)` sharing a common `K` and `Q`.
    Order is (vac, ex, shells…); see the file header.
    - `weights::AbstractVector{<:Real}`, length `K`:
    [`partial_wave_weights`](@ref)`(lMax)`.

# Returns
- `Array{Float64,3}` of size `(n, n, Q)`, symmetric in its first two axes.
"""
function gram(
    Bs::AbstractVector{<:AbstractArray{<:Complex,3}}, weights::AbstractVector{<:Real}
)::Array{Float64,3}
    n = length(Bs)
    n ≥ 1 || throw(ArgumentError("gram: need at least one species"))
    K = length(weights)
    Q = size(Bs[1], 3)
    for a in 1:n
        size(Bs[a], 2) == K || throw(ArgumentError(
            "gram: Bs[$a] has $(size(Bs[a], 2)) (l,m) rows but weights has length $K"))
        size(Bs[a], 3) == Q || throw(ArgumentError(
            "gram: Bs[$a] has Q=$(size(Bs[a], 3)) but Bs[1] has Q=$Q"))
    end

    G = Array{Float64,3}(undef, n, n, Q)
    @inbounds for a in 1:n
        G[a, a, :] .= self_scatter(Bs[a], weights)
        for b in (a + 1):n
            s = cross_scatter(Bs[a], Bs[b], weights)
            G[a, b, :] .= s
            G[b, a, :] .= s
        end
    end
    return G
end

gram(Bs::Tuple, weights::AbstractVector{<:Real}) = gram(collect(Bs), weights)

"""
    intensity(G, v) -> Vector
    intensity(G, V) -> Vector

The orientationally-averaged absolute model intensity

    I(q) = v' * G(:, :, q) * v = Σ_a Σ_b v_a v_b G[a, b, q]

the contrast-weighted contraction of the species Gram matrix from [`gram`](@ref).
`v` is the contrast vector `[1, -dns, dro…]` (see [`contrast_vector`](@ref)).
`I(q) ≥ 0` because every `G(:, :, q)` is positive semidefinite.

The second method takes a **q-dependent** contrast `V::(n, Q)`, whose column `k`
is the contrast vector at `qvals[k]`.

# Arguments
    -   `G::AbstractArray{<:Real,3}`, `(n, n, Q)`: species Gram matrix.
    -   `v::AbstractVector{<:Real}`, length `n`: species contrast coefficients; or
        `V::AbstractMatrix{<:Real}`, `(n, Q)`: per-`q` contrast coefficients.

# Returns
- `Vector` of length `Q`, `eltype` promoted from `G` and the contrast.
"""
function intensity(G::AbstractArray{<:Real,3}, v::AbstractVector{<:Real})
    n = length(v)
    (size(G, 1) == n && size(G, 2) == n) || throw(DimensionMismatch(
        "intensity: G is $(size(G, 1))×$(size(G, 2)) in its species axes but v has length $n"))
    Q = size(G, 3)
    T = promote_type(eltype(G), eltype(v))
    out = Vector{T}(undef, Q)
    @inbounds for k in 1:Q
        acc = zero(T)
        for b in 1:n
            vb = v[b]
            for a in 1:n
                acc += v[a] * G[a, b, k] * vb
            end
        end
        out[k] = acc
    end
    return out
end

function intensity(G::AbstractArray{<:Real,3}, V::AbstractMatrix{<:Real})
    n = size(V, 1)
    (size(G, 1) == n && size(G, 2) == n) || throw(DimensionMismatch(
        "intensity: G is $(size(G, 1))×$(size(G, 2)) in its species axes but V has $n rows"))
    Q = size(G, 3)
    size(V, 2) == Q || throw(DimensionMismatch(
        "intensity: V has $(size(V, 2)) columns but G has Q=$Q"))
    T = promote_type(eltype(G), eltype(V))
    out = Vector{T}(undef, Q)
    @inbounds for k in 1:Q
        acc = zero(T)
        for b in 1:n
            vb = V[b, k]
            for a in 1:n
                acc += V[a, k] * G[a, b, k] * vb
            end
        end
        out[k] = acc
    end
    return out
end

"""
    intensity_calc(I, m, c) -> Vector{Float64}

Put the absolute model intensity `I` (from [`intensity`](@ref)) onto the
detector's scale: `I_calc(q) = m * I(q) + c`, with `m` the overall scale
between calculated (electrons²) and measured (arbitrary-unit) intensity and `c`
a flat background left by imperfect buffer subtraction.

# Arguments
- `I::AbstractVector{<:Real}`, length `Q`: absolute intensity.
- `m::Real`: overall scale.
- `c::Real`: constant background.

# Returns
- `Vector{Float64}` of length `Q`.
"""
intensity_calc(I::AbstractVector{<:Real}, m::Real, c::Real) = m .* I .+ c

"""
    contrast_vector(dns, ρ::NTuple{3,<:Real}) -> Vector{Float64}
    contrast_vector(dns, ρ::Real)             -> Vector{Float64}

The species contrast vector `v` that [`intensity`](@ref) contracts a [`gram`](@ref)
against. Species order matches the file header:

    3-species:  v = [1, -dns, dro]                    dro   = DRO_UNIT * ρ
    5-species:  v = [1, -dns, dro_1, dro_2, dro_3]    dro_k = DRO_UNIT * ρ_k

`dns` rescales the excluded-volume term to the true mean electron density of
the displaced bulk solvent. To supply a shell contrast in raw e·Å⁻³ instead,
build `[1.0, -dns, dro]` directly rather than routing through `ρ`.

# Arguments
    - `dns::Real`:  excluded-volume density scale.
    - `ρ`:          dimensionless shell contrast(s).

# Returns
- `Vector{Float64}` of length 3 or 5.
"""
function contrast_vector(dns::Real, ρ::NTuple{3,<:Real})
    T = promote_type(typeof(dns), eltype(ρ), typeof(DRO_UNIT))
    return T[one(T), -dns, DRO_UNIT * ρ[1], DRO_UNIT * ρ[2], DRO_UNIT * ρ[3]]
end

function contrast_vector(dns::Real, ρ::Real)
    T = promote_type(typeof(dns), typeof(ρ), typeof(DRO_UNIT))
    return T[one(T), -dns, DRO_UNIT * ρ]
end

# ---------------------------------------------------------------------------
# CRYSOL's r₀: the fitted excluded-volume radius
# ---------------------------------------------------------------------------

"""
Exponent coefficient of [`excluded_volume_factor`](@ref): `(4π/3)^(2/3) / 4π`.

Converts CRYSOL's *radius* parameterisation into the *volume* parameterisation
[`_gaussian_dummy`](@ref) is written in, via `V = (4π/3) r³` (which is exactly
`Molecules.sphere_volume`, so `r_m` and the dummy volumes stay consistent).
"""
const _EV_EXP_COEFF = (4π / 3)^(2 / 3) / (4π)

"""
    excluded_volume_factor(qvals, r_m, r0) -> Vector

CRYSOL's excluded-volume envelope `G(q)`: the factor multiplying the `ex`
species when every dummy atom's radius is expanded from the structure's mean
`r_m` to the fitted `r₀`.

    c₁ = r₀ / r_m
    G(q) = c₁³ · exp(−q² (c₁² − 1) (4π/3)^(2/3) r_m² / 4π)

The approximation is exact for an atom of radius `r_m` and degrades with the
spread of radii about it. Absorbs systemic biases introduced by atomic radii table. 

`r₀ == r_m` returns exactly `1.0` at every `q`, i.e. the uncorrected model.

# Arguments
- `qvals::AbstractVector{<:Real}`, length `Q`: momentum transfer in Å⁻¹.
- `r_m::Real`: the structure's mean atomic radius in Å; `> 0`.
- `r0::Real`: the fitted excluded-volume radius in Å; `> 0`. CRYSOL's `r₀`.

# Returns
- `Vector` of length `Q`, `eltype` promoted from the arguments.
"""
function excluded_volume_factor(qvals::AbstractVector{<:Real}, r_m::Real, r0::Real)
    r_m > 0 || throw(DomainError(r_m, "excluded_volume_factor: r_m must be > 0"))
    r0 > 0 || throw(DomainError(r0, "excluded_volume_factor: r0 must be > 0"))
    c1 = r0 / r_m
    k = (c1^2 - 1) * _EV_EXP_COEFF * r_m^2
    return @. c1^3 * exp(-(qvals^2) * k)
end

"""
    contrast_matrix(dns, ρ, g_ex) -> Matrix

The `q`-dependent contrast that [`intensity`](@ref) contracts a [`gram`](@ref)
against once `r₀` is fitted: [`contrast_vector`](@ref)`(dns, ρ)` with the `ex`
entry scaled by the [`excluded_volume_factor`](@ref) envelope `g_ex`.

Column `k` is the contrast vector at `qvals[k]`:

    V[:, k] = [1, -dns · g_ex[k], dro₁, dro₂, dro₃]

Only species 2 (`ex`) carries the envelope .

# Arguments
- `dns::Real`: bulk solvent electron density in e·Å⁻³ (CRYSOL fixes `0.334`).
- `ρ`: dimensionless shell contrast(s), as [`contrast_vector`](@ref).
- `g_ex::AbstractVector{<:Real}`, length `Q`: from [`excluded_volume_factor`](@ref).

# Returns
- `Matrix` of size `(n, Q)`, `n` = 3 or 5.
"""
function contrast_matrix(dns::Real, ρ, g_ex::AbstractVector{<:Real})
    v = contrast_vector(dns, ρ)
    T = promote_type(eltype(v), eltype(g_ex))
    V = Matrix{T}(undef, length(v), length(g_ex))
    @inbounds for k in eachindex(g_ex), a in eachindex(v)
        V[a, k] = a == 2 ? v[a] * g_ex[k] : v[a]
    end
    return V
end
