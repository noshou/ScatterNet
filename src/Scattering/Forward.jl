# The forward model, assembled. This is the module-level entry point: given a
# molecule, a q grid, a band limit and the beam energy, build the five species'
# multipoles, reduce them to the species Gram matrix `G(q)`, and contract that
# against a contrast vector to a detector-scale intensity.
#
#     mol ─► (B_vac, B_ex, B_sh_convex, B_sh_concave, B_sh_cavity)   species_multipoles
#         ─► G(q) ∈ ℝ^{5×5×Q}                                        gram_matrix
#         ─► ForwardCache(G, qvals, r_m)                             forward_cache
#         ─► I_calc(q) = m·(v(q)ᵀ G v(q)) + c                        forward
#
# `G(q)` depends only on geometry and beam, never on `(m, c, dns, ρ, r0)`. Build the
# cache once per structure with `forward_cache`; every likelihood evaluation is then
# the O(Q) `forward(cache, m, c, dns, ρ; r0)`.

"""
    species_multipoles(mol, qvals, lMax, energy; kwargs...)
        -> NTuple{5, Array{ComplexF64,3}}

The five CRYSOL-3 species' multipoles `B_lm`, in the fixed order

    (vac, ex, sh_convex, sh_concave, sh_cavity)

that [`gram`](@ref) and [`contrast_vector`](@ref) assume. Each is `(C, K, Q)` in
[`compute_B_lm`](@ref)'s packed layout (`C = 2` for `vac` near an absorption
edge, `1` otherwise; `C = 1` for the four dummy species). Geometry and beam
only -- no contrast parameters enter here.

# Arguments
- `mol::Molecule`.
- `qvals::AbstractVector{<:Real}`, length `Q`: momentum transfer in Å⁻¹.
- `lMax::Integer`: spherical-harmonic band limit.
- `energy::Real`: photon energy in eV (for the `vac` anomalous form factors).

# Keywords
    - `ions::Vector{String} = elms(mol)`: ion/element string per atom.
    - `chunk::Unsigned = B_LM_CHUNK`: `compute_B_lm` batch size (results invariant).
    - `form_factor_source::FormFactorSource = FORM_FACTOR_SOURCE`.
    - `thickness::Real = SHELL_THICKNESS`, `probe::Real = PROBE_RADIUS`,
    `n_target = SHELL_N_TARGET`, `classes = SHELL_CLASSES`: forwarded to
    [`hydration`](@ref). A class not in `classes` comes back an all-zero `B_lm`
    (≡ its `dro_k = 0`), so the return is always a 5-tuple.
"""
function species_multipoles(
    mol::Molecule, qvals::AbstractVector{<:Real}, lMax::Integer, energy::Real;
    ions::Vector{String}                 = Molecules.elms(mol),
    chunk::Unsigned                      = B_LM_CHUNK,
    form_factor_source::FormFactorSource = FORM_FACTOR_SOURCE,
    thickness::Real                      = SHELL_THICKNESS,
    probe::Real                          = PROBE_RADIUS,
    n_target::Union{Nothing,Integer}     = SHELL_N_TARGET,
    classes                              = SHELL_CLASSES,
)
    _CHUNK = UInt64(chunk)
    b_vac = vacuo(
        mol, 
        qvals, 
        lMax, 
        ions, 
        Float64(energy), 
        _CHUNK;
        form_factor_source = form_factor_source
    )
    b_ex  = excluded(mol, qvals, lMax, _CHUNK)
    sh    = hydration(mol, qvals, lMax, _CHUNK;
                    thickness = Float64(thickness), probe = Float64(probe),
                    n_target  = n_target, classes = classes)
    return (b_vac, b_ex, sh.convex, sh.concave, sh.cavity)
end

"""
    gram_matrix(mol, qvals, lMax, energy; kwargs...) -> Array{Float64,3}

The five-species Gram matrix `G_ab(q) = S_ab(q)`, shape `(5, 5, Q)`, symmetric
and positive-semidefinite in its species axes. Depends only on geometry and beam.

`kwargs` are those of [`species_multipoles`](@ref).
"""
gram_matrix(mol::Molecule, qvals::AbstractVector{<:Real}, lMax::Integer, energy::Real;
            kwargs...)::Array{Float64,3} =
    gram(collect(species_multipoles(mol, qvals, lMax, energy; kwargs...)),
        partial_wave_weights(lMax))

"""
    mean_atomic_radius(mol) -> Float64

The structure's mean atomic radius `r_m` in Å, which is the reference radius CRYSOL's
fitted `r₀` is measured against (`c₁ = r₀ / r_m`), and the point at which
[`excluded_volume_factor`](@ref)'s single-envelope approximation is exact.
Geometry only, so it caches alongside `G`.
"""
function mean_atomic_radius(mol::Molecule)::Float64
    r = Molecules.radii(mol)
    isempty(r) && throw(ArgumentError("mean_atomic_radius: molecule has no atoms"))
    return sum(r) / length(r)
end

"""
    ForwardCache

Everything about a structure the forward model needs that does **not** depend on
the fit parameters: the species Gram matrix `G`, the `q` grid it was built on,
and the mean atomic radius `r_m` that [`excluded_volume_factor`](@ref) measures
`r₀` against. Build one per structure with [`forward_cache`](@ref); every
likelihood evaluation is then the O(Q) [`forward`](@ref)`(cache, …)`.

# Fields
- `G::Array{Float64,3}`, `(n, n, Q)`: from [`gram_matrix`](@ref).
- `qvals::Vector{Float64}`, length `Q`: the grid `G` was built on.
- `r_m::Float64`: mean atomic radius in Å.
"""
struct ForwardCache
    G::Array{Float64,3}
    qvals::Vector{Float64}
    r_m::Float64
end

"""
    forward_cache(mol, qvals, lMax, energy; kwargs...) -> ForwardCache

The geometry-only pass: [`gram_matrix`](@ref) plus the `q` grid and
[`mean_atomic_radius`](@ref) that fitting `r₀` needs. **Cache once per
structure**; `kwargs` are those of [`species_multipoles`](@ref).
"""
forward_cache(
    mol::Molecule, 
    qvals::AbstractVector{<:Real}, 
    lMax::Integer,
    energy::Real; 
    kwargs...
)::ForwardCache =
    ForwardCache(
        gram_matrix(mol, qvals, lMax, energy; kwargs...),
        collect(Float64, qvals), 
        mean_atomic_radius(mol)
    )

"""
    forward(G, m, c, dns, ρ) -> Vector

The detector-scale model intensity `I_calc(q) = m · (vᵀ G(q) v) + c`, with the
contrast vector `v = contrast_vector(dns, ρ)`. `G` is from [`gram_matrix`](@ref)
(or [`gram`](@ref)); `ρ` is `(ρ_convex, ρ_concave, ρ_cavity)` for a five-species
`G` or a scalar for a three-species one.

This method holds the excluded-volume radius at its tabulated value. To fit
CRYSOL's `r₀` as well, use the [`ForwardCache`](@ref) method below.
"""
forward(G::AbstractArray{<:Real,3}, m::Real, c::Real, dns::Real, ρ) =
    intensity_calc(intensity(G, contrast_vector(dns, ρ)), m, c)

"""
    forward(cache, m, c, dns, ρ; r0 = nothing) -> Vector

The full CRYSOL-parameter forward model:

    I_calc(q) = m · (v(q)ᵀ G(q) v(q)) + c
    v(q)      = [1, -dns · G_ex(q; r₀), dro₁, dro₂, dro₃]

`r0` is CRYSOL's fitted excluded-volume radius in Å, and `dro_k = DRO_UNIT · ρ_k`
its fitted shell contrasts. `r0 = nothing` (the default) or `r0 == cache.r_m`
holds the excluded volume at its tabulated value and reduces exactly to the
5-argument [`forward`](@ref) above.

O(Q) in the fit parameters -- `cache.G` is never rebuilt, because `r₀` enters as
a `q`-dependent reweighting of the contrast rather than of the geometry (see
[`excluded_volume_factor`](@ref)). The whole path promotes its element type, so
`ForwardDiff`/`Zygote` can differentiate `(m, c, dns, ρ, r0)` through it.

# A note on `dns` and `r₀`

CRYSOL fixes the bulk solvent density at `dns = 0.334 e·Å⁻³` and fits `r₀`.
Both are exposed here, but they are strongly degenerate -- `dns` scales the `ex`
species and `r₀` scales it by `c₁³` at `q = 0` -- so fitting both leaves a
near-flat direction. Fix one (CRYSOL's choice: `dns`) unless the sampler is
explicitly meant to explore that ridge.
"""
function forward(cache::ForwardCache, m::Real, c::Real, dns::Real, ρ;
                 r0::Union{Nothing,Real} = nothing)
    (r0 === nothing || r0 == cache.r_m) &&
        return forward(cache.G, m, c, dns, ρ)
    g_ex = excluded_volume_factor(cache.qvals, cache.r_m, r0)
    return intensity_calc(intensity(cache.G, contrast_matrix(dns, ρ, g_ex)), m, c)
end

"""
    forward(mol, qvals, lMax, energy; m, c, dns, ρ, r0, kwargs...) -> Vector

The whole forward model from a molecule in one call: build the cache, then
evaluate. A convenience for one-offs; for repeated evaluation at fixed geometry
call [`forward_cache`](@ref) once and the [`ForwardCache`](@ref) method of
[`forward`](@ref) per parameter set.

`kwargs` beyond `m, c, dns, ρ, r0` are those of [`species_multipoles`](@ref).
"""
function forward(
    mol::Molecule, qvals::AbstractVector{<:Real}, lMax::Integer, energy::Real;
    m::Real, c::Real, dns::Real, ρ, r0::Union{Nothing,Real} = nothing, kwargs...
)
    return forward(forward_cache(mol, qvals, lMax, energy; kwargs...),
                   m, c, dns, ρ; r0 = r0)
end
