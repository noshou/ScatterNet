# Per-species partial-wave terms: one function per scatterer species in
# `A_total(q) = A_vac(q) - dns*A_ex(q) + Σ_k dro_k*A_sh_k(q)` (see the
# `Scattering` docstring). Each builds that species' `(N, Q)` amplitude, hands it
# to the shared `compute_B_lm`, and returns just that species' multipoles. The
# `S_ab` reduction -- diagonals and cross terms alike -- is assembled downstream
# from the per-species `B_lm`, so nothing here calls `self_scatter`.
using ..Interfaces: Interfaces, FormFactorSource
using ..Molecule.SASA: SASA
using ..Molecule.Molecules: Molecules, Molecule

# `SHELL_THICKNESS`, `PROBE_RADIUS`, `SHELL_N_TARGET`, `SHELL_CLASSES` and
# `FORM_FACTOR_SOURCE` are defined at the `Scattering` module level (see
# `Scattering.jl`); this file only reads them as call defaults.

"""
    _gaussian_dummy(vols, qvals) -> Matrix{Float64}

Per-dummy amplitude `f[i,k] = v_i * exp(-q_k² * v_i^(2/3) / 4π)`, the
Fraser/MacRae/Suzuki (CRYSOL) form used by both dummy species.

A uniform sphere of volume `v_i`, approximated by the Gaussian of equal volume:
at `q -> 0` it scatters as `v_i`, decaying as a Gaussian whose width is set by
`v_i^(1/3)`.

# Arguments
- `vols::AbstractVector{<:Real}`, length `N`: per-dummy volume in Å³, `>= 0`.
- `qvals::AbstractVector{<:Real}`, length `Q`: momentum-transfer grid in Å⁻¹.

# Returns
-   `Matrix{Float64}`, `(N, Q)`, in [`compute_B_lm`](@ref)'s `f_atoms` layout.
"""
function _gaussian_dummy(
    vols::AbstractVector{<:Real}, qvals::AbstractVector{<:Real}
)::Matrix{Float64}
    any(<(0), vols) && throw(ArgumentError("_gaussian_dummy: volumes must be >= 0"))
    # (N,) against (1, Q) broadcasts to the (N, Q) f_atoms layout.
    return vols .* exp.(.-(qvals' .^ 2) .* (vols .^ (2 / 3)) ./ (4π))
end

"""
    vacuo(mol, qvals, lMax, ions, energy, _CHUNK; form_factor_source) ->
    AbstractArray{<:Complex,3}

Vacuum term: the real atoms of `mol` with no solvent at all.

`B_vac`, the `(C, K, Q)` multipoles feeding `S_vac,·` downstream. The amplitude
is the true X-ray form factor per atom, pulled through the `Interfaces` facade
at photon energy `energy`; near an absorption edge it is complex
(`f0 + f' + i*f''`), so `compute_B_lm` returns two channels here where the dummy
species return one.

# Arguments
-   `mol`: the molecule; only its spherical coordinates are used.
-   `qvals::AbstractVector{<:Real}`, length `Q`: momentum-transfer grid in Å⁻¹.
-   `lMax::Integer`: maximum spherical harmonic degree.
-   `ions::Vector{String}`: ion/element string per atom, i.e. `elms(mol)`.
-   `energy::Float64`: photon energy in eV.
-   `_CHUNK::UInt64`: atoms processed per pass inside [`compute_B_lm`](@ref).

# Keywords
-   `form_factor_source::FormFactorSource = FORM_FACTOR_SOURCE`: the
    form-factor backend, mirroring `Molecule.create`'s `radii_source`. The
    default reads the bundled tables; a stub subtype lets this be exercised
    without them.

# Returns
-   `B_lm::AbstractArray{<:Complex,3}`, `(C, K, Q)` in [`compute_B_lm`](@ref).
"""
function vacuo(
    mol::Molecule,
    qvals::AbstractVector{<:Real},
    lMax::Integer,
    ions::Vector{String},
    energy::Float64,
    _CHUNK::UInt64;
    form_factor_source::FormFactorSource = FORM_FACTOR_SOURCE,
)::AbstractArray{<:Complex,3}
    crd = Molecules.coords_spherical(mol)
    tbl = Interfaces.form_factor_table(form_factor_source, energy, ions, qvals)
    amp = Interfaces.form_factors(tbl, ions, qvals)
    return compute_B_lm(crd, qvals, amp, lMax, _CHUNK)
end

"""
    excluded(mol, qvals, lMax, _CHUNK) -> AbstractArray{<:Complex,3}

Excluded-volume term: one Gaussian dummy per atom, at the atom's own position.
`B_ex`, feeding `S_ex,·` downstream. Bulk solvent cannot occupy the space an
atom already fills, so the volume each atom displaces has to be subtracted from
the amplitude before squaring. Every atom displaces solvent whether or not it is
on the surface, so this runs over the full molecule, unlike [`hydration`](@ref).

# Arguments
- `mol`: the molecule; its spherical coordinates and per-atom volumes are used.
- `qvals::AbstractVector{<:Real}`, length `Q`: momentum-transfer grid in Å⁻¹.
- `lMax::Integer`: maximum spherical harmonic degree.
- `_CHUNK::UInt64`: atoms processed per pass inside [`compute_B_lm`](@ref).

# Returns
-   `B_lm`: as [`vacuo`](@ref), but with `C = 1` -- a dummy sphere's amplitude is
    real, so there is no `f''` channel.
"""
function excluded(
    mol::Molecule,
    qvals::AbstractVector{<:Real},
    lMax::Integer,
    _CHUNK::UInt64,
)::AbstractArray{<:Complex,3}
    crd = Molecules.coords_spherical(mol)
    amp = _gaussian_dummy(Molecules.vols(mol), qvals)
    return compute_B_lm(crd, qvals, amp, lMax, _CHUNK)
end

"""
    hydration(mol, qvals, lMax, _CHUNK; thickness, probe, n_target, classes)
        -> @NamedTuple{convex::Array{ComplexF64,3}, concave::Array{ComplexF64,3}, cavity::Array{ComplexF64,3}}

Hydration-shell term, split into CRYSOL 3's three border-layer populations.

[`SASA.shell_points`](@ref) is run once; its beads are partitioned by
[`SASA.BeadClass`](@ref) and each class gets its own `B_lm` via
[`compute_B_lm`](@ref), every bead carrying the `area * thickness` slab of shell
it stands for, so the cloud tiles the layer rather than approximating it with an
envelope. The three arrays are the `sh_convex` / `sh_concave` / `sh_cavity`
species of the five-species expansion
    `A_total = A_vac - dns·A_ex + Σ_k dro_k·A_sh_k`; each takes its own fitted
contrast `dro_k` downstream (CRYSOL's `ρ = (1, 1, 0)` defaults).

# Arguments
- `mol`: the molecule; its accessible surface is used, not its atom positions.
- `qvals::AbstractVector{<:Real}`, length `Q`: momentum-transfer grid in Å⁻¹.
- `lMax::Integer`: maximum spherical harmonic degree.
- `_CHUNK::UInt64`: dummies processed per pass inside [`compute_B_lm`](@ref).

# Keywords
    - `thickness::Float64 = SHELL_THICKNESS`: shell thickness in Å; `> 0`.
    - `probe::Float64 = PROBE_RADIUS`: solvent probe radius, forwarded to `shell_points`.
    - `n_target::Union{Nothing,Int} = SHELL_N_TARGET`: total shell dummies before the class
    split, the analogue of CRYSOL's `--fb`. `nothing` lets `SASA.shell_points`
    size the cloud from the accessible area (`≈ area / SASA.SHELL_AREA_PER_POINT`,
    floored at `SASA.SHELL_MIN_POINTS`); pass an `Int` to pin it. Runtime is
    linear in it.
    - `classes = SHELL_CLASSES`: which populations to actually build. A class left
    out still appears in the result as an all-zero `B_lm` (equivalent to
    `dro_k = 0`); omitting it only skips its `compute_B_lm` pass. Must be
    non-empty.

# Returns
-   `@NamedTuple{convex, concave, cavity}` of `(C, K, Q)` `Array{ComplexF64,3}`,
    `C = 1` (dummy amplitudes are real).
"""
function hydration(
    mol::Molecule,
    qvals::AbstractVector{<:Real},
    lMax::Integer,
    _CHUNK::UInt64;
    thickness::Float64           = SHELL_THICKNESS,
    probe::Float64               = PROBE_RADIUS,
    n_target::Union{Nothing,Int} = SHELL_N_TARGET,
    classes                      = SHELL_CLASSES,
)::@NamedTuple{convex::Array{ComplexF64,3}, concave::Array{ComplexF64,3}, cavity::Array{ComplexF64,3}}
    thickness > 0.0 || throw(ArgumentError("hydration: thickness must be > 0"))
    isempty(classes) && throw(ArgumentError("hydration: classes must be non-empty"))

    pts, area, class = SASA.shell_points(mol; probe = probe, n_target = n_target)

    _shell(want) = begin
        sel = want in classes ? findall(==(want), class) : Int[]
        crd = Molecules.to_spherical(pts[:, sel])
        amp = _gaussian_dummy(area[sel] .* thickness, qvals)
        compute_B_lm(crd, qvals, amp, lMax, _CHUNK)
    end

    return (convex  = _shell(SASA.CONVEX),
            concave = _shell(SASA.CONCAVE),
            cavity  = _shell(SASA.CAVITY))
end
