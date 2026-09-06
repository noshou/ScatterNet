# Per-species partial-wave terms: one function per scatterer species in
# `A_total(q) = A_vac(q) - dns*A_ex(q) + dro*A_sh(q)` (see the `Scattering`
# docstring, section 5/6). Each builds that species' `(N, Q)` amplitude, hands
# it to the shared `compute_B_lm`, and returns just that species' multipoles.
using ..Interfaces: Interfaces, FormFactorSource, FormFactorSourceXrayDB
using ..Molecule.SASA: SASA
using ..Molecule.Molecules: Molecules, Molecule

"""
Hydration-shell thickness in Å: how far the perturbed-density water layer
extends beyond the solvent-accessible surface. (CRYSOL default).
"""
const SHELL_THICKNESS = 3.0

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
at photon energy `energy`; near an absorption edge it is complex (`f0 + f' + i*f''`),
so `compute_B_lm` returns two channels here where the dummy species return one.

# Arguments
-   `mol`: the molecule; only its spherical coordinates are used.
-   `qvals::AbstractVector{<:Real}`, length `Q`: momentum-transfer grid in Å⁻¹.
-   `lMax::Integer`: maximum spherical harmonic degree.
-   `ions::Vector{String}`: ion/element string per atom, i.e. `elms(mol)`.
-   `energy::Float64`: photon energy in eV.
-   `_CHUNK::UInt64`: atoms processed per pass inside [`compute_B_lm`](@ref).

# Keywords
-   `form_factor_source::FormFactorSource = FormFactorSourceXrayDB()`: the
    form-factor backend, mirroring `Molecule.create`'s `radii_source`. The
    default needs the live xraydb extension; a stub subtype lets this be
    exercised without one.

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
    form_factor_source::FormFactorSource = FormFactorSourceXrayDB()
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
-   `B_lm`: as [`vacuo`](@ref), but with `C = 1` a dummy sphere's
    amplitude is real, so there is no `f''` channel.
"""
function excluded(
    mol::Molecule,
    qvals::AbstractVector{<:Real},
    lMax::Integer,
    _CHUNK::UInt64
)::AbstractArray{<:Complex,3}
    crd = Molecules.coords_spherical(mol)
    amp = _gaussian_dummy(Molecules.vols(mol), qvals)
    return compute_B_lm(crd, qvals, amp, lMax, _CHUNK)
end

"""
    hydration(mol, qvals, lMax, _CHUNK; thickness, probe, n_target, classes) ->
    AbstractArray{<:Complex,3}

Hydration-shell term: a Gaussian dummy on every accessible surface point.

Dummies come from [`SASA.shell_points`](@ref), each carrying the
`area * thickness` slab of shell it stands for, so the cloud tiles the layer
rather than generalizing an envlope around it.


# Arguments
- `mol`: the molecule; its accessible surface is used, not its atom positions.
- `qvals::AbstractVector{<:Real}`, length `Q`: momentum-transfer grid in Å⁻¹.
- `lMax::Integer`: maximum spherical harmonic degree.
- `_CHUNK::UInt64`: dummies processed per pass inside [`compute_B_lm`](@ref).

# Keywords
- `thickness::Float64 = SHELL_THICKNESS`: shell thickness in Å; `> 0`.
- `probe::Float64 = 1.4`: solvent probe radius, forwarded to `shell_points`.
- `n_target::Union{Nothing,Int} = nothing`: total shell dummies, the analogue of
CRYSOL's `--fb`. `nothing` lets `SASA.shell_points` size the cloud from the
accessible area (`≈ area / SASA.SHELL_AREA_PER_POINT`, floored at
`SASA.SHELL_MIN_POINTS`); pass an `Int` to pin it. Runtime is linear in it.

# Returns
-   `B_lm`: as [`excluded`](@ref) (`C = 1`). A molecule with no accessible
    surface yields an all-zero `B_lm`, rather than an error.
"""
function hydration(
    mol::Molecule,
    qvals::AbstractVector{<:Real},
    lMax::Integer,
    _CHUNK::UInt64;
    thickness::Float64           = SHELL_THICKNESS,
    probe::Float64               = 1.4,
    n_target::Union{Nothing,Int} = nothing,
    classes                      = (SASA.CONVEX, SASA.CONCAVE)
)::AbstractArray{<:Complex,3}
    thickness > 0.0 || throw(ArgumentError("hydration: thickness must be > 0"))
    isempty(classes) && throw(ArgumentError("hydration: classes must be non-empty"))

    pts, area, class = SASA.shell_points(mol; probe = probe, n_target = n_target)
    keep = findall(c -> c in classes, class)
    crd = Molecules.to_spherical(pts[:, keep])
    amp = _gaussian_dummy(area[keep] .* thickness, qvals)
    return compute_B_lm(crd, qvals, amp, lMax, _CHUNK)
end
