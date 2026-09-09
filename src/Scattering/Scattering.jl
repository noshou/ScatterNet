"""
The SAXS/SANS forward model: a molecule and a `q` grid in, the orientationally
averaged detector intensity `I_calc(q)` out. 
# Module layout

- `SphFuncs`    -   `Y_lm`, `j_l`, normalised Legendre.
- `PartialWave` -   `compute_B_lm` (the multipole moments), plus
                    `self_scatter` / `cross_scatter` / `partial_wave_weights`
                    (the reductions to `S_ab(q)`).
- `Scatterers`  -   one builder per species: `vacuo`, `excluded`, `hydration`
                    (the last returns one `B_lm` per `SHELL_CLASSES` entry).
- `Intensity`   -  `gram` (the `S_ab` matrix `G`), `intensity` (`vᵀ G v`),
                    `intensity_calc` (`m·I + c`), `contrast_vector` /
                    `contrast_matrix`, and `excluded_volume_factor` (CRYSOL's
                    fitted excluded-volume radius `r₀`).
- `Forward`     -   the assembled model: `species_multipoles`, `gram_matrix`,
                    `forward_cache` and `forward`. Build the geometry-only
                    `ForwardCache` once with `forward_cache`, then call
                    `forward(cache, m, c, dns, ρ; r0)` per parameter set; or
                    `forward(mol, qvals, lMax, energy; m, c, dns, ρ, r0)` for a
                    one-off.

# Background

1.  For a fixed orientation, the coherent scattering amplitude of
    N point-like scatterers with form factors f_i(q) at positions r_i
    is a Fourier sum:

        `A(q) = Σ_i f_i(q) * exp(i q·r_i)`

    where q is the momentum-transfer vector. In solution scattering
    (SAXS/SANS), molecules tumble freely, so the measured intensity
    is the square of this amplitude averaged over every possible orientation:

        `I(q) = < |A(q)|² >_orientations`

    which is computationally intractible to calculate exactly for large molecules.
\\
2.  We can instead expand the plane wave in spherical harmonics, since
    `exp(i q·r)` has an exact expansion (the Rayleigh expansion) in terms
    of spherical Bessel functions j_l and spherical harmonics Y_lm:

        `exp(i q·r) = 4π * Σ_l Σ_m i^l * j_l(q r) * Y_lm(q_hat) * conj(Y_lm(r_hat))`

    where q_hat, r_hat are unit vectors and r = |r|.

    Substitute this into A(q) and swap the order of the atom-sum and the (l,m)-sum:

        `A(q) = 4π * Σ_l Σ_m i^l * Y_lm(q_hat) * B_lm(q)`

    where:

        `B_lm(q) = Σ_i f_i(q) * j_l(q r_i) * conj(Y_lm(θ_i, φ_i))`

    B_lm(q) is therefore the multipole moment of degree (l,m) of the
    scattering amplitude A(q). Squaring A(q) gives:

        `|A(q)|² =  (4π)² * Σ_{l,m} Σ_{l',m'} i^l*(-i)^l'
                    * Y_lm(q_hat) * conj(Y_l'm'(q_hat))
                    * B_lm(q) * conj(B_l'm'(q))`

    Averaging over orientation means averaging over the direction q_hat.
    Y_lm is orthonormal on the sphere, so `∫dΩ_q_hat Y_lm(q_hat)*conj(Y_l'm'(q_hat))`
    is 1 when (l,m)=(l',m') and 0 otherwise, collapsing the double sum to a single sum.
    The surviving l=l' diagonal's  phase factor becomes `i^l*(-i)^l = (i*(-i))^l = 1^l = 1`,
    cancelling. What survives both is:

        `I(q) = 4π * Σ_l Σ_{m=-l}^{l} |B_lm(q)|²`

    an orientational average computed once per atom, in closed form,
    instead of by averaging over rotations.
\\
3.  `B_lm` only needs to be stored for `m = 0, ..., l`, not the full
    `-l, ..., l` range, since spherical harmonics satisfy
    `Y_{l,-m} = (-1)^m * conj(Y_lm)`, so for a real `f_i(q)` the same
    identity forces `B_{l,-m} = (-1)^m * conj(B_lm)`, hence `|B_{l,-m}|² = |B_lm|²`:

        `Σ_{m=-l}^{l} |B_lm|² = 1·|B_l0|² + Σ_{m=1}^{l} 2·|B_lm|²`

    a weighting given `1`-for-`m=0` and  `2`-for-`m>0`  over the half that computed.

    The derivation needs `f_i(q)` real. Near an absorption edge, atomic
    form factors are complex (`f = f0 + f' + i*f''`, the anomalous term),
    and the identity in (3) breaks exactly where it needed `conj(f_i) = f_i`.
    The fix: split `f_i = Re(f_i) + i*Im(f_i)` and build `B_lm` separately
    for each real-valued piece (two "channels"). Each channel individually
    satisfies (3) exactly, and the cross term between channels is odd in
    `m` and cancels once summed over the full `-l..l` range, so summing
    the channels' `|B_lm|²` incoherently reproduces the exact total with
    no approximation.
\\
4.  A real molecule scatters as several **species** superposed at the
    amplitude level (coherently), then squared:

        `A_total(q) = Σ_a c_a A_a(q)`

    so, after the orientational average of (2),

        `I(q) = Σ_a Σ_b c_a c_b S_ab(q)`,
        `S_ab(q) = 4π Σ_lm w_lm Re(B^a_lm(q) conj(B^b_lm(q)))`

    `S_aa` (a species against itself) is the self term `4π Σ_lm w_lm |B_lm|²`;
    `S_ab` with `a ≠ b` is the cross term. The species:

        - `vac`         molecule as if no solvent existed.
        - `ex`          one Gaussian excluded-volume dummy per atom: the bulk
                        solvent each atom displaces (a negative contrast).
        - `sh_convex`  ┐ hydration-shell dummies on the solvent-accessible
        - `sh_concave` ├ surface, split by local geometry into CRYSOL 3's
        - `sh_cavity`  ┘ three border-layer populations (`SHELL_CLASSES`)

    Classic single-shell CRYSOL is the `n = 3` reduction `(vac, ex, sh)` with
    the three shell classes merged.
\\
5.  The dummy species do not carry the true local electron density

        `v = (1, -dns, dro_1, dro_2, dro_3)`,     dro_k = DRO_UNIT * ρ_k

    - `dns` rescales `A_ex` to the mean electron density of the displaced
            bulk solvent (`≈ 0.334 e·Å⁻³`).
    - `ρ_k` dimensionless shell contrast per class (CRYSOL's `--dro`
            multiple, default `(1, 1, 0)`); `dro_k` is the class's excess
            electron density over bulk.

    With the geometry-only Gram matrix `G_ab(q) = S_ab(q)` the whole
    expansion is the bilinear form

        `I(q) = vᵀ G(q) v`,     G(q) ⪰ 0

    `n(n+1)/2` distinct `S_ab(q)` curves (15 for `n = 5`, 6 for `n = 3`).
    `G` depends only on geometry and beam; `v` only on the fit parameters, so
    `G` is built once per structure and reused across every parameter set.
\\
6.  A table of atomic radii gets the *displaced* volume wrong -- how much bulk
    water an atom excludes depends on its chemical environment, not just its
    element. CRYSOL's fix is one global expansion factor `c_1 = r_0/r_m` over
    all dummies, `r_0` fitted and `r_m` the structure's mean atomic radius.
    Expanding a dummy's radius sends `V_j -> c_1^3 V_j` in the Gaussian of (4),
    i.e.

        `f_j(q) -> c_1^3 * f_j(q) * exp(-q^2 (c_1^2 - 1) V_j^(2/3) / 4π)`

    The residual envelope still carries `V_j`, so exactly it does not leave the
    atom sum and `B_ex` (hence `G`) would have to be rebuilt per `r_0`. CRYSOL's
    standard approximation -- kept here -- replaces the per-atom `V_j^(2/3)` by
    the mean-radius value `V_m^(2/3) = (4π/3)^(2/3) r_m^2`, making it one scalar
    function of `q` that leaves the sum entirely:

        `G_ex(q) = c_1^3 exp(-q^2 (c_1^2 - 1) (4π/3)^(2/3) r_m^2 / 4π)`

    So `r_0` reweights the *contrast*, not the geometry: `v` simply becomes
    `q`-dependent in its `ex` entry, and `G` stays cached.

        `v(q) = (1, -dns*G_ex(q), dro_1, dro_2, dro_3)`,   `I(q) = v(q)ᵀ G(q) v(q)`

    Exact at `q = 0` (total excluded volume scales by `c_1^3`) and for an atom
    of radius `r_m`; it degrades with the spread of radii about `r_m`, which for
    protein heavy atoms is small. `r_0 = r_m` gives `G_ex ≡ 1`, the uncorrected
    model. Note `dns` and `r_0` are strongly degenerate -- CRYSOL fixes `dns` at
    `0.334` and fits `r_0`.
\\
7.  A real detector reads an arbitrary scale over an imperfect buffer
    subtraction, so the reported intensity is

        `I_calc(q) = m * I(q) + c`

    with `m` the overall scale (absolute → arbitrary units) and `c` a flat
    background.
"""
module Scattering

using ..Interfaces: Interfaces, FormFactorSource, FormFactorSourceTables
using ..Molecule.SASA: SASA

# The public surface. `forward` is the forward model; `gram_matrix` is the
# geometry-only pass to cache when sweeping fit parameters. Everything the
# `include`s below bring in (`compute_B_lm`, `vacuo`/`excluded`/`hydration`,
# `gram`/`intensity`/…) is the machinery those two compose -- reachable by
# qualified name for tests and advanced callers, but not part of the API.
export forward, gram_matrix, forward_cache

# ---------------------------------------------------------------------------
# Module-level configuration -- the single source of truth for every knob the
# forward model takes. `vacuo` / `hydration` / `Intensity.jl` / `Forward.jl`
# read these as their defaults; override per call where needed.
# ---------------------------------------------------------------------------

"""
Hydration-shell thickness in Å: how far the perturbed-density water layer
extends beyond the solvent-accessible surface. `3.0` is CRYSOL's border-layer
default.
"""
const SHELL_THICKNESS = 3.0

"Solvent probe radius in Å (water), forwarded to `SASA.shell_points`."
const PROBE_RADIUS = 1.4

"""
Default shell-dummy budget (`hydration`'s `n_target`). `nothing` lets
`SASA.shell_points` size the cloud from the accessible area
(`≈ area / SASA.SHELL_AREA_PER_POINT`, floored at `SASA.SHELL_MIN_POINTS`);
an `Int` pins it.
"""
const SHELL_N_TARGET::Union{Nothing,Int} = nothing

"""
CRYSOL 3's three border-layer populations, in the field order of the
`NamedTuple` [`hydration`](@ref) returns. Each carries its own fitted contrast
`dro_k = DRO_UNIT * ρ_k` downstream; CRYSOL's defaults are `ρ = (1, 1, 0)`.
"""
const SHELL_CLASSES = (SASA.CONVEX, SASA.CONCAVE, SASA.CAVITY)

"Shell-contrast unit in e·Å⁻³ (CRYSOL's `--dro`); `dro_k = DRO_UNIT * ρ_k`."
const DRO_UNIT = 0.03

"Default X-ray form-factor backend (bundled Waasmaier-Kirfel + Chantler tables)."
const FORM_FACTOR_SOURCE = FormFactorSourceTables()

"Atoms/dummies per pass in `compute_B_lm`."
const B_LM_CHUNK = UInt64(2048)

include("SphFuncs.jl")
include("PartialWave.jl")
include("Scatterers.jl")
include("Intensity.jl")
include("Forward.jl")

using .SphFuncs: SphFuncs

end # module
