"""
Scattering primitives implemented so far: spherical harmonics / Bessel
functions (`SphFuncs`) and X-ray form factors (`FormFactorXrayDB`); the
`B_lm`/`S(q)` combination described below is the model this module is being
built toward, not yet code that lives here. Depends on the top-level
`Interfaces` markers.

# Background

1.  For a fixed orientation, the coherent scattering amplitude of 
    N point-like scatterers with form factors f_i(q) at positions r_i 
    is  a Fourier sum: 

        `A(q) = Σ_i f_i(q) * exp(i q·r_i)`

    where q  is the momentum-transfer vector. In solution scattering 
    (SAXS/SANS), molecules tumble freely, so the measured intensity 
    is the square of this amplitude averaged over every possible orientation: 

        `I(q) = < |A(q)|² >_orientations`

    which is computationally intractible to calculate exactly. 
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
    scattering amplitude A(q). Squaring A(q) gives a DOUBLE sum, one
    (l,m) from A and one (l',m') from conj(A):

        `|A(q)|² =  (4π)² * Σ_{l,m} Σ_{l',m'} i^l*(-i)^l' 
                    * Y_lm(q_hat) * conj(Y_l'm'(q_hat)) 
                    * B_lm(q) * conj(B_l'm'(q))`

    Averaging over orientation means averaging over the direction q_hat
    (B_lm depends only on the molecule's own coordinates, so it rides
    along unchanged). Two things happen to that bracketed double sum at
    once: Y_lm is orthonormal on the sphere, so
    `∫dΩ_q_hat Y_lm(q_hat)*conj(Y_l'm'(q_hat))` is 1 when (l,m)=(l',m')
    and 0 otherwise, collapsing the double sum to a single one; and on
    that surviving l=l' diagonal, the phase factor becomes
    `i^l*(-i)^l = (i*(-i))^l = 1^l = 1`, cancelling exactly rather than
    merely averaging away. What survives both is:

        `I(q) = 4π * Σ_l Σ_{m=-l}^{l} |B_lm(q)|²`

    -- an orientational average computed once per atom, in closed form,
    instead of by averaging over rotations.
\\
3.  `B_lm` only needs to be stored for `m = 0, ..., l`, not the full
    `-l, ..., l` range: spherical harmonics satisfy
    `Y_{l,-m} = (-1)^m * conj(Y_lm)`, so for a real `f_i(q)` the same
    identity forces `B_{l,-m} = (-1)^m * conj(B_lm)`, hence
    `|B_{l,-m}|² = |B_lm|²`:

        `Σ_{m=-l}^{l} |B_lm|² = 1·|B_l0|² + Σ_{m=1}^{l} 2·|B_lm|²`

    a `1`-for-`m=0`, `2`-for-`m>0` weighting over the half that is
    actually computed.
\\
4.  That derivation needs `f_i(q)` real. Near an absorption edge, atomic
    form factors are complex (`f = f0 + f' + i*f''`, the anomalous term),
    and the identity in (3) breaks exactly where it needed `conj(f_i) = f_i`.
    The fix: split `f_i = Re(f_i) + i*Im(f_i)` and build `B_lm` separately
    for each real-valued piece (two "channels"). Each channel individually
    satisfies (3) exactly, and the cross term between channels is odd in
    `m` and cancels once summed over the full `-l..l` range, so summing
    the channels' `|B_lm|²` incoherently reproduces the exact total with
    no approximation.
\\
5.  A real molecule is modelled as several distinct species of scatterer
    at once, which combine into one total amplitude before squaring, not
    just added as separate intensities:

        `A_total(q) = A_vac(q) - A_ex(q) + A_sh(q)`

    expanding into:

        `I(q)   = <|A_total(q)|²>
                =   S_vac,vac + S_ex,ex + S_sh,sh
                    - 2*S_vac,ex + 2*S_vac,sh - 2*S_ex,sh`

    - `vac`: the real atoms in vacuo (ie: no solvent).
    - `ex`: Gaussian excluded-volume dummies, one per atom. A negative
    contrast correction subtracting the bulk solvent that each atom's own
    volume displaces (solvent can't occupy space an atom already fills).
    - `sh`: hydration-shell dummies, one per solvent-exposed atom (from the
    SASA patch geometry); a positive contrast contribution from the
    perturbed-density solvent layer coating the molecule's actual surface.

    A diagonal term (`S_vac,vac`, ...) is `4π * Σ_lm w_lm * |B_lm|²`, one
    species against itself. A cross term (`S_vac,ex`, ...) is
    `4π * Σ_lm w_lm * Re(B_a * conj(B_b))` between two different species'
    `B_lm`; the `+1`/`-1` in front of each term above is just each
    species' own sign in `A_total`, squared for the diagonal and
    multiplied pairwise for the cross terms.
\\
6.  Real molecules don't fit that combination perfectly: the baseline
    excluded-volume and hydration-shell terms are only generalised
    estimates (a nominal dummy-sphere density, an assumed shell
    thickness), not the true local electron density. Two fit parameters,
    `dns` and `dro`, correct for that by rescaling each species'
    amplitude rather than being separate terms of their own:

        `A_total(q) = A_vac(q) - dns*A_ex(q) + dro*A_sh(q)`

    - `dns` rescales the excluded-volume term to the true mean electron
    density of the displaced bulk solvent, since the dummy-sphere
    model's own assumed density is only nominal.
    - `dro` is the hydration shell's excess electron density over bulk
    solvent, since ordered/perturbed water at the surface is measurably
    denser than bulk water by an amount not known a priori.

    This expands (same squaring as (5), now with `dns`/`dro` folded into
    each species' sign) into:

        `I(q)   = <|A_total(q)|²>
                =   S_vac,vac
                    - 2*dns*S_vac,ex
                    + 2*dro*S_vac,sh
                    + dns^2*S_ex,ex
                    + dro^2*S_sh,sh
                    - 2*dns*dro*S_ex,sh`
\\
7.  `dns`/`dro` correct the model's own approximations, but a real
    detector measures on an arbitrary intensity scale, and real buffer
    subtraction is rarely perfect. Two more parameters absorb that:

        `I_calc(q) = m*I(q) + c`

    - `m`: overall scale between calculated (absolute) and measured (arbitrary-unit) intensity.
    - `c`: constant background left over from imperfect buffer subtraction.

    The full parameter set actually fit to data is therefore `[m, c, dns, dro]`, on top of 
    whatever geometry (SH coefficients) is being fit at the same time.
"""
module Scattering

include("SphFuncs.jl")
include("FormFactorXrayDB.jl")

using .SphFuncs: SphFuncs
using .FormFactorXrayDB: FormFactorXrayDB

end # module
