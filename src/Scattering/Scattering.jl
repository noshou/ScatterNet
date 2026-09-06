"""
Scattering primitives.

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
4.  A real molecule is modelled as several distinct scatterers
    at once, which combine into one total amplitude before squaring:

        `A_total(q) = A_vac(q) - A_ex(q) + A_sh(q)`

    expanding into:

        `I(q)   = <|A_total(q)|²>
                =   S_vac,vac + S_ex,ex + S_sh,sh
                    - 2*S_vac,ex + 2*S_vac,sh - 2*S_ex,sh`
    where: 
        - `vac`: the real atoms in vacuo (ie: no solvent).
        - `ex`: Gaussian excluded-volume dummies, one per atom. A negative
                contrast correction subtracting the bulk solvent that each 
                atom's own volume displaces.
        - `sh`: hydration-shell dummies, one per solvent-exposed atom. A 
                positive contrast contribution from the perturbed-density 
                solvent layer coating the molecule's actual surface.

    A diagonal term (`S_vac,vac`, ...) is `4π * Σ_lm w_lm * |B_lm|²`, which calculates
    the self-interaction of a scatterer. A cross term (`S_vac,ex`, ...) is
    `4π * Σ_lm w_lm * Re(B_a * conj(B_b))` between two different scatterer's `B_lm`.
\\
5.  Real molecules don't perfectly fit the model described above. The baseline
    excluded-volume and hydration-shell terms are only generalised
    estimates (a nominal dummy-sphere density, and an assumed shell
    thickness), not the true local electron density. Two fit parameters,
    `dns` and `dro`, correct for that by rescaling each scatterer's amplitudes:

        `A_total(q) = A_vac(q) - dns*A_ex(q) + dro*A_sh(q)`
    
    where:
            - `dns`:    scaling factor for the excluded-volume term. Recales to the 
                        true mean electron density of the displaced bulk solvent.
            - `dro`:    the hydration shell's excess electron density over bulk
                        solvent, since ordered/perturbed water at the surface is 
                        denser than bulk water by an amount not known a priori.

    This expands into:

        `I(q)   = <|A_total(q)|²>
                =   S_vac,vac
                    - 2*dns*S_vac,ex
                    + 2*dro*S_vac,sh
                    + dns^2*S_ex,ex
                    + dro^2*S_sh,sh
                    - 2*dns*dro*S_ex,sh`
\\
6.  `dns`/`dro` correct the model's own approximations, but a real detector measures on 
    an arbitrary intensity scale, and real buffer subtraction is imperfect. Two more 
    parameters account for this:

        `I_calc(q) = m*I(q) + c`
    
    where:
        - `m`:  a scaling factor to rescale the calculated (absolute) intensity to the
                measured (arbitrary-unit) intensity.
        - `c`:  a constant background left over from imperfect buffer subtraction.
"""
module Scattering

include("SphFuncs.jl")
include("PartialWave.jl")
include("Scatterers.jl")

using .SphFuncs: SphFuncs

end # module
