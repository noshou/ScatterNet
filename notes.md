# SASA calculation

Purpose: per-atom solvent-accessible surface area + exposed surface points,
used to place the hydration / border layer term B(q) in the forward model.
Runs offline (once per structure), not inside the sampler, so no gradients
needed. Shell contrast is a fitted nuisance param -> few-% accuracy is fine.

## Method: Shrake-Rupley (numerical point sampling)

- Primary ref: Shrake & Rupley (1973), J. Mol. Biol. 79(2):351-371,
  doi:10.1016/0022-2836(73)90011-9
- Defining surface: Lee & Richards (1971), J. Mol. Biol. 55(3):379-400
- Reference impl to validate against: Biopython `Bio.PDB.SASA` (ShrakeRupley),
  and FreeSASA (Mitternacht 2016, F1000Research 5:189) as a test oracle
- Practical writeup: https://pmc.ncbi.nlm.nih.gov/articles/PMC2712621/

Per atom i:
  Ri = r_vdW_i + probe        (probe = 1.4 A, effective water radius; convention)
  place N points at center_i + Ri * unit_sphere_point
  point is buried if within Rj of any neighbor j
  f_i = (# non-buried) / N
  SASA_i = f_i * 4*pi*Ri^2
Neighbor j of i if dist(i,j) < Ri + Rj. Cell-list/grid for O(N).

## Point distribution

- Fibonacci / golden-spiral sphere: near-optimal uniformity for fixed N,
  deterministic, what Biopython uses. Closed set: N baked in, not extensible.
- Deterministic-random (seeded, e.g. normalized 3D Gaussian): valid & uniform,
  but high discrepancy (clumps/gaps) -> ~3-5x more points for same accuracy,
  and a single frozen set gives structured per-atom bias (bad, since we use
  relative per-atom exposure to weight the shell). Rejected.
- Plastic / R2 sequence (Martin Roberts,
  https://extremelearning.net/en/blog/the-unreasonable-effectiveness-of-quasirandom-sequences/):
  t_k = frac(k * (1/rho, 1/rho^2)), rho = plastic number (x^3 = x + 1),
  mapped to sphere via EQUAL-AREA map (x->longitude, y->z=2y-1).
  Slightly higher discrepancy than fixed-N Fibonacci, but OPEN/extensible:
  every prefix is well-distributed, so points can be appended and prior
  work reused.

Default: Fibonacci for a plain fixed-N impl (better distribution, matches
reference impls). Use R2 only if doing adaptive refinement (below), which
needs incremental extension.

## Adaptive refinement

Idea: atoms need unequal accuracy. Buried and 
atoms converge at tiny N; only occlusion-boundary atoms (mid-range f, error
~ sqrt(f(1-f)/N)) need many points. Start small per atom, extend only the
ambiguous ones, reusing points via the R2 sequence.

Guards:

- Classify with a confidence bound, not the point estimate: with k exposed of
  n, k=0 still allows true f up to ~3/n (rule of three). Reject as buried only
  if f_upper * 4*pi*Ri^2 < area_tol. n_min ~ 50-100 keeps this safe.
- Free pre-filter from neighbor list: if no neighbor reaches Ri's surface
  (dist(i,j) >= Ri + Rj for all j) then f_i = 1 exactly, no sampling.
- Stop boundary atoms when SE(f) = sqrt(f(1-f)/n), times 4*pi*Ri^2, < tol,
  or n hits n_cap (~960). Double n each round (64 -> 128 -> 256 ...).
- Scope: "reject buried" applies ONLY to SASA / shell. Excluded-volume term
  C(q) still uses every atom's full volume. Buried atom -> zero surface points;
  shell construction must handle an empty set.

Pseudocode:

    build neighbor list (cell list / grid)              # shared, once

    for each atom i:
        neighbors_i = cell_list.query(i, Ri + Rmax)
        if no neighbor reaches Ri surface:
            f_i = 1; continue                            # pre-filter

        n, k = 0, 0
        loop:
            add next block of points (n -> 2n) from R2 sequence
            k += count(points in block not inside any neighbor_i)
            f = k / n
            if f_upper(k, n) * 4*pi*Ri^2 < area_tol:
                f_i = 0; break                           # buried
            if SE(f, n) * 4*pi*Ri^2 < tol or n >= n_cap:
                f_i = f; break                           # converged

        SASA_i = f_i * 4*pi*Ri^2
        # keep surviving (non-buried) points for hydration-shell placement

## Codebase hooks

- `Molecule.coords` : (3,n) centered xyz  (added for this)
- `Molecule.radii`  : (n,) per-atom vdW/ionic radius  (added for this)
- validate total + per-atom SASA against FreeSASA on 1CRN (crambin)
