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
  SASA_i = f_i * 4*π*Ri^2
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
  if f_upper * 4*π*Ri^2 < area_tol. n_min ~ 50-100 keeps this safe.
- Free pre-filter from neighbor list: if no neighbor reaches Ri's surface
  (dist(i,j) >= Ri + Rj for all j) then f_i = 1 exactly, no sampling.
- Stop boundary atoms when SE(f) = sqrt(f(1-f)/n), times 4*π*Ri^2, < tol,
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
            if f_upper(k, n) * 4*π*Ri^2 < area_tol:
                f_i = 0; break                           # buried
            if SE(f, n) * 4*π*Ri^2 < tol or n >= n_cap:
                f_i = f; break                           # converged

        SASA_i = f_i * 4*π*Ri^2
        # keep surviving (non-buried) points for hydration-shell placement

## Codebase hooks

- `Molecule.coords` : (3,n) centered xyz  (added for this)
- `Molecule.radii`  : (n,) per-atom vdW/ionic radius  (added for this)
- validate total + per-atom SASA against FreeSASA on 1CRN (crambin)

---

---

# Form-factor backend, test coverage, and design audits

Later working notes, appended to the SASA notes above. The plastic/R2 rationale
in "Point distribution" is revisited and measured in section 4.

## 1. Replacing the Python `xraydb` bridge with bundled tables

### What was there

`f(q,E) = f0(s) + f1(E) + i·f2(E)` came from the Python `xraydb` package,
reached over PythonCall + CondaPkg through a package extension:

```
src/Interfaces/FormFactorXrayDB/FormFactorXrayDB.jl   Julia side, FF container, tier logic docs
src/Interfaces/FormFactorXrayDB/py/FormFact_py.py     the actual computation (4 xraydb calls)
ext/FormFactorXrayDBExt.jl                            the PythonCall crossing
CondaPkg.toml                                         provisioned python + numpy + xraydb
```

The Python surface was four calls: `f0`, `f1_chantler`, `f2_chantler`,
`chantler_energies`.

### Why replace it

Measured before touching anything (aarch64, best-of-N):


|                               |          |
| ------------------------------- | ---------- |
| `import PythonCall`           | 6.23 s   |
| full load with ScatterNet     | 11.13 s  |
| first form-factor call (cold) | 2.92 s   |
| warm, 10 ions × 101 q        | 17.33 ms |

The scaling is the interesting part: **1.76 ms per unique ion, and essentially
independent of the q grid** — 11 q points cost 17.32 ms, 1001 q points cost
17.57 ms. A 91× increase in real work for 1.4% more time. So essentially none of
it was arithmetic.

Breaking it down inside the conda Python, with no PythonCall involved:


| component                                  | time         | share |
| -------------------------------------------- | -------------- | ------- |
| `f1_chantler` + `f2_chantler`, 10 elements | 8.87 ms      | 59%   |
| `chantler_energies`, 10 elements           | 3.87 ms      | 26%   |
| `f0`, 10 ions × 101 q                     | 0.62 ms      | 4%    |
| rest of the Python function                | ~1.7 ms      | 11%   |
| **whole Python function**                  | **15.03 ms** |       |
| PythonCall marshalling (17.33 − 15.03)    | 2.30 ms      | 13%   |

### What was pulled, and from where

Extracted by `src/Interfaces/FormFactor/extract.py` from **xraydb 4.5.8's
`xraydb.sqlite`**, whose LICENSE places that file and `data_sources/` under
**CC0 1.0**. Full citations in `src/Interfaces/FormFactor/README.md`.

- `Waasmaier` table → `waasmaier`: 211 species. Waasmaier & Kirfel (1995),
  *Acta Cryst.* A51, 416. Verified byte-identical against DABAX
  `f0_WaasKirf.dat` (max coefficient difference 0.0).
- `Chantler` table → `chantler`: 92 elements, 133,950 grid points. Chantler
  (1995) *JPCRD* 24, 71 and (2000) 29, 597, as distributed by NIST FFAST. Per
  XrayDB's README these values were supplied directly by C.T. Chantler and are
  a finer grid than the public FFAST web tabulation.

**Licence trap avoided:** DABAX also ships `f1f2_Chantler.dat`, whose header
reads *"The present license has been purchased by the ESRF Programming Group. No
use of these data is allowed from outside ESRF."* The MIT licence on the
DabaxFiles repo covers the repo, not that embedded restriction. Taking the
Chantler data from xraydb (CC0) is clean; taking it from DABAX would not be.

Result: **3.36 MB** `form_factors.sqlite3`, replacing a 10.3 MB sqlite plus a
511 MB conda environment.

### Implementation details that are load-bearing

**Floating-point association in f0.** NumPy evaluates `-e*q*q` as `((-e)*q)*q`.
Writing the mathematically identical `exp(-(b*s^2))` in Julia drifts 1–3 ulp on
8 of 110 reference points. Writing `exp(((-b)*s)*s)` is bit-identical. That is
why `FormFactor.f0` looks the way it does; do not "clean it up".

**Interpolation must be reproduced exactly.** `f1` is a *local 7-point*
not-a-knot cubic spline over `max(1, j−3) … min(n, j+3)`; `f2` is linear in
log–log space with a `1e-99` clamp. At 8 keV the grid is coarse enough (1%
spacing) that the scheme barely matters, but near an edge it dominates: at
7112.5 eV five plausible alternatives span 0.11 e in f₁ and 7% in f₂, and plain
linear-in-E is 7.8e-4 off — failing the suite's 1e-6 tolerance outright.

One agent recommended switching f₁ to linear-in-E for cleaner AD derivatives.
**Rejected**, for two reasons: it breaks the pinned reference values, and the AD
argument does not hold — stage 1 fits `[m, c, dns, dro]` plus SH coefficients,
and form factors are *constants* with respect to every one of those. AD through
f₁/f₂ would only matter if photon energy became a fit parameter.

**Validation.** Bit-compatibility was the acceptance criterion, and it held
where it matters: the pinned physics values are reproduced exactly and
`test/test_formfactor.jl`'s reference values needed **no regeneration**.

```
fe3+ @8keV q=0.1: 21.73320334 + 3.20285267im   ← pinned value, exact
fe3+ @8keV q=0.2: 21.71503439 + 3.20285267im   ← pinned value, exact
```

An early 110-point spot-check read as fully bit-exact on `f0`/`f2`. **That was
too small a sample and the claim was wrong.** Widening it to the 1200-point
fixture set (`test/fixtures/`) showed:

```
f0  : 347/348 bit-exact  (u6+ at s = 1.989 is 1 ulp out), worst rel 2.0e-16
f1  : 852 points,        worst |Δ| = 5.7e-14
f2  : ~99% bit-exact,    worst rel 1.2e-16
```

Two distinct causes, worth separating:

- **A real bug of mine, found by the wider fixture and fixed.** NumPy's `interp`
  computes `slope = (y2−y1)/(x2−x1)` and then `slope*(x−x1) + y1`. The first
  implementation folded that as `(x−x1)*(y2−y1)/(x2−x1)` — same value
  mathematically, 1 ulp away on ~2% of grid points. Fixing the association took
  `f2` deviations from 14 to 11.
- **A residue that is not matchable**: NumPy's and Julia's `log`/`exp` differ by
  ≤1 ulp on some arguments. Chasing it would pin the implementation to a NumPy
  detail forever, for a difference ~12 orders below the 0.4% spread between
  independent form-factor tabulations. The tests assert a scale-aware 2-ulp
  bound and say why, rather than claiming an exactness that cannot hold.

The lesson worth keeping: the original spot-check was drawn from too narrow a
species/energy range to expose either. The fixtures deliberately reach anions,
the `cval`/`siva` valence states, Z up to 98, both ends of Chantler's energy
range, and 500 points at 1 eV spacing across the Fe K edge.

**Speed after:**


| workload          | Julia     | Python bridge | speedup |
| ------------------- | ----------- | --------------- | --------- |
| 10 ions × 11 q   | 10.7 µs  | 17.32 ms      | 1619×  |
| 10 ions × 101 q  | 38.9 µs  | 17.33 ms      | 446×   |
| 10 ions × 1001 q | 284.5 µs | 17.57 ms      | 62×    |
| 50 ions × 1001 q | 1.43 ms   | —            | —      |

Unlike the Python path it now scales with the q grid, because it is doing
arithmetic rather than paying fixed per-ion overhead.

### Deliberate deviations from upstream

1. **Cs grid deduplicated**, 1504 → 1502 points. Upstream has two exact duplicate
   energies (11.4, 13.1 eV) and `UnivariateSpline(s=0)` *raises* on them —
   `f1_chantler('Cs', 11.4)` throws in Python. Ours does not.
2. **Unused columns dropped** (`mu_photo`, `mu_incoh`, `mu_total`, `sigma_mu`,
   `mue_f2`, `corr_*`). Nothing read them; the `corr_*` values are already folded
   into the stored `f1`.
3. **`NEUTRAL <ion>` is now logged.** The Python caught an unknown-ion lookup
   failure and silently retried with the bare element, so `fe4+` scattered as 26
   electrons instead of 22 with no record. The fallback is kept — usually the
   right answer — but it is no longer silent.
4. **`s` range is enforced** (`0 ≤ s ≤ 6`). WK's ionic fits carry large negative
   constants (`fe3+`: c = −61.93) and go negative above s ≈ 9. Python
   extrapolated silently.
5. **`f1`'s "3/5 CL" relativistic convention is documented.** Inherited from
   XrayDB, not Chantler's H82 default; a constant per-element offset (U differs
   by 1.08 e between conventions). It is what the pinned values encode.

### Test guard

`test/test_formfactor.jl`'s anti-vacuity check originally asserted that the
xraydb package extension was loaded and owned the answering method — a guard on
*who* answered. With the backend in-package there is no stub to fall through to,
so that mechanism could not be written at all. The first replacement (checking
`S_MAX == 6.0` and that `f0` returned something plausible) was **weaker than what
it replaced**, which was the wrong trade.

It is now a guard on *what* is answered: every one of the ~1200 fixture points
asserted to within 2 ulp, plus a guard-on-the-guard that the fixture files are
not truncated or swapped (species count, presence of specific anions/valence
states/Z>92, both ends of the element range, and the near-edge sample density).
That is strictly stronger than the provenance check, and it is what caught the
`interp` association bug above.

### Files

```
added:   src/Interfaces/FormFactor/FormFactor.jl        the module
         src/Interfaces/FormFactor/form_factors.sqlite3 3.36 MB, two tables + provenance
         src/Interfaces/FormFactor/extract.py           reproducible derivation
         src/Interfaces/FormFactor/README.md            provenance, conventions, licensing
         test/fixtures/fx_f0.csv, fx_f1f2.csv           1200 oracle points from the old backend
         test/fixtures/README.md                        provenance + why each tolerance differs
removed: src/Interfaces/FormFactorXrayDB/  ext/  CondaPkg.toml
changed: Project.toml        weakdeps/extensions gone; +LinearAlgebra
         test/Project.toml   PythonCall/CondaPkg gone
         src/Interfaces/Interfaces.jl, src/Scattering/Scatterers.jl   FormFactorSourceXrayDB -> FormFactorSourceTables
         test/runtests.jl, test/test_formfactor.jl
```

---

## 2. Test coverage added beforehand

The partial-wave and three-species code had zero correctness tests. Two files
were added, taking the suite from 6645 to 7000+ assertions.

**`test/test_partialwave.jl`** — unit contracts for `PartialWave.jl`:
`compute_B_lm` against an independent naive transcription of
`B_lm = Σ_i f_i j_l(qr_i) conj(Y_lm)`; the two-channel claim verified directly
(channel 1 == `B(real.(f))`, channel 2 == `B(imag.(f))`); chunk invariance;
`cross_scatter(B,B,w) == self_scatter(B,w)`; the full error contract.

**`test/test_scatterers.jl`** — ground truth, not self-consistency. The backbone
is the **Debye formula** `I(q) = Re Σ_i Σ_j f_i conj(f_j) j_0(q r_ij)`, computed
directly over cartesian coordinates, sharing no code with `SphFuncs`/
`PartialWave` — no `Y_lm`, no packing, no chunking, no channels.

The oracle side takes **no input from the pipeline**: radii, coordinates, q grid
and form factors are hardcoded literals, and volumes are recomputed by hand as
`(4/3)πr³`. The pipeline side runs `create` → `radii`/`vols` →
`_gaussian_dummy` → `compute_B_lm` → `self_scatter`. That asymmetry is the point:
a bug in the radii lookup or the volume formula surfaces as a mismatch instead of
being inherited by both sides and cancelling.

Convergence is asserted as a **ladder** rather than at one tolerance (error must
fall at every rung), plus a deliberately-truncated `lMax = 1` case required to
disagree by >1% so the tolerance cannot pass anything. Measured for `excluded` on
a 20-atom asymmetric fixture: `lMax = 1, 2, 4, 8, 12` → `0.150, 1.60e-2, 2.03e-4, 2.82e-9, 3.46e-15`.

Two known gaps are covered separately rather than left to the oracle: the Debye
sum is invariant under rigid motion so it cannot see a bug in `_center` or the
cartesian→spherical conversion (tested directly against hand arithmetic); and
`hydration`'s positions come from `SASA.shell_points`, which cannot be
meaningfully hardcoded, so those tests check only the
`shell_points → B_lm → I(q)` leg.

## 3. Chunking in `compute_B_lm` — checked, valid

Asked whether chunking by atom is even legitimate. It is, for a precise reason:
`B_lm = Σ_i (…)` is **linear** in the atom set, so partitioning the sum is
associativity. Verified — chunked vs single-shot agrees to 1.4e-14, and
`chunk ≥ N` is bit-identical.

What *would* be invalid is chunking the intensity: `S = 4π Σ w|B|²` is quadratic,
and summing per-chunk `self_scatter` gives 57–83% errors. The code gets this
right — it never squares until the atom sum completes — but the two are one
refactor apart and nothing said so. Related: `n_chan` is computed globally before
the loop, which it must be; per-chunk it would assign different channel counts to
different chunks.

Chunking buys **no speed** (chunk=512 and no-chunking are within noise at
N=2000), and small chunks actively hurt (chunk=1 is 1.8× slower, 7× more
allocation, because `sphHarm` rebuilds its cache per call). Its value is peak
memory, which matters only at scale: N=2000/lMax=16/Q=50 peaks at ~20 MB
(irrelevant), but N=10⁵/lMax=32/Q=500 would need ~0.9 GB for `Y` and ~13 GB for
`j`. Keep it.

Two defects noted, not yet fixed:

- The docstring's `O(N·Q·lMax)` bound is wrong for `Y`, which is
  `O(N·lMax²)` — there are `(lMax+1)(lMax+2)/2` packed rows. True bound is
  `O(chunk·(Q·lMax + lMax²))`.
- `_CHUNK` has no default and is threaded through every species function, so
  each caller invents one — and anything below ~64 is strictly worse.

---

## 4. Plastic points vs Fibonacci — where progressivity is actually used

`PlasticMap` (2-D plastic/R₂ sequence, Lambert equal-area lift) is used at three
sites in `SASA.jl`: `sasa` (`n_exp = 4096`), `shell_points` (`_SHELL_SAMPLE = 256`), `_class_loop` (`_BEAD_RAY_DIRS = 64`).

Only **one** genuinely needs the sequence's progressivity: `_sasa_loop!` consumes
`pmap[1:n_occ]` for the witness pass and then `n_occ+1:n_exp` for the confirm
pass, from one array, with `_check_sasa_args` enforcing `n_occ ≤ n_exp`. That
prefix is load-bearing for **correctness, not precision** — the `cnt == 0` bail
decides burial from the prefix alone. `_prefix_thin` uses it in a weaker,
filtered form (a prefix of the occlusion-filtered subsequence). The other two
sites consume their sets whole and use no such property.

Measured spherical-cap discrepancy (4000 random caps, lower is better):


| N    | plastic | Fibonacci  | random |
| ------ | --------- | ------------ | -------- |
| 256  | 0.0239  | 0.0242     | 0.0901 |
| 4096 | 0.0037  | **0.0024** | 0.0239 |

Fibonacci is ~1.4× better at fixed N. But taking the first `k` of an N=4096 set:


| k    | plastic | Fibonacci  | plastic ‖centroid‖ | Fib ‖centroid‖ |
| ------ | --------- | ------------ | ---------------------- | ------------------ |
| 64   | 0.0749  | **0.9622** | 0.0405               | **0.9844**       |
| 1024 | 0.0094  | **0.7091** | 0.0024               | **0.7500**       |

A Fibonacci prefix is a **polar cap** — `z_i = 1 − (2i−1)/N` marches
monotonically pole to pole. With that, the witness pass would probe one cap and
be structurally blind to exposure on the opposite side, silently declaring atoms
buried. (The measured plastic prefix centroid at k=64, 0.0405, reproduces the
`0.04` already documented in `shell_points`' docstring.)

Real drawbacks, worth knowing: convergence on the actual cap-indicator integrand
is ~N^−0.74, between Monte Carlo and smooth QMC, because indicator functions have
unbounded Hardy–Krause variation — this caps how much the sequence choice can
buy. The set has a fixed lab-frame axis and is never rotated, so results are
orientation-dependent (~7e-4 relative at N=4096) and errors correlate across
atoms rather than averaging out. A deterministic per-atom rotation would fix both
at no cost while preserving the reproducibility the tests rely on.

Note CRYSOL — which this module explicitly calibrates against (`--fb 17`,
`F(17) = 1597`; `SHELL_MIN_POINTS = 55 = F(10)`) — uses Fibonacci grids. The
deviation is deliberate; the comments should say so.
