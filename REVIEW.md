# OCaml → Julia migration review for ScatterNet

Date: 2026-09-02. Reviews the OCaml tree at commit `b2c69e4` (preserved on the
`master` branch) plus the planned pipeline in `CLAUDE.md`. An F#/.NET port was
also prototyped and evaluated before Julia was chosen; this review is
Julia-specific.

## TL;DR

Julia is the **right target for this project**, and the reason is the *planned
pipeline*, not the current code:

- The forward model is ~1,500 lines of arithmetic wrapping GSL and `xraydb`.
  Any language ports it. On its own it doesn't justify a migration.
- The pipeline (`CLAUDE.md`) is: special functions → forward model →
  **gradient-based NUTS** → **conditional VAE**. Julia is one of only two
  ecosystems (the other is Python/JAX) where all four are first-class *in one
  language*: `SphericalHarmonics.jl`/`Bessels.jl`, `AdvancedHMC.jl`/`Turing.jl`
  (best-in-class NUTS, no C++ compile), `Lux.jl` + `Zygote`/`Enzyme`.
- Critically: stage-1 HMC needs `∇` of the likelihood through the forward
  model. In Julia the forward model is **differentiable in-language**, so no
  reimplementation in a separate autodiff framework. OCaml, F#, and Haskell all
  fail this — they'd force a second implementation of the physics in JAX.
- The OCaml design already routed the sampler and VAE through Python (`pyml`).
  Julia **removes that boundary** instead of preserving it. After the port the
  only non-Julia code is `FormFact_py.py`.

The cost is real but bounded — see "What you give up".

## What the port simplifies or deletes

| OCaml pain | Julia |
|---|---|
| `_patches/owl-1.2-exponpow-args.patch` + `opam pin` of a hand-patched Owl clone + `setup.sh` step 1 | **gone** — no Owl. The array ops here are elementwise on (≤3, n); plain `Array`s. |
| `ocamllex` `Lexer.mll` + `Toks.ml` + `Parser.ml` (≈90 LOC, codegen step) for strings like `"fe3+"` | one regex, ≈15 LOC |
| `helpers/Cache.ml` — `Lazy` + `Mutex` for multi-domain safety | `Cache.Lazy{T}`, ~20 LOC |
| `Db.ml` + `sqlite3` opam pkg + `Domain.DLS` connection-per-domain dance | 542 static rows read once at `__init__` into `const Dict`s; SQLite.jl only as the loader (bundled native lib, no system libsqlite3) |
| GSL bindings (`gsl` opam pkg, system libgsl, GPL) | `SphericalHarmonics.jl` + `Bessels.jl` + `LegendrePolynomials.jl` — registered, pure-Julia, permissive, aarch64-clean. `computeYlm`'s default is *already* the CS-phased QM `Y_l^m`. |
| `pyml` — embed libpython, hand-manage the GIL, `sys.path` hack via `Sys.executable_name` | `PythonCall` + `CondaPkg` — env auto-provisioned, GIL automatic. `pkgdir()` for the path. |
| `opam` + `dune` + generated `.opam` + 5 `dune` files + the `(include_subdirs no)` workaround | one `Project.toml`, one `test/Project.toml` |
| `ppx_inline_test` + full Jane Street ppx stack + a hand-rolled assertion runner in `test_molecule.ml` | `Test` + `Aqua` + `JET` |
| aarch64/Asahi: OpenBLAS via `conf-openblas`, libgsl, libsqlite3, libpython all must exist as system libs | Julia is Tier-1 on `aarch64-linux`; the math packages are pure Julia; SQLite/Python natives are bundled by their JLLs |

**Pipeline-level:** the OCaml plan is OCaml + Python (sampler) + Python (VAE) —
three toolchains, two languages, and a forward model that must be written twice
(once in OCaml, once in JAX for stage-1 gradients). The Julia plan is one
language, one autodiff stack, one process, forward model written once.

## What you give up (be honest)

1. **Sum types + exhaustive pattern matching.** The single biggest loss coming
   from ML. `Union{Ion,Nothing}` etc. work and union-split efficiently, but
   there's no compiler exhaustiveness check on a `match`. The ion parser and the
   4-tier fallback chain lose that safety net.
2. **No functors / ML module system.** `RadSrc.RadiiSource` and
   `FFSrc.FormFactorSource` become abstract types + a *documented* method
   contract (`Interfaces.jl`) — nothing enforces that an implementer provides
   `ff_log`/`ff_lookup` until you call them.
3. **Dynamic typing + the type-stability tax.** Julia code is fast *if*
   type-stable and can silently degrade 10–100× when it isn't. This port is
   written for stability (concrete struct fields, typed accumulators, a
   function-barrier `Lazy{T}`, `@inferred` assertions in the tests) — but it's a
   discipline you now carry forever, and `@code_warntype`/`JET` become part of
   the workflow. OCaml/F#/Haskell never have this failure mode.
4. **`.mli` / signature-file enforcement.** Julia docstrings + `export` lists
   are convention, not a checked interface. `molecule.mli` was a compiler-checked
   contract; `Molecules.jl`'s accessor set is not.
5. **Compile latency** (time-to-first-run) and a heavier first-call cost than a
   native OCaml binary.
6. **`ppx_inline_test`'s "tests next to the code, stripped from release"** — gone;
   tests live in `test/`. (The OCaml repo was already half-migrating to this.)
7. **JET is coupled to the Julia compiler internals** and its supported version
   moves with each Julia minor — the static-analysis check needs version
   guarding.

## Numerical-fidelity notes

- `SphericalHarmonics.computeYlm` default = QM complex `Y_l^m`, ∫|Y|²dΩ = 1,
  **Condon–Shortley phase on** — verified against the package's own tests
  (`Y_1^1 = -√(3/8π) sinθ e^{iφ}`, negative). Matches GSL
  `legendre_sphPlm · e^{imφ}` for m ≥ 0. The port keeps the OCaml's m ≥ 0
  triangular packing (`idx = l*(l+1)/2 + m`) so `test_sphfuncs` holds.
- `LegendrePolynomials.Plm(x,l,m; norm=Val(:normalized)) / √(2π)` is *exactly*
  GSL `legendre_sphPlm(l,m,x)` for m ≥ 0 (both carry the CS phase). Used for the
  bare-`P̄_l^m` accessor.
- `Bessels.sphericalbesselj` uses a small-argument power series in the `l > x`
  regime (the SAXS regime), so it's stable where naive upward recurrence isn't —
  this is the case the OCaml hand-rolled a Miller recurrence for. `j_l(0)`
  returns 1 for l = 0, 0 otherwise. No batched all-`l` call; loop per `l`.
- `xraydb` path: numbers are unchanged (same `FormFact_py.py`, same Chantler
  `f1_chantler`/`f2_chantler` + Waasmaier–Kirfel `f0`). Regression constants
  (`fe3+` @ 8000 eV → `21.73320334 + 3.20285267i` at q=0.1) carried into
  `test_formfactor.jl`. **No native Julia replacement exists** — NeXLCore.jl has
  no `f0(q)` and only mass-absorption coefficients; FFAST.jl gives raw
  (uncorrected) Chantler columns and is archived. Keeping Python here is correct.
- The exact-float q-grid hash lookup in `ff_lookup` is reproduced bug-for-bug
  (1 ULP off after any round-trip raises). Flagged for a deliberate keep-vs-fix.

## Bugs found in the OCaml on the way through

Independent of the migration; worth fixing in either language:

- `setup.sh` step 3 and the error string at `FormFactorXrayDB.ml:38` use a
  stale pre-refactor path (`scattering/form_factor_xraydb/…`).
- `_patches/README.md` has a stale path, and the live `opam pin` points at
  `~/.opam-patches/owl-1.2` — divergent from what `setup.sh` writes, undocumented.
- **Negative ionic radii not filtered.** `ionic_radii` has `h1+` (−38 pm),
  `n5+` (−10.4 pm), `c4+` (−8 pm). `DATA.md` says treat `radius ≤ 0` as absent;
  `Db.ion_radius` returns them and `compute_vols` would cube them into a
  negative volume. **Fixed in this port** (dropped at load).
- `Db.ml` dead code: `close_db`, `element_radii_batch`, `charges_for` are never
  called.
- `sphBess`'s comment claims one GSL call yields every `l`; the code calls
  `bessel_jl` per `l`. (The Julia version also loops per `l` — `Bessels.jl` has
  no batched form — so the comment would now be accurate if kept.)
- `scattering.ml` / `.mli` are empty — the library exposes nothing.

## Recommended plan

1. **Target Julia 1.10 LTS** (the aarch64/16K-page story is clean — Julia is
   Tier-1 there, unlike .NET).
2. **Port order** (done in this scaffold, low → high risk):
   `Cache`/`Interfaces` → `AtomicRadii` → `Molecules` → `SphFuncs` (validate
   against a GSL/scipy grid — the OCaml tests only cover l ≤ 1) →
   `FormFactorXrayDB`.
3. **Keep the Python boundary only at `xraydb`.** Move the sampler and VAE
   *into* Julia (`AdvancedHMC.jl`, `Lux.jl`) rather than to Python — this is the
   whole point, and `CLAUDE.md` is updated to say so.
4. **Adopt `JET` + `@inferred` in CI from day one** — type stability is now a
   correctness-adjacent property, not a nicety.
5. **Fix the bugs above** with the ported tests (plus a new high-`l` SphFuncs
   oracle grid) as the acceptance bar.
6. Consider `Unitful` / units-of-measure on the physical quantities (Å, eV,
   Å⁻¹) while rewriting signatures — Julia supports it well and scattering code
   benefits.
