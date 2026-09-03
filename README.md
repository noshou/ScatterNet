# ScatterNet.jl

Forward model for X-ray / neutron scattering shape reconstruction. Ported from
an OCaml original (preserved on the `master` branch); rationale in
[`docs/WHY_JULIA.md`](docs/WHY_JULIA.md), trade-offs in [`REVIEW.md`](REVIEW.md).
Planned pipeline on top: [`CLAUDE.md`](CLAUDE.md).

## Layout

```
src/
  ScatterNet.jl        top module (includes + re-exports)
  Cache.jl             Lazy{T} + force
  Interfaces.jl        RadiiSource / FormFactorSource markers
  AtomicRadii.jl       ion parsing + 542-row table (loaded once) + fallback chain
  SphFuncs.jl          sphHarm, sphBess           (SphericalHarmonics.jl, Bessels.jl)
  Molecules.jl         create + r/theta/phi/coords/radii/vols/elms/name
  FormFactorXrayDB.jl  compute_form_factors       (xraydb via PythonCall)
py/FormFact_py.py      Python form-factor tiers (unchanged)
data/atomic_radii.sqlite3
test/                  ports of the five original test suites + Aqua/JET
```

## Modules

Each submodule is a `module … end` in its own file, `include`d by
`src/ScatterNet.jl` in dependency order. `export` lists the public surface;
`_`-prefixed names are internal (not enforced — Julia has no `.mli`). Swappable
backends (`RadiiSource`, `FormFactorSource`) are abstract types + a documented
method contract, in place of OCaml functors.

## Run

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'      # xraydb tests self-skip without the CondaPkg env
```

`PythonCall`/`CondaPkg` build the `xraydb`+`numpy` env on first `pyimport`; force
it with `julia --project=. -e 'using CondaPkg; CondaPkg.resolve()'`.

## Not verified here (no Julia toolchain on the authoring machine)

- Allocation-freeness of the `sphHarm` cache loop after warmup (`@allocated`).
- Exact `DBInterface`/`Tables` row-accessor types in `AtomicRadii._load!`.
- `@inferred` / `JET` assertions in `test/test_quality.jl` encode the intended
  type stability; run them to confirm.
