# ScatterNet.jl

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
py/FormFact_py.py      Python form-factor tiers
data/atomic_radii.sqlite3
test/                  test suites + Aqua/JET
```

## Modules

Each submodule is a `module … end` in its own file, `include`d by `src/ScatterNet.jl` in dependency order. `export` lists the public surface;  `_`-prefixed names are internal. Swappable  backends (`RadiiSource`, `FormFactorSource`) are abstract types + a documented method contract, in place of OCaml functors.

## Run

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'      # xraydb tests self-skip without the CondaPkg env
```

`PythonCall`/`CondaPkg` build the `xraydb`+`numpy` env on first `pyimport`; force
it with `julia --project=. -e 'using CondaPkg; CondaPkg.resolve()'`.
