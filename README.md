# ScatterNet.jl

## Layout

```
src/
  ScatterNet.jl          thin aggregator (includes Interfaces + the two groups)
  Interfaces.jl          RadiiSource / FormFactorSource markers (shared)
  Structure/
    Structure.jl         module Structure
    Cache.jl             Lazy{T} + force
    AtomicRadii.jl       ion parsing + table (loaded once) + fallback chain
    Molecules.jl         create + r/theta/phi/coords/radii/vols/elms/name
  Scattering/
    Scattering.jl        module Scattering
    SphFuncs.jl          sphHarm, sphBess          (SphericalHarmonics.jl, Bessels.jl)
    FormFactorXrayDB.jl  compute_form_factors      (xraydb via PythonCall)
py/FormFact_py.py        Python form-factor tiers
data/atomic_radii.sqlite3
test/                    test suites + Aqua/JET
```

## Modules

`src/ScatterNet.jl` includes `Interfaces` then two grouped submodules, each its
own folder with a parent `module … end` that `include`s its files in dependency
order: `Structure` (`Cache`, `AtomicRadii`, `Molecules`) and `Scattering`
(`SphFuncs`, `FormFactorXrayDB`). Reach a leaf as
`ScatterNet.Structure.AtomicRadii` etc. `export` lists the public surface;
`_`-prefixed names are internal. Swappable backends (`RadiiSource`,
`FormFactorSource`) are abstract types + a documented method contract, in place
of OCaml functors.

## Run

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'      # xraydb tests self-skip without the CondaPkg env
```

`PythonCall`/`CondaPkg` build the `xraydb`+`numpy` env on first `pyimport`; force
it with `julia --project=. -e 'using CondaPkg; CondaPkg.resolve()'`.
