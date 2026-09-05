# ScatterNet.jl

## Layout

```
src/
  ScatterNet.jl          thin aggregator (includes Interfaces + AtomicRadii + the two groups)
  Interfaces.jl          RadiiSource / FormFactorSource markers (shared)
  AtomicRadii/
    AtomicRadii.jl       ion parsing + table (loaded once) + fallback chain
    atomic_radii.sqlite3
  Molecule/
    Molecule.jl          module Molecule
    Cache.jl             Lazy{T} + force
    Molecules.jl         create + r/theta/phi/coords/radii/vols/elms/name
    SASA.jl              Shrake-Rupley SASA + hydration-shell patch geometry
  Scattering/
    Scattering.jl        module Scattering
    SphFuncs.jl          sphHarm, sphBess          (SphericalHarmonics.jl, Bessels.jl)
    FormFactorXrayDB.jl  compute_form_factors      (xraydb via PythonCall)
py/FormFact_py.py        Python form-factor tiers
test/                    test suites + Aqua/JET
```

## Modules

`src/ScatterNet.jl` includes `Interfaces`, `AtomicRadii`, then two grouped
submodules, each its own folder with a parent `module … end` that `include`s
its files in dependency order: `Molecule` (`Cache`, `Molecules`, `SASA`) and
`Scattering` (`SphFuncs`, `FormFactorXrayDB`). `AtomicRadii` is a sibling of
`Molecule`, not nested inside it: it only implements `Interfaces.RadiiSource`,
the abstract backend `Molecule.Molecules.create` consults, so it has no
dependency on `Molecule` in the other direction and ships its own bundled
`atomic_radii.sqlite3` next to it. Reach a leaf as `ScatterNet.Molecule.SASA`
etc. `export` lists the public surface; `_`-prefixed names are internal.
Swappable backends (`RadiiSource`, `FormFactorSource`) are abstract types + a
documented method contract, in place of OCaml functors.

## Run

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'      # xraydb tests self-skip without the CondaPkg env
```

`PythonCall`/`CondaPkg` build the `xraydb`+`numpy` env on first `pyimport`; force
it with `julia --project=. -e 'using CondaPkg; CondaPkg.resolve()'`.
