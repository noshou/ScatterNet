# Why migrate ScatterNet from OCaml to Julia

Status: assessment. Written to be argued with.

## TL;DR

- **The forward model must be differentiable in the same language the sampler runs in.** Stage 1 of the planned pipeline is NUTS/HMC over spherical-harmonic coefficients, and HMC needs the gradient of the likelihood *through* the physics. Julia gives that in-language (`ForwardDiff`/`Zygote`/`Enzyme`). OCaml does not — the current plan reimplements the physics in Python/JAX. Julia removes the second implementation.
- **Both planned stages move in-process.** `AdvancedHMC.jl`/`Turing.jl` (stage 1) and `Lux.jl`+`Zygote` (stage 2) run in one process against one forward model. The OCaml plan ran both in Python via `pyml` because OCaml has no usable NUTS and no ML stack.
- **The native-C surface mostly disappears.** GSL, the `owl` source patch + `opam pin`, `conf-openblas`, and `ocamllex` are replaced by registered, pure-Julia, permissive packages. `libsqlite3` stays, bundled inside `SQLite.jl`.
- **One Python boundary remains — the same one OCaml had.** `xraydb` has no native-Julia equivalent (`NeXLCore.jl` lacks `f0(q)`; `FFAST.jl` is archived and raw). `PythonCall`+`CondaPkg` replace `pyml`. Parity, not regression.
- **aarch64/Asahi is a first-class Julia target.** Tier-1, pure-Julia math packages, no 16K-page concern. Not true for the OCaml native stack, and actively hostile for .NET.

## What the migration deletes

| OCaml pain | Why it exists | Julia replacement |
|---|---|---|
| `_patches/owl-1.2-exponpow-args.patch` + `opam pin` in `setup.sh` | `owl` 1.2 doesn't compile as released | No `owl`; arrays are stdlib, distributions from `Distributions.jl` |
| `setup.sh` (pin owl → `opam install` → `pip install`) | build can't bootstrap itself | `Pkg.instantiate()` from `Manifest.toml`; `CondaPkg` for the Python side |
| `Lexer.mll` + `Parser.ml` + `Toks.ml` (`ocamllex`) for `"fe3+"` | a lexer generator for a 2-token grammar | one regex, ~20 lines, no codegen |
| `helpers/Cache.ml` (`Lazy` + `Mutex` per derived field) | OCaml 5 `Lazy.force` races across domains | eager compute, or one `@lock`ed cell; `Threads.@spawn`/`fetch` idioms |
| `Helpers.repo_root` (parse `Sys.executable_name`, stop at `_build`) | find the data file regardless of `dune` cwd | `pkgdir(@__MODULE__)` |
| `Domain.DLS`-keyed sqlite handle per domain | a `~mutex:NO` handle isn't domain-safe | load 542 rows once into `const Dict`s at `__init__`; no live handle |
| `ppx_inline_test` + `(preprocess (pps …))` per stanza | inline `let%test` needs a PPX | `Test` stdlib + `@testset`, run by `Pkg.test()` |
| system **GSL**, **OpenBLAS** (`conf-openblas`), **libpython**, **libsqlite3** | four separate prerequisites | GSL → `SphericalHarmonics.jl`/`LegendrePolynomials.jl`/`Bessels.jl`; OpenBLAS ships with Julia; libsqlite3 inside `SQLite.jl`; only libpython remains, `CondaPkg`-provisioned |
| `opam` + `dune` + generated `.opam` | two-layer build, non-lockfile | `Project.toml` + `Manifest.toml` (the manifest *is* the lockfile) |

## The pipeline argument

Stage 1 (NUTS/HMC) needs `∇(log-posterior)` w.r.t. the SH coefficients, and the log-posterior contains the forward model. So HMC needs the gradient *through* the spherical-harmonic sum, the Bessel-weighted radial integral, the form-factor-weighted atomic sum.

- **Julia:** that forward model is ordinary code; `ForwardDiff`/`Zygote`/`Enzyme` differentiate it; `AdvancedHMC.jl` consumes the gradient. One implementation serves both.
- **OCaml:** no autodiff through `owl`+GSL, no NUTS. The plan's answer is a *second* forward model in Python/JAX, kept in sync by hand — a permanent tax.
- **F#/.NET:** no mature NUTS; `TorchSharp`/libtorch has no first-party arm64-Linux native package. Stage 2 goes back to Python and stage 1 has nothing to stand on.
- **Haskell:** `ad` won't reach through an FFI numerics layer; no production NUTS.

Julia is the only option here where the forward model, sampler, and VAE live in one language and process, differentiable end to end. That is the reason to move.

## aarch64 / Asahi

Julia is Tier-1 on `aarch64-linux`; every replacement math package is pure Julia, so the 16K-page kernel never enters the picture and there is no per-arch native build to fail. The OCaml native surface (GSL + OpenBLAS + libpython + libsqlite3 + patched `owl`) is five separate arch/page-size risks, and `owl` needs manual intervention before arch even matters.

## What you give up (honest)

- **Compiler-enforced encapsulation.** OCaml's `.mli` is a wall; Julia's `export`/`_`-prefix is a convention. `molecule.mli`, `RadSrc.RadiiSource` as a checked signature — no Julia equivalent.
- **No functors / ML module system.** `module Radii : RadSrc.RadiiSource = …` becomes a duck-typed contract: documented and tested, not compiler-verified.
- **Sum types without exhaustiveness checking.** OCaml warns on a non-exhaustive `match`; the radii fallback chain is currently compiler-checked. In Julia it's `if`/`nothing` and nothing tells you a case is missing.
- **Dynamic typing + a performance discipline.** Type-unstable code runs 10–100× slower with no error. `@code_warntype` and `JET.jl` become routine, not optional. This port is written for stability (concrete fields, typed accumulators, function-barrier `Lazy{T}`, `@inferred` in tests) — that discipline is now permanent.
- **`JET.jl` couples to compiler internals** and can break on a Julia minor bump.
- **Compile latency.** First call in a session pays JIT + package load — seconds, more with Turing. `dune exec` on a built binary starts instantly. Amortizes to nothing over a sampler run; a friction tax for interactive use.

## Net recommendation

Migrate. The port is ~1,500 lines and its hardest parts become library calls. Against that one-time cost: HMC-with-gradients-through-the-physics is native in Julia and impossible in OCaml without a duplicate forward model; the migration deletes the `owl` patch, `opam pin`, `setup.sh`, the `ocamllex` trio, the `Lazy`+`Mutex` cache, the `repo_root` hack, the `Domain.DLS` dance, `ppx_inline_test`, and three of four system C libraries; aarch64/Asahi goes from a risk surface to Tier-1. You accept a weaker type system and a type-stability discipline. For a project about to grow a sampler and a neural net, that trade is worth making. Keep Python for `xraydb`.
