(* Inline tests (ppx_inline_test) for Form_factor_xraydb. Both this module
and FormFact_py.py are internal to the scattering library (see
scattering.mli), so these tests live inside the library rather than as an
external test executable.

Requires python3 with numpy/xraydb importable (see requirements.txt) -
these tests genuinely call out to Python via pyml, they don't stub it. *)

module F = Form_factor_xraydb.FormFactorSourceXrayDB
module Nd = Owl_dense_ndarray_d
module Cd = Owl_dense_ndarray_z

let tol = 1e-6

let check_complex (expected : Complex.t) (actual : Complex.t) =
    Complex.norm (Complex.sub expected actual) < tol

let qvals = Nd.of_array [| 0.1; 0.2 |] [| 2 |]

(* known xraydb reference values: fe3+ at 8000 eV, q = 0.1 / 0.2 *)
let%test "known fe3+ values at 8000 eV" =
    let t = F.create 8000.0 [ "fe3+" ] qvals in
    match F.lookup t [ "fe3+" ] qvals with
    | [ (ion, arr) ] ->
        ion = "fe3+"
        && Cd.numel arr = 2
        && check_complex { re = 21.73320334; im = 3.20285267 } (Cd.get arr [| 0 |])
        && check_complex { re = 21.71503439; im = 3.20285267 } (Cd.get arr [| 1 |])
    | _ -> false

let%test "dummy ion is logged and silently dropped from lookup" =
    let t = F.create 8000.0 [ "fe3+"; "xx" ] qvals in
    Array.exists (( = ) "DUMMY   xx") (F.log t)
    && List.length (F.lookup t [ "fe3+"; "xx" ] qvals) = 1

let%test "lookup silently drops an ion that was never in the container" =
    let t = F.create 8000.0 [ "fe3+" ] qvals in
    F.lookup t [ "not_built" ] qvals = []

let%test "lookup at a q not in the original grid raises" =
    let t = F.create 8000.0 [ "fe3+" ] qvals in
    let off_grid = Nd.of_array [| 0.15 |] [| 1 |] in
    match F.lookup t [ "fe3+" ] off_grid with
    | _ -> false
    | exception Form_factor_xraydb.Form_factor_xraydb_error _ -> true

let%test "compute_form_factors: empty ions raises" =
    match Form_factor_xraydb.compute_form_factors [||] 8000.0 qvals with
    | _ -> false
    | exception Form_factor_xraydb.Form_factor_xraydb_error _ -> true

let%test "compute_form_factors: empty qvals raises" =
    let empty = Nd.of_array [||] [| 0 |] in
    match Form_factor_xraydb.compute_form_factors [| "fe3+" |] 8000.0 empty with
    | _ -> false
    | exception Form_factor_xraydb.Form_factor_xraydb_error _ -> true

let%test "compute_form_factors: energy <= 0 raises" =
    match Form_factor_xraydb.compute_form_factors [| "fe3+" |] 0.0 qvals with
    | _ -> false
    | exception Form_factor_xraydb.Form_factor_xraydb_error _ -> true

let%test "compute_form_factors: negative qval raises" =
    let bad_q = Nd.of_array [| -0.1 |] [| 1 |] in
    match Form_factor_xraydb.compute_form_factors [| "fe3+" |] 8000.0 bad_q with
    | _ -> false
    | exception Form_factor_xraydb.Form_factor_xraydb_error _ -> true

let%test "log is empty when every ion resolves fully (no dummy/f0-only)" =
    let t = F.create 8000.0 [ "fe3+" ] qvals in
    Array.length (F.log t) = 0
