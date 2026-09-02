(* Sanity checks for SphFuncs.sphHarm / SphFuncs.sphBess against known
closed-form values, plus exception-contract checks.
Prints nothing on success. On failure, prints one line per failing
check and exits with a non-zero status. *)
open Printf

module S = Scattering.SphFuncs

let tol = 1e-9

let failures = ref 0

let fail fmt = ksprintf (fun s -> incr failures; print_endline ("FAIL " ^ s)) fmt

let check_float name expected actual =
    let d = Float.abs (expected -. actual) in
    if d > tol then
        fail "%s: expected %.15f, got %.15f (diff %.3e)" name expected actual d

let check_complex name (expected : Complex.t) (actual : Complex.t) =
    let d = Complex.norm (Complex.sub expected actual) in
    if d > tol then
        fail "%s: expected (%.15f + %.15fi), got (%.15f + %.15fi) (diff %.3e)"
            name expected.re expected.im actual.re actual.im d

let check_raises name f =
    match f () with
    | _ -> fail "%s: expected an exception, none was raised" name
    | exception _ -> ()

let pi = Float.pi

(* ---- sphHarm: known closed-form Y_l^m values ---- *)

let y00 = 1.0 /. (2.0 *. sqrt pi)
let y10 theta = sqrt (3.0 /. (4.0 *. pi)) *. cos theta
let y11 theta phi = Complex.polar (-.(sqrt (3.0 /. (8.0 *. pi))) *. sin theta) phi

let idx l m = (l * (l + 1) / 2) + m

let test_sphHarm_known_values () =
    let cases =
        [ (pi /. 2.0, 0.0); (0.0, 0.7) ]
    in
    List.iter (fun (theta_v, phi_v) ->
        let theta = Owl_dense_ndarray_d.of_array [| theta_v |] [| 1 |] in
        let phi = Owl_dense_ndarray_d.of_array [| phi_v |] [| 1 |] in
        let y = S.sphHarm 1 theta phi in
        let tag = sprintf "sphHarm(theta=%.6f, phi=%.6f)" theta_v phi_v in
        check_complex (tag ^ " Y_0^0") { re = y00; im = 0.0 } (
            Owl_dense_ndarray_z.get y [| idx 0 0; 0 |]
        );
        check_complex (tag ^ " Y_1^0") { re = y10 theta_v; im = 0.0 } (
            Owl_dense_ndarray_z.get y [| idx 1 0; 0 |]
        );
        check_complex (tag ^ " Y_1^1") (y11 theta_v phi_v) (
            Owl_dense_ndarray_z.get y [| idx 1 1; 0 |]
        )
    ) cases

(* ---- sphHarm: exception contract ---- *)

let test_sphHarm_exceptions () =
    let theta1 = Owl_dense_ndarray_d.of_array [| 1.0 |] [| 1 |] in
    let phi1 = Owl_dense_ndarray_d.of_array [| 1.0 |] [| 1 |] in
    let theta2 = Owl_dense_ndarray_d.of_array [| 1.0; 2.0 |] [| 2 |] in
    let empty = Owl_dense_ndarray_d.of_array [||] [| 0 |] in
    let theta_2d = Owl_dense_ndarray_d.zeros [| 2; 2 |] in
    check_raises "sphHarm negative lMax" (fun () -> S.sphHarm (-1) theta1 phi1);
    check_raises "sphHarm mismatched lengths" (fun () -> S.sphHarm 2 theta1 theta2);
    check_raises "sphHarm empty theta/phi" (fun () -> S.sphHarm 2 empty empty);
    check_raises "sphHarm non-1D theta" (fun () -> S.sphHarm 2 theta_2d phi1)

(* ---- sphBess: known closed-form j_l values ---- *)

let j0 x = if x = 0.0 then 1.0 else sin x /. x
let j1 x = if x = 0.0 then 0.0 else (sin x /. (x *. x)) -. (cos x /. x)

let test_sphBess_known_values () =
    let cases = [ 1.0; 0.0 ] in
    List.iter
        (fun x ->
            let r = Owl_dense_ndarray_d.of_array [| x |] [| 1 |] in
            let q = Owl_dense_ndarray_d.of_array [| 1.0 |] [| 1 |] in
            let j = S.sphBess r q 1 in
            let tag = sprintf "sphBess(x=%.6f)" x in
            check_float (tag ^ " j_0") (j0 x) (Owl_dense_ndarray_d.get j [| 0; 0; 0 |]);
            check_float (tag ^ " j_1") (j1 x) (Owl_dense_ndarray_d.get j [| 1; 0; 0 |]))
        cases

(* ---- sphBess: exception contract ---- *)

let test_sphBess_exceptions () =
    let r1 = Owl_dense_ndarray_d.of_array [| 1.0 |] [| 1 |] in
    let q1 = Owl_dense_ndarray_d.of_array [| 1.0 |] [| 1 |] in
    let neg = Owl_dense_ndarray_d.of_array [| -1.0 |] [| 1 |] in
    let empty = Owl_dense_ndarray_d.of_array [||] [| 0 |] in
    check_raises "sphBess empty radii" (fun () -> S.sphBess empty q1 1);
    check_raises "sphBess empty q grid" (fun () -> S.sphBess r1 empty 1);
    check_raises "sphBess negative radii" (fun () -> S.sphBess neg q1 1);
    check_raises "sphBess negative q" (fun () -> S.sphBess r1 neg 1);
    check_raises "sphBess negative lMax" (fun () -> S.sphBess r1 q1 (-1))

let () =
    test_sphHarm_known_values ();
    test_sphHarm_exceptions ();
    test_sphBess_known_values ();
    test_sphBess_exceptions ();
    if !failures > 0 then exit 1
