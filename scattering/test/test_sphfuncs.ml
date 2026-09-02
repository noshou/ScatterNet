(* Inline tests (ppx_inline_test) for SphFuncs: known closed-form values,
plus exception contracts. SphFuncs is internal to the scattering library
(see scattering.mli), so these tests live inside the library rather than
as an external test executable. *)

let tol = 1e-9

let check_float expected actual = Float.abs (expected -. actual) < tol

let check_complex (expected : Complex.t) (actual : Complex.t) =
    Complex.norm (Complex.sub expected actual) < tol

let check_raises f = match f () with _ -> false | exception _ -> true

let pi = Float.pi

(* ---- sphHarm: known closed-form Y_l^m values ---- *)

let y00 = 1.0 /. (2.0 *. sqrt pi)
let y10 theta = sqrt (3.0 /. (4.0 *. pi)) *. cos theta
let y11 theta phi = Complex.polar (-.(sqrt (3.0 /. (8.0 *. pi))) *. sin theta) phi

let idx l m = (l * (l + 1) / 2) + m

let check_sphHarm_case (theta_v, phi_v) =
    let theta = Owl_dense_ndarray_d.of_array [| theta_v |] [| 1 |] in
    let phi = Owl_dense_ndarray_d.of_array [| phi_v |] [| 1 |] in
    let y = SphFuncs.sphHarm 1 theta phi in
    check_complex { re = y00; im = 0.0 } (Owl_dense_ndarray_z.get y [| idx 0 0; 0 |])
    && check_complex { re = y10 theta_v; im = 0.0 }
        (Owl_dense_ndarray_z.get y [| idx 1 0; 0 |])
    && check_complex (y11 theta_v phi_v) (Owl_dense_ndarray_z.get y [| idx 1 1; 0 |])

let%test "sphHarm known values" =
    List.for_all check_sphHarm_case [ (pi /. 2.0, 0.0); (0.0, 0.7) ]

(* ---- sphHarm: exception contract ---- *)

let%test "sphHarm exception contract" =
    let theta1 = Owl_dense_ndarray_d.of_array [| 1.0 |] [| 1 |] in
    let phi1 = Owl_dense_ndarray_d.of_array [| 1.0 |] [| 1 |] in
    let theta2 = Owl_dense_ndarray_d.of_array [| 1.0; 2.0 |] [| 2 |] in
    let empty = Owl_dense_ndarray_d.of_array [||] [| 0 |] in
    let theta_2d = Owl_dense_ndarray_d.zeros [| 2; 2 |] in
    check_raises (fun () -> SphFuncs.sphHarm (-1) theta1 phi1)
    && check_raises (fun () -> SphFuncs.sphHarm 2 theta1 theta2)
    && check_raises (fun () -> SphFuncs.sphHarm 2 empty empty)
    && check_raises (fun () -> SphFuncs.sphHarm 2 theta_2d phi1)

(* ---- sphBess: known closed-form j_l values ---- *)

let j0 x = if x = 0.0 then 1.0 else sin x /. x
let j1 x = if x = 0.0 then 0.0 else (sin x /. (x *. x)) -. (cos x /. x)

let check_sphBess_case x =
    let r = Owl_dense_ndarray_d.of_array [| x |] [| 1 |] in
    let q = Owl_dense_ndarray_d.of_array [| 1.0 |] [| 1 |] in
    let j = SphFuncs.sphBess r q 1 in
    check_float (j0 x) (Owl_dense_ndarray_d.get j [| 0; 0; 0 |])
    && check_float (j1 x) (Owl_dense_ndarray_d.get j [| 1; 0; 0 |])

let%test "sphBess known values" = List.for_all check_sphBess_case [ 1.0; 0.0 ]

(* ---- sphBess: exception contract ---- *)

let%test "sphBess exception contract" =
    let r1 = Owl_dense_ndarray_d.of_array [| 1.0 |] [| 1 |] in
    let q1 = Owl_dense_ndarray_d.of_array [| 1.0 |] [| 1 |] in
    let neg = Owl_dense_ndarray_d.of_array [| -1.0 |] [| 1 |] in
    let empty = Owl_dense_ndarray_d.of_array [||] [| 0 |] in
    check_raises (fun () -> SphFuncs.sphBess empty q1 1)
    && check_raises (fun () -> SphFuncs.sphBess r1 empty 1)
    && check_raises (fun () -> SphFuncs.sphBess neg q1 1)
    && check_raises (fun () -> SphFuncs.sphBess r1 neg 1)
    && check_raises (fun () -> SphFuncs.sphBess r1 q1 (-1))
