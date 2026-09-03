(* Checks for Molecule's create/theta/phi/r/vols against known values.
Prints nothing on success; on failure prints one line per failing check
and exits non-zero. *)

module M = Molecule

let failures = ref 0

let fail fmt = Printf.ksprintf (fun s -> incr failures; print_endline ("FAIL " ^ s)) fmt

let check_shape name expected (arr : Owl_dense_ndarray_d.arr) =
    let got = Owl_dense_ndarray_d.shape arr in
    if got <> expected then
        fail "%s: expected shape [%s], got [%s]" name
            (String.concat ";" (Array.to_list (Array.map string_of_int expected)))
            (String.concat ";" (Array.to_list (Array.map string_of_int got)))

let check_float name expected actual =
    let d = Float.abs (expected -. actual) in
    if d > 1e-9 then fail "%s: expected %f, got %f (diff %.3e)" name expected actual d

let test_two_atoms_on_x_axis () =
    (* atoms at (1,0,0) and (-1,0,0): centroid is already the origin *)
    let m = M.create "test" [| "h"; "h" |] [| (1.0, 0.0, 0.0); (-1.0, 0.0, 0.0) |] in

    let r = M.r m in
    check_shape "r" [| 2 |] r;
    check_float "atom0 r" 1.0 (Owl_dense_ndarray_d.get r [| 0 |]);
    check_float "atom1 r" 1.0 (Owl_dense_ndarray_d.get r [| 1 |]);

    let theta = M.theta m and phi = M.phi m in
    check_shape "theta" [| 2 |] theta;
    check_shape "phi" [| 2 |] phi;
    (* atom0 at (1,0,0): theta = acos(0/1) = pi/2, phi = atan2(0,1) = 0 *)
    check_float "atom0 theta" (Float.pi /. 2.0) (Owl_dense_ndarray_d.get theta [| 0 |]);
    check_float "atom0 phi" 0.0 (Owl_dense_ndarray_d.get phi [| 0 |]);
    (* atom1 at (-1,0,0): theta = pi/2, phi = atan2(0,-1) = pi *)
    check_float "atom1 theta" (Float.pi /. 2.0) (Owl_dense_ndarray_d.get theta [| 1 |]);
    check_float "atom1 phi" Float.pi (Owl_dense_ndarray_d.get phi [| 1 |])

let test_centering_shifts_to_centroid () =
    (* atoms at (0,0,0) and (2,0,0): centroid is (1,0,0), so both atoms end
    up 1 unit from it. If centering hadn't happened, r would be [0, 2]
    instead of [1, 1]. *)
    let m = M.create "test" [| "h"; "h" |] [| (0.0, 0.0, 0.0); (2.0, 0.0, 0.0) |] in
    let r = M.r m in
    check_float "atom0 r after centering" 1.0 (Owl_dense_ndarray_d.get r [| 0 |]);
    check_float "atom1 r after centering" 1.0 (Owl_dense_ndarray_d.get r [| 1 |])

let test_single_atom_r_zero_no_nan () =
    (* a single atom is its own centroid, so r = 0; must not produce NaN *)
    let m = M.create "test" [| "h" |] [| (5.0, 5.0, 5.0) |] in
    check_float "single atom r" 0.0 (Owl_dense_ndarray_d.get (M.r m) [| 0 |]);
    let theta = Owl_dense_ndarray_d.get (M.theta m) [| 0 |] in
    if Float.is_nan theta then fail "single atom theta: got NaN"

let test_empty_coords_raises () =
    match M.create "empty" [||] [||] with
    | _ -> fail "empty coords: expected Molecule_error, none was raised"
    | exception M.Molecule_error _ -> ()

let test_vols_known_elements () =
    let m =
        M.create "test" [| "fe"; "o"; "rn" |]
            [| (0.0, 0.0, 0.0); (1.0, 0.0, 0.0); (2.0, 0.0, 0.0) |]
    in
    let v = M.vols m in
    check_shape "vols" [| 3 |] v;
    let expected_volume (r : float) : float = (4.0 /. 3.0) *. Float.pi *. (r ** 3.0) in
    check_float "fe volume" (expected_volume 1.274) (Owl_dense_ndarray_d.get v [| 0 |]);
    check_float "rn volume" (expected_volume 2.24) (Owl_dense_ndarray_d.get v [| 2 |])

let test_vols_unknown_element_raises () =
    (* elements are only looked up when vols is forced, not at create time *)
    let m = M.create "test" [| "zzzz" |] [| (0.0, 0.0, 0.0) |] in
    match M.vols m with
    | _ -> fail "unknown element: expected Molecule_error, none was raised"
    | exception M.Molecule_error _ -> ()

let test_create_computes_r_theta_phi_vols () =
    let m =
        M.create "test" [| "o"; "h"; "h" |]
            [| (0.0, 0.0, 0.0); (1.0, 0.0, 0.0); (0.0, 1.0, 0.0) |]
    in
    let r = M.r m in
    check_shape "r" [| 3 |] r;
    check_float "r atom0" 0.4714045208 (Owl_dense_ndarray_d.get r [| 0 |]);
    check_shape "theta" [| 3 |] (M.theta m);
    check_shape "phi" [| 3 |] (M.phi m);
    check_shape "vols" [| 3 |] (M.vols m);
    (* forcing again must return the same value *)
    let r2 = M.r m in
    check_float "r stable on repeat access" 
    (Owl_dense_ndarray_d.get r [| 0 |]) 
    (Owl_dense_ndarray_d.get r2 [| 0 |])

let test_create_length_mismatch_raises () =
    match M.create "bad" [| "o"; "h" |] [| (0.0, 0.0, 0.0) |] with
    | _ -> fail "length mismatch: expected Molecule_error, none was raised"
    | exception M.Molecule_error _ -> ()

let () =
    test_two_atoms_on_x_axis ();
    test_centering_shifts_to_centroid ();
    test_single_atom_r_zero_no_nan ();
    test_empty_coords_raises ();
    test_vols_known_elements ();
    test_vols_unknown_element_raises ();
    test_create_computes_r_theta_phi_vols ();
    test_create_length_mismatch_raises ();
    if !failures > 0 then begin
        Printf.printf "%d check(s) failed\n" !failures;
        exit 1
    end else
        print_endline "all checks passed"
