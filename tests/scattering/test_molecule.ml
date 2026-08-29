(* Checks for Scattering.Molecule's center_coords/radii/angles against
known values. Prints nothing on success; on failure prints one line per
failing check and exits non-zero. *)

module Nd = Owl_base_dense_ndarray_d
module M = Scattering.Molecule

let failures = ref 0

let fail fmt = Printf.ksprintf (fun s -> incr failures; print_endline ("FAIL " ^ s)) fmt

let check_shape name expected (arr : Nd.arr) =
    let got = Nd.shape arr in
    if got <> expected then
        fail "%s: expected shape [%s], got [%s]" name
            (String.concat ";" (Array.to_list (Array.map string_of_int expected)))
            (String.concat ";" (Array.to_list (Array.map string_of_int got)))

let check_float name expected actual =
    let d = Float.abs (expected -. actual) in
    if d > 1e-9 then fail "%s: expected %f, got %f (diff %.3e)" name expected actual d

let test_two_atoms_on_x_axis () =
    (* atoms at (1,0,0) and (-1,0,0): centroid is already the origin *)
    let coords_in = [| (1.0, 0.0, 0.0); (-1.0, 0.0, 0.0) |] in
    let coords = M.center_coords coords_in in
    check_shape "coords" [| 3; 2 |] coords;
    check_float "atom0 x" 1.0 (Nd.get coords [| 0; 0 |]);
    check_float "atom1 x" (-1.0) (Nd.get coords [| 0; 1 |]);

    let r = M.radii coords in
    check_shape "r" [| 1; 2 |] r;
    check_float "atom0 r" 1.0 (Nd.get r [| 0; 0 |]);
    check_float "atom1 r" 1.0 (Nd.get r [| 0; 1 |]);

    let angles = M.angles coords r in
    check_shape "angles" [| 2; 2 |] angles;
    (* atom0 at (1,0,0): theta = acos(0/1) = pi/2, phi = atan2(0,1) = 0 *)
    check_float "atom0 theta" (Float.pi /. 2.0) (Nd.get angles [| 0; 0 |]);
    check_float "atom0 phi" 0.0 (Nd.get angles [| 1; 0 |]);
    (* atom1 at (-1,0,0): theta = pi/2, phi = atan2(0,-1) = pi *)
    check_float "atom1 theta" (Float.pi /. 2.0) (Nd.get angles [| 0; 1 |]);
    check_float "atom1 phi" Float.pi (Nd.get angles [| 1; 1 |])

let test_centering_shifts_to_centroid () =
    (* atoms at (0,0,0) and (2,0,0): centroid is (1,0,0), so centered
       positions should be (-1,0,0) and (1,0,0) *)
    let coords_in = [| (0.0, 0.0, 0.0); (2.0, 0.0, 0.0) |] in
    let coords = M.center_coords coords_in in
    check_float "atom0 x after centering" (-1.0) (Nd.get coords [| 0; 0 |]);
    check_float "atom1 x after centering" 1.0 (Nd.get coords [| 0; 1 |])

let test_single_atom_r_zero_no_nan () =
    (* a single atom is its own centroid, so r = 0; must not produce NaN *)
    let coords_in = [| (5.0, 5.0, 5.0) |] in
    let coords = M.center_coords coords_in in
    check_float "single atom x" 0.0 (Nd.get coords [| 0; 0 |]);
    let r = M.radii coords in
    check_float "single atom r" 0.0 (Nd.get r [| 0; 0 |]);
    let angles = M.angles coords r in
    let theta = Nd.get angles [| 0; 0 |] in
    if Float.is_nan theta then fail "single atom theta: got NaN"

let test_empty_coords_raises () =
    match M.center_coords [||] with
    | _ -> fail "empty coords: expected Molecule_error, none was raised"
    | exception M.Molecule_error _ -> ()

let test_vols_known_elements () =
    let v = M.vols [| "fe"; "o"; "rn" |] in
    check_shape "vols" [| 1; 3 |] v;
    let expected_volume (r : float) : float = (4.0 /. 3.0) *. Float.pi *. (r ** 3.0) in
    check_float "fe volume" (expected_volume 1.274) (Nd.get v [| 0; 0 |]);
    check_float "rn volume" (expected_volume 2.24) (Nd.get v [| 0; 2 |])

let test_vols_unknown_element_raises () =
    match M.vols [| "zzzz" |] with
    | _ -> fail "unknown element: expected Molecule_error, none was raised"
    | exception M.Molecule_error _ -> ()

let test_vols_empty_raises () =
    match M.vols [||] with
    | _ -> fail "empty elms: expected Molecule_error, none was raised"
    | exception M.Molecule_error _ -> ()

let test_create_wires_lazy_fields_correctly () =
    let m =
        M.create "test" [| "o"; "h"; "h" |]
            [| (0.0, 0.0, 0.0); (1.0, 0.0, 0.0); (0.0, 1.0, 0.0) |]
    in
    check_shape "create coords" [| 3; 3 |] m.coords;
    let r = Helpers.Cache.force m.r in
    check_shape "create r" [| 1; 3 |] r;
    check_float "create r atom0" 0.4714045208 (Nd.get r [| 0; 0 |]);
    let angles = Helpers.Cache.force m.angles in
    check_shape "create angles" [| 2; 3 |] angles;
    let vols = Helpers.Cache.force m.vols in
    check_shape "create vols" [| 1; 3 |] vols;
    (* forcing again must return the same cached values, not recompute *)
    let r2 = Helpers.Cache.force m.r in
    check_float "create r stable on repeat force" (Nd.get r [| 0; 0 |]) (Nd.get r2 [| 0; 0 |])

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
    test_vols_empty_raises ();
    test_create_wires_lazy_fields_correctly ();
    test_create_length_mismatch_raises ();
    if !failures > 0 then begin
        Printf.printf "%d check(s) failed\n" !failures;
        exit 1
    end else
        print_endline "all checks passed"
