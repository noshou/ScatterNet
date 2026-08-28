(* Checks for Scattering.AtomicVolume against a dummy RadiiSource (no db
involved). Prints nothing on success; on failure prints one line per
failing check and exits non-zero. *)

module A = Scattering.AtomicVolume

(* dummy source: fe -> 2.0, h -> 1.0, everything else -> not found *)
module Dummy : A.RadiiSource = struct
    let lookup (ions : string Seq.t) : (string * float option) Seq.t =
        Seq.map
            (fun ion ->
                match ion with
                | "fe" -> (ion, Some 2.0)
                | "h" -> (ion, Some 1.0)
                | _ -> (ion, None))
            ions
end

module V = A.AtomicVolume (Dummy)

let failures = ref 0

let fail fmt = Printf.ksprintf (fun s -> incr failures; print_endline ("FAIL " ^ s)) fmt

let sphere_volume (r : float) : float = (4.0 /. 3.0) *. Float.pi *. (r ** 3.0)

let test_found_element_gets_radius_and_volume () =
    let results = List.of_seq (V.get_vols (List.to_seq [ "fe" ])) in
    match results with
    | [ ("fe", Some (r, v)) ] ->
        if r <> 2.0 then fail "fe: expected radius 2.0, got %f" r;
        let expected_v = sphere_volume 2.0 in
        if Float.abs (v -. expected_v) > 1e-9 then
            fail "fe: expected volume %f, got %f" expected_v v
    | _ -> fail "fe: unexpected result shape"

let test_missing_element_gets_none () =
    let results = List.of_seq (V.get_vols (List.to_seq [ "zz" ])) in
    match results with
    | [ ("zz", None) ] -> ()
    | _ -> fail "zz: expected (zz, None)"

let test_batch_preserves_order_and_count () =
    let input = [ "fe"; "zz"; "h" ] in
    let results = List.of_seq (V.get_vols (List.to_seq input)) in
    let elems = List.map fst results in
    if elems <> input then
        fail "batch: expected order [%s], got [%s]"
            (String.concat ";" input) (String.concat ";" elems);
    if List.length results <> 3 then
        fail "batch: expected 3 results, got %d" (List.length results)

let test_empty_batch () =
    let results = List.of_seq (V.get_vols Seq.empty) in
    if results <> [] then fail "empty batch: expected [], got %d results" (List.length results)

let () =
    test_found_element_gets_radius_and_volume ();
    test_missing_element_gets_none ();
    test_batch_preserves_order_and_count ();
    test_empty_batch ();
    if !failures > 0 then begin
        Printf.printf "%d check(s) failed\n" !failures;
        exit 1
    end else
        print_endline "all checks passed"
