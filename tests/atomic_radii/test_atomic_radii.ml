(* Checks for Atomic_radii.Parser / Lexer / Db against known values and
the real atomic_radii.sqlite3. Run from the repo root (db path is
relative). Prints nothing on success; on failure prints one line per
failing check and exits non-zero. *)

module P = Atomic_radii.Parser
module D = Atomic_radii.Db

let failures = ref 0

let fail fmt = Printf.ksprintf (fun s -> incr failures; print_endline ("FAIL " ^ s)) fmt

let check_int name expected actual =
    if expected <> actual then fail "%s: expected %d, got %d" name expected actual

let check_string name expected actual =
    if expected <> actual then fail "%s: expected %S, got %S" name expected actual

let check_float name expected actual =
    let d = Float.abs (expected -. actual) in
    if d > 1e-9 then fail "%s: expected %f, got %f (diff %.3e)" name expected actual d

let check_ion name expected (actual : P.ion) =
    check_string (name ^ " element") expected.P.element actual.P.element;
    check_int (name ^ " charge") expected.P.charge actual.P.charge

let check_some name = function
    | Some x -> x
    | None -> fail "%s: expected Some, got None" name; assert false

let check_none name = function
    | None -> ()
    | Some _ -> fail "%s: expected None, got Some" name

let check_raises name f =
    match f () with
    | _ -> fail "%s: expected an exception, none was raised" name
    | exception _ -> ()

(* ---- Parser ---- *)

let test_parser_orderings () =
    check_ion "fe3+ (elem charge sign)"
        { element = "fe"; charge = 3 } (P.parse_string "fe3+");
    check_ion "au1- (elem charge sign)"
        { element = "au"; charge = -1 } (P.parse_string "au1-");
    check_ion "na+ (bare sign, elem sign)"
        { element = "na"; charge = 1 } (P.parse_string "na+");
    check_ion "cl- (bare sign, elem sign)"
        { element = "cl"; charge = -1 } (P.parse_string "cl-")

let test_parser_malformed () =
    check_raises "empty string" (fun () -> P.parse_string "");
    check_raises "no charge or sign" (fun () -> P.parse_string "fe");
    check_raises "digits with no sign" (fun () -> P.parse_string "fe3");
    check_raises "digits first, no element" (fun () -> P.parse_string "3+");
    check_raises "trailing garbage" (fun () -> P.parse_string "fe3+x")

let test_parser_round_trip () =
    List.iter
        (fun s -> check_string ("round trip " ^ s) s (P.to_string (P.parse_string s)))
        [ "fe3+"; "au1-"; "na1+"; "cl1-"; "fe4+" ];
    check_raises "to_string of charge 0"
        (fun () -> P.to_string { element = "fe"; charge = 0 })

(* ---- Db (against the real atomic_radii.sqlite3) ---- *)

let test_db_ion_radius db =
    check_float "ion_radius fe3+" 49.0
        (check_some "ion_radius fe3+" (D.ion_radius db "fe3+"));
    check_none "ion_radius zz9+ (nonexistent)"
        (D.ion_radius db "zz9+")

let test_db_ion_radius_of db =
    let r =
        check_some "ion_radius_of au1-" (D.ion_radius_of db (P.parse_string "au1-"))
    in
    check_float "ion_radius_of au1-" 220.0 r

let test_db_element_radius db =
    let r, t = check_some "element_radius rn" (D.element_radius db "rn") in
    check_float "element_radius rn value" 2.24 r;
    check_string "element_radius rn type" "vdw" t;
    let r2, t2 = check_some "element_radius fe" (D.element_radius db "fe") in
    check_float "element_radius fe value" 1.274 r2;
    check_string "element_radius fe type" "metallic" t2;
    check_none "element_radius nonexistent" (D.element_radius db "zz")

let test_db_nearest_ion db =
    check_string "nearest_ion fe target 5"
        "fe4+" (check_some "nearest_ion fe target 5" (
            D.nearest_ion db ~element:"fe" ~charge:5)
            )
        ;
    check_string "nearest_ion fe target 3 (exact)"
        "fe3+" (check_some "nearest_ion fe target 3" (
            D.nearest_ion db ~element:"fe" ~charge:3)
            )
        ;
    check_none "nearest_ion nonexistent element"
        (D.nearest_ion db ~element:"zz" ~charge:1)

let test_db_charges_for db =
    let cs = D.charges_for db "fe" in
    if cs <> [ 2; 3; 4; 6 ] then
        fail "charges_for fe: expected [2;3;4;6], got [%s]"
            (String.concat ";" (List.map string_of_int cs));
    if D.charges_for db "zz" <> [] then fail "charges_for zz: expected []"

let test_db_element_radii_batch db =
    let results = D.element_radii_batch db [ "fe"; "o"; "zz"; "fe" ] in
    let elts = List.sort compare (List.map fst results) in
    (* "zz" omitted (not found); "fe" queried twice but that's the
       caller's input, not deduped by this function -- so two "fe" rows. *)
    if elts <> [ "fe"; "fe"; "o" ] then
        fail "element_radii_batch: expected [fe;fe;o], got [%s]" (String.concat ";" elts)

let () =
    test_parser_orderings ();
    test_parser_malformed ();
    test_parser_round_trip ();
    let db = D.open_db () in
    test_db_ion_radius db;
    test_db_ion_radius_of db;
    test_db_element_radius db;
    test_db_nearest_ion db;
    test_db_charges_for db;
    test_db_element_radii_batch db;
    D.close_db db;
    if !failures > 0 then begin
        Printf.printf "%d check(s) failed\n" !failures;
        exit 1
    end else
        print_endline "all checks passed"
