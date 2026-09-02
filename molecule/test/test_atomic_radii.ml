(* Checks for Atomic_radii.AtomicRadiiSqlite3Source -- the library's only public
surface (Db/Parser/Toks/Lexer are sealed, see atomic_radii.mli). Run from
the repo root (db path is relative). Prints nothing on success; on
failure prints one line per failing check and exits non-zero. *)

module S = Atomic_radii_sqlite3.AtomicRadiiSqlite3Source

let failures = ref 0

let fail fmt = Printf.ksprintf (fun s -> incr failures; print_endline ("FAIL " ^ s)) fmt

let check_float name expected actual =
    let d = Float.abs (expected -. actual) in
    if d > 1e-9 then fail "%s: expected %f, got %f (diff %.3e)" name expected actual d

let lookup_one (ion : string) : float option =
    match S.lookup [| ion |] with
    | [| (got_ion, r) |] ->
        if got_ion <> ion then fail "lookup_one %s: echoed ion was %S" ion got_ion;
        r
    | results ->
        fail "lookup_one %s: expected exactly 1 result, got %d"
            ion (Array.length results);
        None

let check_some name expected = function
    | Some r -> check_float name expected r
    | None -> fail "%s: expected Some %f, got None" name expected

let check_none name = function
    | None -> ()
    | Some r -> fail "%s: expected None, got Some %f" name r

(* ---- exact ion match, both token orderings ---- *)

let test_exact_ion_both_orderings () =
    check_some "fe3+" 0.49 (lookup_one "fe3+");
    check_some "fe+3" 0.49 (lookup_one "fe+3");
    check_some "au1-" 2.2 (lookup_one "au1-")

(* ---- bare element falls through to atomic_radii ---- *)

let test_bare_element_fallback () =
    check_some "fe" 1.274 (lookup_one "fe");
    check_some "rn" 2.24 (lookup_one "rn")

(* ---- nearest-charge fallback: fe has 2+/3+/4+/6+ on file, not 5+ ---- *)

let test_nearest_charge_fallback () =
    let r5 = lookup_one "fe5+" in
    let r4 = lookup_one "fe4+" in
    match r5, r4 with
    | Some v5, Some v4 -> check_float "fe5+ (nearest -> fe4+)" v4 v5
    | _ ->  fail "fe5+/fe4+: expected both to resolve, got %b/%b" 
            (Option.is_some r5) (Option.is_some r4)

(* ---- total miss ---- *)

let test_total_miss () = check_none "zzzz9+ (garbage)" (lookup_one "zzzz9+")

(* ---- batch: order/count preserved, repeats deduped to the same result ---- *)

let test_batch_order_and_dedup () =
    let input = [| "fe3+"; "zzzz9+"; "fe3+"; "rn" |] in
    let results = S.lookup input in
    let elems = Array.map fst results in
    if elems <> input then
        fail "batch: expected order [%s], got [%s]"
            (String.concat ";" (Array.to_list input))
            (String.concat ";" (Array.to_list elems));
    match results with
    | [| (_, Some r1); (_, None); (_, Some r2); (_, Some _) |] ->
            check_float "batch: repeated fe3+" r1 r2
    | _ -> fail "batch: unexpected result shape"

let () =
    test_exact_ion_both_orderings ();
    test_bare_element_fallback ();
    test_nearest_charge_fallback ();
    test_total_miss ();
    test_batch_order_and_dedup ();
    if !failures > 0 then begin
        Printf.printf "%d check(s) failed\n" !failures;
        exit 1
    end else
        print_endline "all checks passed"
