(* Inline tests (ppx_inline_test) for AtomicRadiiSqlite3Source: fallback
chain, batch order/dedup. This module is internal to the molecule library
(see molecule.mli), so these tests live inside the library rather than as
an external test executable. *)

let check_float expected actual = Float.abs (expected -. actual) < 1e-9

let lookup_one (ion : string) : float option =
    match AtomicRadiiSqlite3.AtomicRadiiSqlite3Source.lookup [| ion |] with
    | [| (_, r) |] -> r
    | _ -> assert false

let%test "exact ion match, both token orderings" =
    check_float 0.49 (Option.get (lookup_one "fe3+"))
    && check_float 0.49 (Option.get (lookup_one "fe+3"))
    && check_float 2.2 (Option.get (lookup_one "au1-"))

let%test "bare element falls through to atomic_radii" =
    check_float 1.274 (Option.get (lookup_one "fe"))
    && check_float 2.24 (Option.get (lookup_one "rn"))

(* fe has 2+/3+/4+/6+ on file, not 5+ *)
let%test "nearest-charge fallback" =
    match lookup_one "fe5+", lookup_one "fe4+" with
    | Some v5, Some v4 -> check_float v4 v5
    | _ -> false

let%test "total miss" = lookup_one "zzzz9+" = None

let%test "batch: order/count preserved, repeats deduped to the same result" =
    let input = [| "fe3+"; "zzzz9+"; "fe3+"; "rn" |] in
    let results = AtomicRadiiSqlite3.AtomicRadiiSqlite3Source.lookup input in
    Array.map fst results = input
    &&
    match results with
    | [| (_, Some r1); (_, None); (_, Some r2); (_, Some _) |] -> check_float r1 r2
    | _ -> false
