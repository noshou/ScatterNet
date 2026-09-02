(* Inline tests (ppx_inline_test) for Cache. *)

let%test "force memoizes: thunk runs exactly once across repeat forces" =
    let calls = ref 0 in
    let c = Cache.make (fun () -> incr calls; 42) in
    let a = Cache.force c in
    let b = Cache.force c in
    a = 42 && b = 42 && !calls = 1

let%test "make does not run the thunk until force is called" =
    let calls = ref 0 in
    let _c = Cache.make (fun () -> incr calls; 1) in
    !calls = 0

let%test "concurrent force from multiple domains still runs the thunk once" =
    let calls = ref 0 in
    let c = Cache.make (fun () -> incr calls; 7) in
    let domains = Array.init 8 (fun _ -> Domain.spawn (fun () -> Cache.force c)) in
    let results = Array.map Domain.join domains in
    Array.for_all (( = ) 7) results && !calls = 1

let%test "distinct caches are independent" =
    let c1 = Cache.make (fun () -> 1) in
    let c2 = Cache.make (fun () -> 2) in
    Cache.force c1 = 1 && Cache.force c2 = 2
