(* sqlite3 query primitives over atomic_radii.db. *)

type db = Sqlite3.db

let open_db (path : string) : db = Sqlite3.db_open path

let close_db (db : db) : unit =
    if not (Sqlite3.db_close db) then
        failwith "Db.close_db: sqlite3 refused to close (statements still open?)"

(* run a prepared statement with the given bound params, expecting at most
one row back; [extract] pulls the caller's value out of that row. *)
let query_one
    (db : db) (sql : string) (params : Sqlite3.Data.t list)
    (extract : Sqlite3.stmt -> 'a) : 'a option =
    let stmt = Sqlite3.prepare db sql in
    List.iteri (fun i p -> ignore (Sqlite3.bind stmt (i + 1) p)) params;
    let result = match Sqlite3.step stmt with
        | Sqlite3.Rc.ROW -> Some (extract stmt)
        | Sqlite3.Rc.DONE -> None
        | rc ->
            failwith (Printf.sprintf "Db.query_one: sqlite3 step failed: %s"
                (Sqlite3.Rc.to_string rc))
    in
    ignore (Sqlite3.finalize stmt);
    result

(* exact charge-aware lookup, ion in db key format, e.g. "fe3+". *)
let ion_radius (db : db) (ion : string) : float option =
    query_one db
        "SELECT radius FROM ionic_radii WHERE ion = ?"
        [ Sqlite3.Data.TEXT ion ]
        (fun stmt -> Sqlite3.Data.to_float_exn (Sqlite3.column stmt 0))

(* same, but takes an already-parsed ion record. *)
let ion_radius_of (db : db) (ion : Parser.ion) : float option =
    ion_radius db (Parser.to_string ion)

(* bare-element fallback radius, plus which physical quantity it is
   ("vdw" | "metallic" | "covalent") so the caller doesn't blend types. *)
let element_radius (db : db) (element : string) : (float * string) option =
    query_one db
        "SELECT radius, radius_type FROM atomic_radii WHERE element = ?"
        [ Sqlite3.Data.TEXT element ]
        (fun stmt ->
            ( Sqlite3.Data.to_float_exn (Sqlite3.column stmt 0),
              Sqlite3.Data.to_string_exn (Sqlite3.column stmt 1) ))

(* nearest charge state on file for [element] to the requested [charge],
   via the element_charges cache -- e.g. asking for fe5+ when only fe2+/
   fe3+/fe4+/fe6+ exist returns "fe4+". Returns the db-format ion string,
   not a radius, so the caller can re-look-up via ion_radius (or inspect
   which charge state it actually got). *)
let nearest_ion (db : db) ~(element : string) ~(charge : int) : string option =
    query_one db
        "SELECT ion FROM element_charges WHERE element = ? \
         ORDER BY ABS(charge - ?) LIMIT 1"
        [ Sqlite3.Data.TEXT element; Sqlite3.Data.INT (Int64.of_int charge) ]
        (fun stmt -> Sqlite3.Data.to_string_exn (Sqlite3.column stmt 0))

(* batch bare-element lookup: one statement, prepared once, reused via
   [reset] for every element rather than re-preparing the SQL text per
   query -- the intended entry point for Xyz.scan's output, where a file
   may reference thousands of atoms but only a handful of distinct
   elements. Returns only the elements that were actually found. *)
let element_radii_batch (db : db) (elements : string list) : (string * (float * string)) list =
    let stmt = Sqlite3.prepare db "SELECT radius, radius_type FROM atomic_radii WHERE element = ?" in
    let lookup_one (elt : string) : (string * (float * string)) option =
        ignore (Sqlite3.reset stmt);
        ignore (Sqlite3.bind stmt 1 (Sqlite3.Data.TEXT elt));
        match Sqlite3.step stmt with
        | Sqlite3.Rc.ROW ->
            let r = Sqlite3.Data.to_float_exn (Sqlite3.column stmt 0) in
            let t = Sqlite3.Data.to_string_exn (Sqlite3.column stmt 1) in
            Some (elt, (r, t))
        | Sqlite3.Rc.DONE -> None
        | rc ->
            failwith (Printf.sprintf "Db.element_radii_batch: step failed: %s"
                (Sqlite3.Rc.to_string rc))
    in
    let results = List.filter_map lookup_one elements in
    ignore (Sqlite3.finalize stmt);
    results

(* single-pass xyz radius lookup: streams atom lines out of [lexbuf] via
   Xyz.fold and queries atomic_radii the moment an element is first seen,
   caching the result for every repeat. Nothing about the file is ever
   materialized as a per-atom list -- the cache holds one entry per
   *distinct* element, not one per atom, and that's the only extra memory
   this uses beyond the lexer's own buffer. *)
let radii_of_xyz (db : db) (lexbuf : Lexing.lexbuf) : (string * (float * string)) list =
    let stmt = Sqlite3.prepare db "SELECT radius, radius_type FROM atomic_radii WHERE element = ?" in
    let cache : (string, (float * string) option) Hashtbl.t = Hashtbl.create 16 in
    let resolve (elt : string) : unit =
        if not (Hashtbl.mem cache elt) then begin
            ignore (Sqlite3.reset stmt);
            ignore (Sqlite3.bind stmt 1 (Sqlite3.Data.TEXT elt));
            let result =
                match Sqlite3.step stmt with
                | Sqlite3.Rc.ROW ->
                    Some ( Sqlite3.Data.to_float_exn (Sqlite3.column stmt 0),
                           Sqlite3.Data.to_string_exn (Sqlite3.column stmt 1) )
                | Sqlite3.Rc.DONE -> None
                | rc ->
                    failwith (Printf.sprintf "Db.radii_of_xyz: step failed: %s"
                        (Sqlite3.Rc.to_string rc))
            in
            Hashtbl.add cache elt result
        end
    in
    Xyz.fold (fun () elt -> resolve elt) () lexbuf;
    ignore (Sqlite3.finalize stmt);
    Hashtbl.fold
        (fun elt result acc -> match result with Some r -> (elt, r) :: acc | None -> acc)
        cache []

(* every charge state on file for [element], sorted, e.g. "fe" ->
   [-2; 2; 3; 4; 6] if that's what's tabulated. Useful for callers that want
   to pick a fallback charge themselves rather than trusting nearest_ion's
   distance metric. *)
let charges_for (db : db) (element : string) : int list =
    let stmt =
        Sqlite3.prepare db
            "SELECT charge FROM element_charges WHERE element = ? ORDER BY charge"
    in
    ignore (Sqlite3.bind stmt 1 (Sqlite3.Data.TEXT element));
    let rec loop acc =
        match Sqlite3.step stmt with
        | Sqlite3.Rc.ROW ->
            let c = Sqlite3.Data.to_int_exn (Sqlite3.column stmt 0) in
            loop (c :: acc)
        | Sqlite3.Rc.DONE -> List.rev acc
        | rc ->
            failwith (Printf.sprintf "Db.charges_for: sqlite3 step failed: %s"
                (Sqlite3.Rc.to_string rc))
    in
    let result = loop [] in
    ignore (Sqlite3.finalize stmt);
    result
