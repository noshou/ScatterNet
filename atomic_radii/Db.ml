(** sqlite3 query primitives over [atomic_radii.sqlite3]. *)
open Sqlite3
open Sqlite3.Data
open Printf

(* path to atomic_radii.sqlite3, relative to repo root *)
let _DB_PATH = "atomic_radii/atomic_radii.sqlite3"

(** opens the dtabase. [~mutex:`NO] is safe here because every caller is
expected to hold its own private connection (see [Atomic_radii]'s
per-domain connection), never a connection shared across domains, so
sqlite3's own internal locking would be pure overhead.
@return an open db handle *)
let open_db () : db = db_open ~mutex:`NO _DB_PATH

(**  closes database.
@param db handle to close
@raise Failure if sqlite3 refuses to close (statements still open) *)
let close_db (db : db) : unit =
    if not (db_close db) then
        failwith "Db.close_db: sqlite3 refused to close (statements still open?)"

(** Run a one-shot prepared statement, expecting at most one row.
@param db handle to query
@param  sql statement text with [?] placeholders
@param  params values bound to the placeholders, in order
@param  extract pulls the caller's value out of the matched row
@return [Some (extract row)], or [None] if no row matched *)
let query_one
    (db : db) (sql : string) (params : t list)
    (extract : stmt -> 'a) : 'a option =
    let res = prepare db sql in
    List.iteri (fun i p -> ignore (bind res (i + 1) p)) params;
    let result =
        match step res with
        | Rc.ROW -> Some (extract res)
        | Rc.DONE -> None
        | rc -> failwith (sprintf "Db.query_one: step failed: %s" (Rc.to_string rc))
    in
    ignore (finalize res);
    result

(** Exact charge-aware lookup. [ionic_radii.radius] is stored in
picometers (Shannon's raw published values, unconverted); this divides
by 100 so every radius-returning function in this module is consistently
in Angstrom, matching [atomic_radii.radius].
@param  db handle to query
@param  ion db-key-format ion string, e.g. ["fe3+"]
@return the ion's radius in Angstrom, or [None] if not on file *)
let ion_radius (db : db) (ion : string) : float option =
    query_one db "SELECT radius FROM ionic_radii WHERE ion = ?"
    [ TEXT ion ]
    (fun res -> to_float_exn (column res 0) /. 100.0)

(** {!ion_radius} from an already-parsed {!Parser.ion}.
@param  db handle to query
@param  ion parsed ion record, reconstructed to db-key format internally
@return the ion's radius in Angstrom, or [None] if not on file *)
let ion_radius_of (db : db) (ion : Parser.ion) : float option =
    ion_radius db (Parser.to_string ion)

(** Bare-element fallback radius.
@param db handle to query
@param element lowercase element symbol, no charge, e.g. ["fe"]
@return [Some (radius, radius_type)] where [radius_type] is
    ["vdw"] | ["metallic"] | ["covalent"], or [None] if not on file *)
let element_radius (db : db) (element : string) : (float * string) option =
    query_one db "SELECT radius, radius_type FROM atomic_radii WHERE element = ?"
    [ TEXT element ]
    (fun res -> (to_float_exn  (column res 0), to_string_exn (column res 1)))

(** Nearest charge state on file for [element] to [charge], via the
[element_charges] cache (e.g. [charge = 5] with only fe2+/3+/4+/6+ on
file returns ["fe4+"]).
@param  db handle to query
@param  element lowercase element symbol
@param  charge desired signed charge
@return the closest available ion string, or [None] if [element] has
        no charge states on file at all *)
let nearest_ion (db : db) ~(element : string) ~(charge : int) : string option =
    let query = 
    "SELECT ion FROM element_charges WHERE element = ? ORDER BY ABS(charge - ?) LIMIT 1"
    in 
    query_one db query
    [ TEXT element; INT (Int64.of_int charge) ]
    (fun res -> to_string_exn (column res 0))

(** Bare-element lookup over a list, one prepared statement reused via
[reset] instead of re-preparing SQL per element.
@param  db handle to query
@param  elements element symbols to look up; duplicates re-query
@return one [(element, (radius, radius_type))] pair per element found;
        elements with no [atomic_radii] row are omitted, not [None]-padded *)
let element_radii_batch 
    (db : db) (elements : string list) : (string * (float * string)) list =
    let query = "SELECT radius, radius_type FROM atomic_radii WHERE element = ?" in
    let res = prepare db query in
    let lookup_one (elt : string) : (string * (float * string)) option =
        ignore (reset res);
        ignore (bind res 1 (TEXT elt));
        match step res with
        | Rc.ROW ->
            let r = to_float_exn (column res 0) in
            let t = to_string_exn (column res 1) in
            Some (elt, (r, t))
        | Rc.DONE -> None
        | rc -> failwith (
            sprintf "Db.element_radii_batch: step failed: %s" (
                Rc.to_string rc
                )
            )
    in
    let results = List.filter_map lookup_one elements in
    ignore (finalize res);
    results

(** Every charge state on file for [element], sorted.
@param  db handle to query
@param  element lowercase element symbol
@return signed charges on file, e.g. [[-2; 2; 3; 4; 6]]; [[]] if none *)
let charges_for (db : db) (element : string) : int list =
    let query = "SELECT charge FROM element_charges WHERE element = ? ORDER BY charge" in
    let res = prepare db query in 
    ignore (bind res 1 (TEXT element));
    let rec loop acc =
        match step res with
        | Rc.ROW -> loop (to_int_exn (column res 0) :: acc)
        | Rc.DONE -> List.rev acc
        | rc -> failwith (sprintf "Db.charges_for: step failed: %s" (Rc.to_string rc))
    in
    let result = loop [] in
    ignore (finalize res);
    result
