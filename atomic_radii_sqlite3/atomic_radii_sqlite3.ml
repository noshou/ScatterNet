open Option
open Db

(** Public surface of this library. Each domain opens and keeps its own
sqlite3 connection, rather than sharing one handle across domains: a
single connection is not safe to use from multiple domains at once
unless the linked libsqlite3 was compiled with mutex support, and
nothing in OCaml can verify that at runtime, so a shared handle isn't
trustworthy. A private connection per domain sidesteps the question
entirely. *)
let db_key : Sqlite3.db Domain.DLS.key = Domain.DLS.new_key open_db

(** Fallback chain for one ion/element string: (1) exact charge-aware
    match; (2) nearest available charge state for that element; (3)
    bare-element radius, dropping [radius_type]; (4) [None]. A string that
    doesn't parse as a charge-bearing ion at all (e.g. bare ["fe"]) skips
    straight to (3). *)
let resolve_one (ion : string) : float option =
    let db = Domain.DLS.get db_key in
    match Parser.parse_string ion with
    | parsed -> (
        match ion_radius_of db parsed with
        | Some r -> Some r
        | None -> (
            match nearest_ion db ~element:parsed.element ~charge:parsed.charge with
            | Some nearest -> (
                match ion_radius db nearest with
                | Some r -> Some r
                | None -> map fst (element_radius db parsed.element))
            | None -> map fst (element_radius db parsed.element)))
    | exception Parser.Parse_error _ -> map fst (element_radius db ion)

module AtomicRadiiSqlite3Source = struct
    (** Batch ion/element lookup, deduped within one call: each distinct input string is 
        resolved via {!resolve_one} at most once, cached for repeats.
        @param  ions atom/ion symbols to look up, e.g. ["fe3+"], ["fe+3"], ["fe"]
        @return each input paired with radius (Angstrom), or [None] if nothing matched. *)
    let lookup (ions : string Seq.t) : (string * float option) Seq.t =
        let cache : (string, float option) Hashtbl.t = Hashtbl.create 16 in
        Seq.map (fun ion -> 
            match Hashtbl.find_opt cache ion with
                | Some cached -> (ion, cached)
                | None ->
                    let result = resolve_one ion in
                    Hashtbl.add cache ion result;
                    (ion, result)
        ) ions
end
