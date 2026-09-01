module AtomicRadiiSqlite3Source : sig
    (** @param ions atom/ion symbols to look up, e.g. ["fe3+"], ["fe+3"], ["fe"]
        @return each input paired with its radius (Angstrom), or [None]
        if nothing matched at any tier of the fallback chain *)
    val lookup : string Seq.t -> (string * float option) Seq.t
end
