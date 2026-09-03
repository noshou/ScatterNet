(** Defines a source for atomic radii with simple batch lookup. *)
module type RadiiSource = sig

    (** returns the atom/ion, paired with radius if found,
        or the atom/ion alone if not found
        @param ions atom/ion symbols to look up
        @return each input paired with its radius (Angstrom), or [None] if not found *)
    val lookup : string array -> (string * float option) array
end
