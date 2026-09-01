(** Defines a source for atomic radii with simple batch lookup. *)
module type RadiiSource = sig

    (** returns the atom/ion, paired with radius if found,
        or the atom/ion alone if not found
        @param ions atom/ion symbols to look up
        @return each input paired with its radius (Angstrom), or [None] if not found *)
    val lookup : string Seq.t -> (string * float option) Seq.t
end

(** Functor that takes a source of atomic radii, a batch of ions, and returns volumes *)
module AtomicVolume (RadSrc : RadiiSource) = struct

    (** returns atom/ion paired with (radius, volume) if found,
        or atom/ion alone if not found.
        @param ions atom/ion symbols to look up
        @return each input paired with [(radius, volume)] (Angstrom, Angstrom^3),
        or [None] if [RadSrc.lookup] found nothing for it *)
    let get_vols (ions : string Seq.t) : (string * (float * float) option) Seq.t =
        let rads = RadSrc.lookup ions in
        Seq.map (fun rad ->
            match rad with
                | (elem, Some r) -> (elem, Some (r, (4. /. 3.) *. Float.pi *. r ** 3.))
                | (elem, None)   -> (elem, None)
        ) rads
end
