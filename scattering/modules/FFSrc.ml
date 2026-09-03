(** Module signature for form factor sources.

    A source builds a per-batch container [t] from a set of ions at one
    energy over a shared q grid, then answers per-ion, per-q-point
    queries against it. *)
module type FormFactorSource = sig

    (** Per-batch form factor container. *)
    type t

    (** Build a container for a batch of ions at one energy, over a q grid.
        @param energy incident energy, eV
        @param ions ion/element symbols to look up
        @param qvals q grid, Angstrom^-1
        @return a new [t] *)
    val create :
        float -> string list -> Owl_dense_ndarray_d.arr -> t

    (** Construction-time log for [t], one line per ion that wasn't a
        full (f0 + anomalous) lookup. *)
    val log : t -> string array

    (** Form factors for the given ions at the given q-points. Ions with
        no data are dropped; check {!log} for details.
        @return one [(ion, form factors)] pair per ion found, each array
        positionally aligned with [qvals] *)
    val lookup :
        t -> string list -> Owl_dense_ndarray_d.arr ->
        (string * Owl_dense_ndarray_z.arr) list
end
