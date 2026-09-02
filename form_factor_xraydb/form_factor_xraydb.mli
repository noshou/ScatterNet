(** Raised when the pyml boundary fails (missing python deps, or input
    FormFact_py.py rejected) or when a {!FormFactorSourceXrayDB.lookup}
    query references a q-point outside the container's grid. *)
exception Form_factor_xraydb_error of string

(** Form factors for a batch of ions over a shared q grid, plus the
    construction-time log (one "TIER ion" line per non-full ion). *)
type ff

(** Form factors for a batch of ions, one entry per unique ion.
    @param ions ion/element symbols to look up
    @param energy incident energy, eV
    @param qvals q grid, Angstrom^-1
    @return a new [ff]
    @raise Form_factor_xraydb_error if python rejects the input *)
val compute_form_factors :
    string array -> float -> Owl_dense_ndarray_d.arr -> ff

(** Concrete {!compute_form_factors}-backed form factor source: builds a container from a 
    batch of ions/energy/qvals, then answers per-ion, per-q-point queries against it. *)
module FormFactorSourceXrayDB : sig

    type t

    (** Build a form factor container for a batch of ions at one energy, over a q grid.
        @param energy incident energy, eV
        @param ions ion/element symbols to look up
        @param qvals q grid, Angstrom^-1
        @return a new [t]
        @raise Form_factor_xraydb_error if the underlying python call fails *)
    val create : float -> string list -> Owl_dense_ndarray_d.arr -> t

    (** The construction-time log for [t]
        @param t the built form factor container
        @return the log recorded when [t] was built *)
    val log : t -> string array

    (** Form factors for the given ions at the given q-points. Invalid elements and
        dummy sites are dropped. Check log for more details. 

        @param t the built form factor container
        @param ions the ions to query
        @param qvals the q-points to query
        @return one [(ion, form factors)] pair per ion actually found,
        each form factor array positionally aligned with [qvals]
        @raise Form_factor_xraydb_error if a queried q isn't one of [t]'s
        q-points *)
    val lookup :
        t -> string list -> Owl_dense_ndarray_d.arr ->
        (string * Owl_dense_ndarray_z.arr) list
end
