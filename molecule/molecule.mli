(** Raised when coordinates, radii, or their shapes are invalid. *)
exception Molecule_error of string

(** A molecule: per-atom Cartesian coordinates centered at their centroid,
    plus derived quantities computed lazily on first access. *)
type t

(** Build a molecule from raw per-atom Cartesian coordinates, element
    symbols, and a name. Coordinates are centered at their centroid;
    {!r}, {!theta}/{!phi}, and {!vols} are computed lazily on first use.
    @param name molecule name
    @param elms element symbol per atom
    @param coords_in per-atom (x, y, z), any coordinate system
    @return a new [t]
    @raise Molecule_error if [coords_in]/[elms] are empty or mismatched *)
val create : string -> string array -> (float * float * float) array -> t

(** Polar angle per atom, flat.
    @param m molecule
    @return (n,), radians *)
val theta : t -> Owl_dense_ndarray_d.arr

(** Azimuthal angle per atom, flat.
    @param m molecule
    @return (n,), radians *)
val phi : t -> Owl_dense_ndarray_d.arr

(** Radial distance from centroid per atom, flat.
    @param m molecule
    @return (n,), Angstrom *)
val r : t -> Owl_dense_ndarray_d.arr

(** Per-atom excluded volume, flat.
    @param m molecule
    @return (n,), Angstrom^3
    @raise Molecule_error if any element in [m] has no radius on file *)
val vols : t -> Owl_dense_ndarray_d.arr

(** Element symbol per atom.
    @param m molecule
    @return elms, as passed to {!create} *)
val elms : t -> string array

(** Molecule name.
    @param m molecule
    @return name, as passed to {!create} *)
val name : t -> string
