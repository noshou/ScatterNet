open Atomic_volume
open Helpers

(** Raised when coordinates, radii, or their shapes are invalid. *)
exception Molecule_error of string

module Nd = Owl_dense_ndarray_d

(** Volume lookup backed by the real {!Atomic_radii_sqlite3.AtomicRadiiSqlite3Source}. *)
module AV = AtomicVolume (Atomic_radii_sqlite3.AtomicRadiiSqlite3Source)

type t = {
    angles : Nd.arr Cache.t; (** (2, n): theta, phi in radians. *)
    r      : Nd.arr Cache.t; (** (1, n): radial distance from the centroid. *)
    vols   : Nd.arr Cache.t; (** (1, n): per-atom excluded volume. *)
    elms   : string array;   (** element symbol per atom. *)
    name   : string;         (** molecule name. *)
}

(** Centers coordinates at their centroid.
    @param coords per-atom (x, y, z)
    @return (3, n), centered
    @raise Molecule_error if [coords] is empty *)
let center_coords (coords : (float * float * float) array) : Nd.arr =
    
    let len = float_of_int (Array.length coords) in
    
    if  len = 0. then
        raise (Molecule_error "Empty coordinates")
    
    else
        
        (*  the idomatic way is below: 
            
            let (sum_x, sum_y, sum_z) =
                Array.fold_left (fun (acc_x, acc_y, acc_z) (x, y, z) ->
                    (acc_x +. x, acc_y +. y, acc_z +. z)
                ) (0.0, 0.0, 0.0) coords
            in
            let (mu_x, mu_y, mu_z) = (sum_x /. len, sum_y /. len, sum_z /. len) in 
            
            but this is is very inefficient for large molecules, since new tuples are
            created for every iteration. Below is the imperative way which is better. *)
        
        begin
            (* initialize mutable variables *)
            let sum_x = ref 0.0 in
            let sum_y = ref 0.0 in
            let sum_z = ref 0.0 in

            (* calculate sum directly over array *)
            Array.iter (fun (x, y, z) ->
            sum_x := !sum_x +. x;
            sum_y := !sum_y +. y;
            sum_z := !sum_z +. z
                        ) coords;
            
            (* calculate mean *)
            let (mu_x, mu_y, mu_z) = (!sum_x /. len, !sum_y /. len, !sum_z /. len) in

            (* center coordinates into a fresh ndarray *)
            let centered = Nd.zeros [| 3; Array.length coords |] in
            Array.iteri (fun i (x, y, z) ->
                Nd.set centered [| 0; i |] (x -. mu_x);
                Nd.set centered [| 1; i |] (y -. mu_y);
                Nd.set centered [| 2; i |] (z -. mu_z)
            ) coords;
            centered
        end

(** Row [i] of a (3, n) or (2, n) array, e.g. [row 0 coords] is x. *)
let row (i : int) (a : Nd.arr) : Nd.arr = Nd.get_slice [ [ i ]; [] ] a

(** Radial distance of each atom from the centroid, given its already-
    sliced x/y/z rows. Shared with {!angles_of_xyz} via {!create}, so
    [coords] only ever gets sliced into x/y/z once per molecule instead
    of once per field. *)
let radii_of_xyz (x : Nd.arr) (y : Nd.arr) (z : Nd.arr) : Nd.arr =
    Nd.sqrt (Nd.add (Nd.add (Nd.mul x x) (Nd.mul y y)) (Nd.mul z z))

(** Polar (theta) and azimuthal (phi) angle of each atom, given already-
    sliced x/y/z rows and [r] from {!radii_of_xyz}. Physics convention:
    theta = arccos(z/r) in \[0, pi\], phi = atan2(y, x).
    r=0 only occurs for single-atom molecules (atom is its own centroid).
    j_l(0)=0 for l>0 so the angle is irrelevant; use r_safe to avoid NaN. *)
let angles_of_xyz (x : Nd.arr) (y : Nd.arr) (z : Nd.arr) (r : Nd.arr) : Nd.arr =
    let r_safe = Nd.map (fun _r -> if _r > 0.0 then _r else 1.0) r in
    let theta = Nd.acos (Nd.clip_by_value ~amin:(-1.) ~amax:1. (Nd.div z r_safe)) in
    let phi = Nd.atan2 y x in
    Nd.concatenate ~axis:0 [| theta; phi |]

(** Per-atom excluded volume, from atomic/ionic radii looked up by
    element symbol.
    @param elms element symbol per atom
    @return (1, n)
    @raise Molecule_error if [elms] is empty, or any element has no
    radius on file at all *)
let compute_vols (elms : string array) : Nd.arr =
    let n = Array.length elms in
    if n = 0 then
        raise (Molecule_error "Empty elements")
    else
        begin
            let results = Array.of_seq (AV.get_vols (Array.to_seq elms)) in
            let out = Nd.zeros [| 1; n |] in
            Array.iteri
                (fun i (elem, result) ->
                    match result with
                    | Some (_, v) -> Nd.set out [| 0; i |] v
                    | None ->
                        raise (Molecule_error (
                            Printf.sprintf "no volume data for element %S" elem))
                        )
                results;
            out
        end

(** Build a molecule from raw per-atom Cartesian coordinates, element
    symbols, and a name. Coordinates are centered at their centroid;
    [r], [angles], and [vols] are computed lazily on first use.
    @param name molecule name
    @param elms element symbol per atom
    @param coords_in per-atom (x, y, z), any coordinate system
    @return a new [t]
    @raise Molecule_error if [coords_in]/[elms] are empty or mismatched *)
let create 
    (name : string) (elms : string array) 
    (coords_in : (float * float * float) array) : t =
    if Array.length coords_in <> Array.length elms then
        raise (Molecule_error "coords and elms must be the same length")
    else
        let coords = center_coords coords_in in
        let x = row 0 coords and y = row 1 coords and z = row 2 coords in
        let r_cache = Cache.make (fun () -> radii_of_xyz x y z) in
        let a_cache = Cache.make (fun () -> angles_of_xyz x y z (Cache.force r_cache)) in
        let vols_cache = Cache.make (fun () -> compute_vols elms) in
        { angles = a_cache; r = r_cache; vols = vols_cache; elms; name }

(** Flatten a (1, n) row into a genuine 1-D (n,) array. *)
let flatten_row (a : Nd.arr) : Nd.arr = Nd.reshape a [| (Nd.shape a).(1) |]

(** Polar angle per atom, flat.
    @param m molecule
    @return (n,), radians *)
let theta (m : t) : Nd.arr = flatten_row (Nd.get_slice [[0]; []] (Cache.force m.angles))

(** Azimuthal angle per atom, flat.
    @param m molecule
    @return (n,), radians *)
let phi (m : t) : Nd.arr = flatten_row (Nd.get_slice [[1]; []] (Cache.force m.angles))

(** Radial distance from centroid per atom, flat.
    @param m molecule
    @return (n,), Angstrom *)
let r (m : t) : Nd.arr = flatten_row (Cache.force m.r)

(** Per-atom excluded volume, flat.
    @param m molecule
    @return (n,), Angstrom^3
    @raise Molecule_error if any element in [m] has no radius on file *)
let vols (m : t) : Nd.arr = flatten_row (Cache.force m.vols)

(** Element symbol per atom.
    @param m molecule
    @return elms, as passed to {!create} *)
let elms (m : t) : string array = m.elms

(** Molecule name.
    @param m molecule
    @return name, as passed to {!create} *)
let name (m : t) : string = m.name