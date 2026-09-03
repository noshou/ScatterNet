open Helpers

exception Molecule_error of string

(** Atomic/ionic radii source, backed by the real
    {!AtomicRadiiSqlite3.AtomicRadiiSqlite3Source}. Swap this binding
    for a different {!RadSrc.RadiiSource} to change where radii come from. *)
module Radii : RadSrc.RadiiSource = AtomicRadiiSqlite3.AtomicRadiiSqlite3Source

(** Excluded volume of a hard sphere of the given radius
    (Angstrom -> Angstrom^3). *)
let sphere_volume (r : float) : float = (4. /. 3.) *. Float.pi *. r ** 3.

type t = {
    coords : Owl_dense_ndarray_d.arr;         (** (3, n): centered x, y, z in Angstrom *)
    angles : Owl_dense_ndarray_d.arr Cache.t; (** (2, n): theta, phi in radians *)
    r      : Owl_dense_ndarray_d.arr Cache.t; (** (1, n): radial distance from centroid *)
    radii  : Owl_dense_ndarray_d.arr Cache.t; (** (1, n): per-atom atomic/ionic radius *)
    vols   : Owl_dense_ndarray_d.arr Cache.t; (** (1, n): per-atom excluded volume *)
    elms   : string array;   (** element symbol per atom. *)
    name   : string;         (** molecule name. *)
}

(** Centers coordinates at their centroid.
    @param coords per-atom (x, y, z)
    @return (3, n), centered
    @raise Molecule_error if [coords] is empty *)
let center_coords (coords : (float * float * float) array) : Owl_dense_ndarray_d.arr =

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
            let centered = Owl_dense_ndarray_d.zeros [| 3; Array.length coords |] in
            Array.iteri (fun i (x, y, z) ->
                Owl_dense_ndarray_d.set centered [| 0; i |] (x -. mu_x);
                Owl_dense_ndarray_d.set centered [| 1; i |] (y -. mu_y);
                Owl_dense_ndarray_d.set centered [| 2; i |] (z -. mu_z)
            ) coords;
            centered
        end

(** Row [i] of a (3, n) or (2, n) array, e.g. [row 0 coords] is x. *)
let row (i : int) (a : Owl_dense_ndarray_d.arr) :
    Owl_dense_ndarray_d.arr = Owl_dense_ndarray_d.get_slice [ [ i ]; [] ] a

(** Radial distance of each atom from the centroid, given its already-
    sliced x/y/z rows. Shared with {!angles_of_xyz} via {!create}, so
    [coords] only ever gets sliced into x/y/z once per molecule instead
    of once per field. *)
let radii_of_xyz (x : Owl_dense_ndarray_d.arr) (y : Owl_dense_ndarray_d.arr)
    (z : Owl_dense_ndarray_d.arr) : Owl_dense_ndarray_d.arr =

    Owl_dense_ndarray_d.sqrt (
        Owl_dense_ndarray_d.add (
            Owl_dense_ndarray_d.add
                (Owl_dense_ndarray_d.mul x x) (Owl_dense_ndarray_d.mul y y)
        ) (Owl_dense_ndarray_d.mul z z)
    )

(** Polar (theta) and azimuthal (phi) angle of each atom, given already-
    sliced x/y/z rows and [r] from {!radii_of_xyz}. Physics convention:
    theta = arccos(z/r) in \[0, pi\], phi = atan2(y, x).
    r=0 only occurs for single-atom molecules (atom is its own centroid).
    j_l(0)=0 for l>0 so the angle is irrelevant; use r_safe to avoid NaN. *)
let angles_of_xyz
    (x : Owl_dense_ndarray_d.arr) (y : Owl_dense_ndarray_d.arr)
    (z : Owl_dense_ndarray_d.arr) (r : Owl_dense_ndarray_d.arr)
    : Owl_dense_ndarray_d.arr =

    let r_safe = Owl_dense_ndarray_d.map (fun _r -> if _r > 0.0 then _r else 1.0) r in
    let theta = Owl_dense_ndarray_d.acos (
        Owl_dense_ndarray_d.clip_by_value ~amin:(-1.) ~amax:1.
        (Owl_dense_ndarray_d.div z r_safe)
    ) in
    let phi = Owl_dense_ndarray_d.atan2 y x in
    Owl_dense_ndarray_d.concatenate ~axis:0 [| theta; phi |]

(** Per-atom atomic/ionic radius, looked up by element symbol.
    @param elms element symbol per atom
    @return (1, n), Angstrom
    @raise Molecule_error if [elms] is empty, or any element has no
    radius on file at all *)
let compute_radii (elms : string array) : Owl_dense_ndarray_d.arr =
    let n = Array.length elms in
    if n = 0 then
        raise (Molecule_error "Empty elements")
    else
        begin
            let radii = Radii.lookup elms in
            let out = Owl_dense_ndarray_d.zeros [| 1; n |] in
            Array.iteri
                (fun i (elem, radius) ->
                    match radius with
                    | Some r -> Owl_dense_ndarray_d.set out [| 0; i |] r
                    | None ->
                        raise (Molecule_error (
                            Printf.sprintf "no radius data for element %S" elem))
                        )
                radii;
            out
        end

(** Per-atom excluded volume, as the volume of a hard sphere of each
    atom's radius.
    @param radii per-atom radius, (1, n), from {!compute_radii}
    @return (1, n), Angstrom^3 *)
let compute_vols (radii : Owl_dense_ndarray_d.arr) : Owl_dense_ndarray_d.arr =
    Owl_dense_ndarray_d.map sphere_volume radii

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
        let radii_cache = Cache.make (fun () -> compute_radii elms) in
        let vols_cache =
            Cache.make (fun () -> compute_vols (Cache.force radii_cache)) in
        { coords; angles = a_cache; r = r_cache; radii = radii_cache;
          vols = vols_cache; elms; name }

(** Flatten a (1, n) row into a genuine 1-D (n,) array. *)
let flatten_row (a : Owl_dense_ndarray_d.arr) : Owl_dense_ndarray_d.arr =
    Owl_dense_ndarray_d.reshape a [| (Owl_dense_ndarray_d.shape a).(1) |]

let theta (m : t) : Owl_dense_ndarray_d.arr =
    flatten_row (Owl_dense_ndarray_d.get_slice [[0]; []] (Cache.force m.angles))

let phi (m : t) : Owl_dense_ndarray_d.arr =
    flatten_row (Owl_dense_ndarray_d.get_slice [[1]; []] (Cache.force m.angles))

let r (m : t) : Owl_dense_ndarray_d.arr = flatten_row (Cache.force m.r)

let coords (m : t) : Owl_dense_ndarray_d.arr = m.coords

let radii (m : t) : Owl_dense_ndarray_d.arr = flatten_row (Cache.force m.radii)

let vols (m : t) : Owl_dense_ndarray_d.arr = flatten_row (Cache.force m.vols)

let elms (m : t) : string array = m.elms

let name (m : t) : string = m.name
