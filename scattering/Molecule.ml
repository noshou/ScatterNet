open Helpers

module Nd = Owl_base_dense_ndarray_d

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
                | (elem, Some r) ->
                    (elem, Some (r, (4. /. 3.) *. Float.pi *. r ** 3.))
                | (elem, None) ->
                    (elem, None)
        ) rads
end

(** Volume lookup backed by the real {!Atomic_radii.AtomicRadiiSource}. *)
module AV = AtomicVolume (Atomic_radii.AtomicRadiiSource)


(** Raised when coordinates, radii, or their shapes are invalid. *)
exception Molecule_error of string

type t = {
    coords : Nd.arr;            (** (3, n): x, y, z, centered at the centroid. *)
    angles : Nd.arr Cache.t;    (** (2, n): theta, phi in radians. *)
    r      : Nd.arr Cache.t;    (** (1, n): radial distance from the centroid. *)
    vols   : Nd.arr Cache.t;    (** (1, n): per-atom excluded volume. *)
    elms   : string array;      (** element symbol per atom. *)
    name   : string;            (** molecule name. *)
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

(** Radial distance of each atom from the centroid.
    @param coords (3, n), centered
    @return (1, n)
    @raise Molecule_error if [coords] is empty or not (3, n) *)
let radii (coords : Nd.arr) : Nd.arr =
    let shape = Nd.shape coords in
    if shape.(0) <> 3 then
        raise (Molecule_error "coords must have exactly 3 rows")
    else if shape.(1) = 0 then
        raise (Molecule_error "Empty coordinates")
    else
        begin
            let x = Nd.get_slice [ [0]; [] ] coords in
            let y = Nd.get_slice [ [1]; [] ] coords in
            let z = Nd.get_slice [ [2]; [] ] coords in
            Nd.sqrt (Nd.add (Nd.add (Nd.mul x x) (Nd.mul y y)) (Nd.mul z z))
        end

(** Polar (theta) and azimuthal (phi) angle of each atom, physics
    convention: theta = arccos(z/r) in \[0, pi\], phi = atan2(y, x).
    @param coords (3, n), centered
    @param r (1, n), radial distance, from {!radii}
    @return (theta, phi), (2, n), radians
    @raise Molecule_error if shapes are empty, wrong, or mismatched *)
let angles (coords : Nd.arr) (r : Nd.arr) : Nd.arr =
    
    (* assertion checks *)
    let (shape, shaper) = (Nd.shape coords, Nd.shape r) in
    if shape.(0) <> 3 then
        raise (Molecule_error "coords must have exactly 3 rows")
    else if shaper.(0) <> 1 then
        raise (Molecule_error "r must have exactly 1 row")
    else if shape.(1) = 0 || shaper.(1) = 0 then
        raise (Molecule_error "Empty coordinates or radii")
    else if shape.(1) <> shaper.(1) then
        raise (Molecule_error "Dimensions of coordinates and radii must match")

    else
        begin
            (*  Physics convention:
                    theta = polar colatitude (arccos(z/r), 0→π),
                    phi   = azimuthal angle  (arctan2(y,x), 0→2π)
                r=0 only occurs for single-atom molecules (atom is its own centroid).
                j_l(0)=0 for l>0 so the angle is irrelevant; use r_safe to avoid NaN. *)
            let r_safe = Nd.map (fun _r -> if _r > 0.0 then _r else 1.0) r in
            let theta =
                Nd.acos(
                    Nd.clip_by_value ~amin:(-1.) ~amax:1.
                    (Nd.div (Nd.get_slice [[2]; []] coords) r_safe)
                )
            in
            let phi = Nd.atan2 (
                Nd.get_slice[[1]; []] coords) (Nd.get_slice[[0]; []] coords
                )
            in
            Nd.concatenate ~axis:0 [| theta; phi |]
        end

(** Per-atom excluded volume, from atomic/ionic radii looked up by
    element symbol.
    @param elms element symbol per atom
    @return (1, n)
    @raise Molecule_error if [elms] is empty, or any element has no
    radius on file at all *)
let vols (elms : string array) : Nd.arr =
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
        let r_cache = Cache.make (fun () -> radii coords) in
        let angles_cache = Cache.make (fun () -> angles coords (Cache.force r_cache)) in
        let vols_cache = Cache.make (fun () -> vols elms) in
        { coords; angles = angles_cache; r = r_cache; vols = vols_cache; elms; name }