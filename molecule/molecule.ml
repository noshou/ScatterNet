
open Helpers

module Nd = Owl_base_dense_ndarray_d

(** Raised when coordinates, radii, or their shapes are invalid. *)
exception Molecule_error of string

type t = {
    coords : Nd.arr;                     (** (n, 3): x, y, z, centered at the centroid. *)
    angles : (Nd.arr * Nd.arr) Cache.t;  (** (theta, phi), each (n, 1), radians. *)
    r      : Nd.arr Cache.t;             (** (n, 1): radial distance from the centroid. *)
    vols   : Nd.arr Cache.t;             (** (n, 1): per-atom excluded volume. *)
    elms   : string array;               (** element symbol per atom. *)
    name   : string;                     (** molecule name. *)
}

(* let init 
(name : string) 
(elms : string array) 
(coords : (float * float * float) array) : t = *)

(* let clip_by_value ~amin ~amax = Nd.map (fun v -> Float.max amin (Float.min amax v))
let atan2 y x = Nd.map2 (fun yi xi -> Float.atan2 yi xi) y x *)


(** Centers coordinates at their centroid.
    @param coords per-atom (x, y, z)
    @return (n, 3), centered
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
            let centered = Nd.zeros [| Array.length coords; 3 |] in
            Array.iteri (fun i (x, y, z) ->
                Nd.set centered [| i; 0 |] (x -. mu_x);
                Nd.set centered [| i; 1 |] (y -. mu_y);
                Nd.set centered [| i; 2 |] (z -. mu_z)
            ) coords;
            centered
        end

(** Radial distance of each atom from the centroid.
    @param coords (n, 3), centered
    @return (n, 1)
    @raise Molecule_error if [coords] is empty or not (n, 3) *)
let radii (coords : Nd.arr) : Nd.arr =
    let shape = Nd.shape coords in
    if  shape.(0) = 0 then
        raise (Molecule_error "Empty coordinates")
    else if shape.(1) <> 3 then
        raise (Molecule_error "coords must have exactly 3 columns")
    else
        begin 
            let x = Nd.get_slice [ []; [0] ] coords in
            let y = Nd.get_slice [ []; [1] ] coords in
            let z = Nd.get_slice [ []; [2] ] coords in
            Nd.sqrt (Nd.add (Nd.add (Nd.mul x x) (Nd.mul y y)) (Nd.mul z z))
        end

(** Polar (theta) and azimuthal (phi) angle of each atom, physics
    convention: theta = arccos(z/r) in \[0, pi\], phi = atan2(y, x).
    @param coords (n, 3), centered
    @param r (n, 1), radial distance, from {!radii}
    @return (theta, phi), each (n, 1), radians
    @raise Molecule_error if shapes are empty, wrong, or mismatched *)
let angles (coords : Nd.arr) (r : Nd.arr) : Nd.arr * Nd.arr =
    
    (* assertion checks *)
    let (shape, shaper) = (Nd.shape coords, Nd.shape r) in
    if shape.(0) = 0 || shaper.(0) = 0 then
        raise (Molecule_error "Empty coordinates or radii")
    else if shape.(1) <> 3 then
        raise (Molecule_error "coords must have exactly 3 columns")
    else if shaper.(1) <> 1 then
        raise (Molecule_error "r must have exactly 1 column")
    else if shape.(0) <> shaper.(0) then
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
                    (Nd.div (Nd.get_slice [[]; [2]] coords) r_safe)
                )
            in    
            let phi = Nd.atan2 (
                Nd.get_slice[[]; [1]] coords) (Nd.get_slice[[]; [0]] coords
                ) 
            in

            (theta, phi)
        end