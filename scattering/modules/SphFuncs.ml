open Gsl
open Printf

exception SphHarm_InvalidLMax of int
exception SphHarm_ThetaPhiInvalidDim of string
exception SphHarm_ThetaPhiUnequal of string
exception SphHarm_EmptyThetaOrPhi of string

exception SphBess_InvalidLMax of string
exception SphBess_EmptyRadii of string
exception SphBess_EmptyQGrid of string
exception SphBess_InvalidRadii of string
exception SphBess_InvalidQ of string

module Owl_dense_ndarray_z = Owl_dense_ndarray_z

(** Computes the complex spherical harmonics Y_l^m(theta, phi)
    @param l Degree of the harmonic (l >= 0)
    @param m Order of the harmonic (-l <= m <= l)
    @param theta Polar angle in [0, pi]
    @param phi Azimuthal angle in [0, 2pi) *)
let sphHarm (lMax : int) (theta : Owl_dense_ndarray_d.arr) 
    (phi : Owl_dense_ndarray_d.arr) : Owl_dense_ndarray_z.arr =

    (* get size of theta/phi and dimensions *)
    let theta_len = Owl_dense_ndarray_d.numel theta in
    let theta_dim = Array.length (Owl_dense_ndarray_d.shape theta) in
    let phi_len = Owl_dense_ndarray_d.numel phi in
    let phi_dim = Array.length (Owl_dense_ndarray_d.shape phi) in

    (* assertion checks *)
    if lMax < 0 then 
        raise(SphHarm_InvalidLMax(lMax))
    ;
    if lMax == 0 then
        print_endline
            "Warning: lMax is set to zero! Only the l=0, m=0 mode will be computed."
    ;
    if theta_dim <> 1 || phi_dim <> 1 then
        let msg =
            sprintf "Theta dim: %d, Phi dim: %d" theta_dim phi_dim
        in
        raise (SphHarm_ThetaPhiInvalidDim msg)
    else ()
    ;
    if theta_len == 0 || phi_len == 0 then
        let msg =
            sprintf "Theta len: %d, Phi len: %d" (theta_len) (phi_len)
        in
        raise(SphHarm_EmptyThetaOrPhi(msg))
    else ()
    ;
    if theta_len <> phi_len then
        let msg = sprintf "Theta len: %d, Phi len: %d" (theta_len) (phi_len) in
        raise (SphHarm_ThetaPhiUnequal msg)
    else ()
    ;

    (* create empty complex array *)
    let y = Owl_dense_ndarray_z.zeros [|((lMax + 1) * (lMax + 2) / 2); theta_len|] in

    let rec helper (y:Owl_dense_ndarray_z.arr) (m : int) (l : int) (max : int) (i : int) 
        (theta: Owl_dense_ndarray_d.arr) (phi : Owl_dense_ndarray_d.arr) : unit =
        
        if l > max then
            ()
        else if m > l then
            helper y 0 (l+1) max i theta phi
        else

            (* get angles *)
            let theta_i = Owl_dense_ndarray_d.get theta [|i|] in
            let phi_i = Owl_dense_ndarray_d.get phi [|i|] in

            (* compute normalized P_l^m amplitude *)
            let ampl = Complex.{
                re = Sf.legendre_sphPlm l m (cos theta_i);
                im = 0.0
            }
            in

            (* compute complex phase exponential: e^(i * m * phi) *)
            let m_float = float_of_int m in
            let phase = Complex.{
                re = cos (m_float *. phi_i);
                im = sin (m_float *. phi_i)
            }
            in

            (* spherical harmonic =  P_l^m(theta) *  e^(i * m * phi) *)
            let sphh = Complex.mul ampl phase in

            (* indexed: [ (l=0,m=0), (l=1,m=0), (l=1,m=1), ..., (l=lMax, m=lMax)] *)
            let idx = (l*(l+1)/2) + m in
            Owl_dense_ndarray_z.set y [|idx; i|] sphh;

            (* recurse *)
            helper y (m+1) l max i theta phi 
    in

    for i = 0 to (theta_len - 1) do
        helper y 0 0 lMax i theta phi 
    done;

    y

(** Computes the spherical Bessel functions j_l for all orders l = 0..lMax over a q-grid
and radii. For reference, j_0(x) = sin(x)/x is the Debye kernel. *)
let sphBess (r : Owl_dense_ndarray_d.arr) (q : Owl_dense_ndarray_d.arr) 
    (lMax : int) : Owl_dense_ndarray_d.arr =
    
    if Owl_dense_ndarray_d.numel r == 0 then
        raise(SphBess_EmptyRadii("radii cannot be empty"))
    ;
    if Owl_dense_ndarray_d.numel q == 0 then
        raise(SphBess_EmptyQGrid("q grid cannot be empty"))
    ;
    if Owl_dense_ndarray_d.exists (fun e -> e < 0.0) r then
        raise(SphBess_InvalidRadii("radii must be >= 0"))
    ;
    if Owl_dense_ndarray_d.exists (fun e -> e < 0.0) q then
        raise(SphBess_InvalidQ("q points must be >= 0"))
    ;
    if lMax < 0 then
        raise(SphBess_InvalidLMax("lMax must be >= 0"))
    ;

    let q_len = Owl_dense_ndarray_d.numel q in
    let r_len = Owl_dense_ndarray_d.numel r in

    (* broadcast q * r *)
    let q_col = Owl_dense_ndarray_d.expand q 2 ~hi:true in (* Shape [|Q; 1|] *)
    let r_row = Owl_dense_ndarray_d.expand r 2 in          (* Shape [|1; N|] *)
    let qr    = Owl_dense_ndarray_d.mul q_col r_row in     (* Broadcasts to [|Q; N|] *)

    (* preallocate output: shape [|lMax+1; Q; N|] *)
    let j = Owl_dense_ndarray_d.zeros [|lMax+1; q_len; r_len|] in

    (* for each (q,r) point, one GSL call gives every l = 0..lMax *)
    for q_idx = 0 to q_len - 1 do
        for r_idx = 0 to r_len - 1 do
            let x = Owl_dense_ndarray_d.get qr [|q_idx; r_idx|] in
            for l = 0 to lMax do
                Owl_dense_ndarray_d.set j [|l; q_idx; r_idx|] (Gsl.Sf.bessel_jl l x)
            done
        done
    done;

    j
