open Printf

exception Form_factor_xraydb_error of string

type ff = {
    tbl : (string, Owl_dense_ndarray_z.arr) Hashtbl.t;
    qmp : (float, int) Hashtbl.t;
    log : string array;
}

(** Initialize python *)
let () =
    (* start python interperter *)
    Py.initialize ()

let py_mod : Py.Object.t =
    let sys = Py.import "sys" in
    let path = Py.Module.get sys "path" in
    ignore (Py.Object.call_method path "insert"
        [| Py.Int.of_int 0; Py.String.of_string "form_factor_xraydb" |]);
    try
        Py.import "FormFact_py"
    with
        Py.E (_, errvalue) ->
            raise (Form_factor_xraydb_error (
                sprintf "failed to load FormFact_py.py (%s) - is python3's numpy/xraydb \
                installed? Try: pip install -r form_factor_xraydb/requirements.txt"
            (Py.Object.to_string errvalue)))

let compute_form_factors
    (ions : string array) (energy : float)
    (qvals : Owl_dense_ndarray_d.arr) : ff =

    Py.Gil.with_lock (
        fun () ->
            try
                (* convert ocaml types to python types *)
                let ions_py   = Py.List.of_array_map Py.String.of_string ions in
                let energy_py = Py.Float.of_float energy in
                let qvals_py  = Numpy.of_bigarray qvals in

                (* get the function, call it, and unpack the result *)
                let func = Py.Module.get_function py_mod "compute_form_factors" in
                let res = func [| ions_py; energy_py; qvals_py |] in
                let (tuples_py, log_py) = (Py.Tuple.get res 0, Py.Tuple.get res 1) in

                (* tuples_py is a list[(string, np.ndarray)]; map it into tbl *)
                let (tbl : (string, Owl_dense_ndarray_z.arr) Hashtbl.t) =
                    Hashtbl.create (Py.List.size tuples_py)
                in
                ignore (tuples_py |> Py.List.fold_left (
                    fun t kv ->
                        let (k, v) =
                            (Py.String.to_string (Py.Tuple.get kv 0),
                            Numpy.to_bigarray Bigarray.complex64 Bigarray.c_layout
                                (Py.Tuple.get kv 1))
                        in
                        Hashtbl.replace t k v;
                        t
                ) tbl);

                (* convert log_py to a string array *)
                let log = Py.List.to_array_map Py.String.to_string log_py in

                (* map each qval to its index, so a single point can be picked
                   back out of a form factor array later *)
                let num_qvals = Owl_dense_ndarray_d.numel qvals in
                let (qmp : (float, int) Hashtbl.t) = Hashtbl.create num_qvals in
                Owl_dense_ndarray_d.iteri (fun i q -> Hashtbl.replace qmp q i) qvals;

                { tbl; qmp; log }

            with
                Py.E (_, errvalue) ->
                    raise (Form_factor_xraydb_error (
                        sprintf "compute_form_factors failed (%s)"
                            (Py.Object.to_string errvalue)
                        )
                    )
    )

module FormFactorSourceXrayDB = struct

    type t = ff

    let log (t : t) : string array = t.log

    let create
        (energy : float) (ions : string list)
        (qvals : Owl_dense_ndarray_d.arr) : t =
        compute_form_factors (Array.of_list ions) energy qvals

    let lookup
        (t : t) (ions : string list) (qvals : Owl_dense_ndarray_d.arr)
        : (string * Owl_dense_ndarray_z.arr) list =
        let n = Owl_dense_ndarray_d.numel qvals in
        let indices = Array.init n (fun i ->
            let q = Owl_dense_ndarray_d.get qvals [| i |] in
            match Hashtbl.find_opt t.qmp q with
            | Some idx -> idx
            | None ->
                raise (
                    Form_factor_xraydb_error (
                        sprintf "lookup: q=%f is not one of this container's q-points" q
                    )
                )
        ) in
        List.filter_map (fun ion ->
            match Hashtbl.find_opt t.tbl ion with
            | None -> None
            | Some arr ->
                let out = Owl_dense_ndarray_z.zeros [| n |] in
                Array.iteri (fun j idx ->
                    let v = Owl_dense_ndarray_z.get arr [| idx |] in
                    Owl_dense_ndarray_z.set out [| j |] v
                ) indices;
                Some (ion, out)
        ) ions
end
