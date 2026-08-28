(* streams an xyz-format text stream one atom line at a time -- nothing
   here ever materializes a full per-atom list. The atom-count header and
   comment line are read and discarded once; each atom line is lexed,
   folded into the caller's accumulator, and forgotten. See Lexer.mll's
   header comment for why this is line-structured rather than token-soup
   lexing. *)

let fold (f : 'acc -> string -> 'acc) (init : 'acc) (lexbuf : Lexing.lexbuf) : 'acc =
    let natoms = Lexer.atom_count_line lexbuf in
    Lexer.skip_line lexbuf;
    let rec loop n acc =
        if n <= 0 then acc
        else
            match Lexer.xyz_atom_line lexbuf with
            | Toks.Elem e -> loop (n - 1) (f acc e)
            | Toks.EOF -> acc (* file ended early; take what we got *)
            | _ -> failwith "Xyz.fold: malformed atom line"
    in
    loop natoms init

(* distinct elements only, via a Hashtbl used as a set: peak memory is
   O(distinct elements actually present) -- typically a handful, even for
   a file with millions of atoms -- never O(natoms). *)
let unique_elements (lexbuf : Lexing.lexbuf) : string list =
    let seen : (string, unit) Hashtbl.t = Hashtbl.create 16 in
    fold (fun () e -> if not (Hashtbl.mem seen e) then Hashtbl.add seen e ()) () lexbuf;
    Hashtbl.fold (fun e () acc -> e :: acc) seen []
