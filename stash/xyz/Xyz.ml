(** Streams an xyz file one atom at a time; never materializes a 
per-atom list. See {!Lexer}'s header for why this is line-structured. *)

(** Folds over every atom's element symbol in file order.
@param f accumulator update, applied once per atom
@param init initial accumulator value
@param lexbuf xyz file content, positioned at the start
@return the final accumulator after every atom line is consumed *)
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

(** Distinct elements present, O(distinct elements) memory.
    @param lexbuf xyz file content, positioned at the start
    @return each element symbol that appears at least once, unordered *)
let unique_elements (lexbuf : Lexing.lexbuf) : string list =
    let seen : (string, unit) Hashtbl.t = Hashtbl.create 16 in
    fold (fun () e -> if not (Hashtbl.mem seen e) then Hashtbl.add seen e ()) () lexbuf;
    Hashtbl.fold (fun e () acc -> e :: acc) seen []
