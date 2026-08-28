(* two independent tokenizers live here:
   - `token`: the ion-string grammar (db key format, e.g. "fe3+"), used by
     Parser.ml/Db.ml for query construction.
   - `atom_count_line`/`skip_line`/`xyz_atom_line`: an xyz-file scanner. An
     xyz file is: <natoms>\n<comment>\n(<elem> <x> <y> <z>\n){natoms}. We
     only ever care about the element field -- the atom count and comment
     line are read and discarded whole, and coordinates are consumed
     without ever being run through the periodic-table grammar (a comment
     line or a coordinate value could otherwise spuriously match `elem`,
     e.g. the word "in" matching Indium). Line structure, not token
     lookahead, is what keeps those separate. *)

{   open Toks   }

(* regex for charge (can be <num>+ or <num>- ) *)
let nonzero = ['1'-'9']
let num     = (nonzero+)
let sign    = '+' | '-'

(* not lexer's job to decide if the charge is valid or not *)
let alpha   = ['a'-'z']
let elem    = alpha | (alpha alpha)

(* proper-case element symbol as written in an xyz file, e.g. "Fe", "H" *)
let alpha_up = ['A'-'Z']
let alpha_lo = ['a'-'z']
let elem_ci  = alpha_up (alpha_lo)?

let digit  = ['0'-'9']
let uint   = digit+

(* control chars *)
let eoln    = '\n' | "\r\n"
let white   = [' ' '\r' '\t']

rule token = parse

    (* skip whitespaces *)
    | [' ' '\r' '\t' '\n']+ {token lexbuf}

    (* match sign of charge *)
    | (sign as s)           {Sign(s)}

    (* match charge *)
    | (num as n)            {Charge(n)}

    (* match valid ion as <elem> *)
    | (elem as e)           {Elem(e)}

    (* EOF *)
    | eof                   {EOF}

    (* Token error (bad input) *)
    | _ as err              {Token_error(err)}

(* line 1 of an xyz file: the atom count, nothing else on the line matters. *)
and atom_count_line = parse
    | white* (uint as n) [^ '\n']* (eoln | eof)  { int_of_string n }
    | white* eof                                 { 0 }

(* consumes one whole line (line 2, the free-text comment/title) and
   discards it -- never tokenized against the element grammar. *)
and skip_line = parse
    | [^ '\n']* (eoln | eof)  { () }

(* one atom line: "<elem> <x> <y> <z>". Only the leading element field is
   captured; everything after it on the line (the coordinates) is consumed
   and thrown away in the same match, so a coordinate like "-3.14" is never
   handed to a rule that expects element/charge syntax. *)
and xyz_atom_line = parse
    | white* (elem_ci as e) [^ '\n']* (eoln | eof)  { Elem (String.lowercase_ascii e) }
    | white* eof                                    { EOF }

{   }
