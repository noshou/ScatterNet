(** WIP, stashed: xyz-file line scanning, pulled out of atomic_radii/Lexer.mll.
    Only reads the element field per atom line; coordinates are consumed and
    discarded, not parsed. Real coordinate parsing (float triples) and the
    polar-coordinate conversion still need to be written before this is
    useful again -- see stash/xyz/README.md.

    To revive: paste these rules (and the regexes they use) back into a
    lexer file alongside atomic_radii's ion-string `token` rule, or give
    this its own dune library once coordinate parsing exists. *)

{   open Toks   }

let alpha_up = ['A'-'Z']
let alpha_lo = ['a'-'z']
let elem_ci  = alpha_up alpha_lo?  (* proper-case symbol, e.g. "Fe" *)

let digit = ['0'-'9']
let uint  = digit+

let eoln  = '\n' | "\r\n"
let white = [' ' '\r' '\t']

(** Line 1 of an xyz file.
    @param lexbuf xyz file content, positioned at the start
    @return the declared atom count; [0] if the line is empty/missing *)
rule atom_count_line = parse
    | white* (uint as n) [^ '\n']* (eoln | eof) { int_of_string n }
    | white* eof                                { 0 }

(** Discards one whole line (the xyz comment/title line).
    @param lexbuf position at the start of the line to discard *)
and skip_line = parse
    | [^ '\n']* (eoln | eof) { () }

(** One atom line, [<elem> <x> <y> <z>]; coordinates are consumed and
    discarded in the same match, never lexed as tokens -- this is exactly
    what needs to change once real coordinate parsing is written: the
    [^ '\n']* here should become a float-triple grammar instead of a
    throwaway match.
    @param lexbuf position at the start of an atom line
    @return [Elem symbol] (lowercased), or [EOF] if no line remains *)
and xyz_atom_line = parse
    | white* (elem_ci as e) [^ '\n']* (eoln | eof) { Elem (String.lowercase_ascii e) }
    | white* eof                                   { EOF }

{   }
