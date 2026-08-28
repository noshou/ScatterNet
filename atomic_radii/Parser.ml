open Toks
open Printf
open Lexing
open Lexer

(* parsed ion: "fe3+" -> { element = "fe"; charge = 3 } *)
type ion = { element : string; charge  : int; }

exception Parse_error of string

let parse_error (msg : string) : 'a = raise (Parse_error msg)

(* combine a charge magnitude ("3", or "" for a bare sign) and anchar into a signed int *)
let combine (n : string) (s : char) : int =
    let magnitude = if n = "" then 1 else int_of_string n in
    match s with
        | '+' -> magnitude
        | '-' -> -magnitude
        | _   -> parse_error (sprintf "invalid sign char %c" s)

let expect_eof (lexbuf : lexbuf) : unit = match token lexbuf with
    | EOF -> ()
    | _ -> parse_error "expected end of input"

(* grammar: elem (charge sign | sign charge | sign) *)
let parse (lexbuf : lexbuf) : ion = match token lexbuf with
    | Elem e -> let charge =
        match token lexbuf with
            | Charge n -> (
                match token lexbuf with
                    | Sign s -> expect_eof lexbuf; combine n s
                    | _ -> parse_error "expected sign after charge"
                )
            | Sign s -> (
                match token lexbuf with
                    | Charge n -> expect_eof lexbuf; combine n s
                    | EOF -> combine "" s (* bare sign; EOF already consumed *)
                    | _ -> parse_error "expected charge or end of input after sign"
                )
            | _ -> parse_error "expected charge or sign after element"
        in
        { element = e; charge }

    | Token_error c -> parse_error (sprintf "unrecognized character %c" c)

    | _ -> parse_error "expected element"

let parse_string (s : string) : ion = parse (from_string s)

(* inverse of parse: {element="fe"; charge=3} -> "fe3+", the ionic_radii.db
   key format. charge=0 has no db representation (every row carries a sign)
   so it's rejected rather than guessed at. *)
let to_string (ion : ion) : string =
    if ion.charge = 0 then
        parse_error "charge 0 has no ion string representation"
    else
        let sign = if ion.charge > 0 then '+' else '-' in
        sprintf "%s%d%c" ion.element (abs ion.charge) sign
