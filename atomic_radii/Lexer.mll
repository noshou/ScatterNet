(** Ion-string grammar (db key format, e.g. ["fe3+"]), used by {!Parser}. *)

{   open Toks   }

let nonzero = ['1'-'9']
let num     = nonzero+
let sign    = '+' | '-'

let alpha   = ['a'-'z']
let elem    = alpha | alpha alpha

(** One token per call.
    @param lexbuf ion-string content
    @return the next [tok]; [Token_error] on an unrecognized character *)
rule token = parse
    | [' ' '\r' '\t' '\n']+ { token lexbuf }
    | (sign as s)           { Sign s }
    | (num as n)            { Charge n }
    | (elem as e)           { Elem e }
    | eof                   { EOF }
    | _ as err              { Token_error err }

{   }
