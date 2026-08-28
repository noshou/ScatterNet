(** Tokens shared by both {!Lexer} entry points. *)
type tok =
    | Token_error of char  (** unrecognized character *)
    | Elem of string       (** element symbol, lowercase *)
    | Sign of char         (** '+' or '-' *)
    | Charge of string     (** charge magnitude digits *)
    | EOF
