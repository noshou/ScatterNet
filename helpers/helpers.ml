(** Small reusable utilities shared across this project's libraries. *)

module Cache = Cache

let repo_root : string =
    let exe = Sys.executable_name in
    let parts = String.split_on_char '/' exe in
    let rec take_until_build acc = function
        | [] -> List.rev acc
        | "_build" :: _    -> List.rev acc
        | part     :: rest -> take_until_build (part :: acc) rest
    in
    String.concat "/" (take_until_build [] parts)
