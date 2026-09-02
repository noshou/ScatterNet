(** Small reusable utilities shared across this project's libraries. *)

(** A memoized value guarded by a mutex, since {!Stdlib.Lazy.force} alone
    is not safe to call concurrently from multiple domains. *)
module Cache : sig

    (** A cache of a value of type ['a]. *) 
    type 'a t

    (** Buildscache from a thunk. The thunk does not run until the first call to {!force}.
        @param f the computation to memoize
        @return a new, unforced cache *)
    val make : (unit -> 'a) -> 'a t

    (** Force a cache. The first call runs the thunk and stores the result; every 
        later call, from any domain, returns that stored result without recomputing it.
        @param c the cache to force
        @return the memoized value *)
    val force : 'a t -> 'a

end

(** Absolute path to the repo root, derived from {!Stdlib.Sys.executable_name}. *)
val repo_root : string

