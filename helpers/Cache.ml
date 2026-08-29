(** A memoized value guarded by a mutex, since {!Stdlib.Lazy.force} is not thread safe. *)

type 'a t = {
    value : 'a Lazy.t;
    lock  : Mutex.t;
}

(** Build a cache from a thunk. The thunk does not run until the first
    call to {!force}.
    @param f the computation to memoize
    @return a new, unforced cache *)
let make (f : unit -> 'a) : 'a t =
    { value = lazy (f ()); lock = Mutex.create () }

(** Force a cache. The first call runs the thunk and stores the result;
    every later call, from any domain, returns that stored result without
    recomputing it.
    @param c the cache to force
    @return the memoized value *)
let force (c : 'a t) : 'a =
    Mutex.protect c.lock (fun () -> Lazy.force c.value)
