type 'a t = {
    value : 'a Lazy.t;
    lock  : Mutex.t;
}

let make (f : unit -> 'a) : 'a t =
    { value = lazy (f ()); lock = Mutex.create () }

let force (c : 'a t) : 'a =
    Mutex.protect c.lock (fun () -> Lazy.force c.value)
