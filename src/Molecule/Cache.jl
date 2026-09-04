"A memoized value guarded by a lock; safe to force concurrently."
module Cache

export Lazy, make, force

"A cache of a value of type `T`. The thunk runs on the first [`force`](@ref)."
mutable struct Lazy{T}
    const f::Any
    const lock::ReentrantLock
    done::Bool
    value::T
    Lazy{T}(f) where {T} = new{T}(f, ReentrantLock(), false)
end

"""
    make(::Type{T}, f) -> Lazy{T}

Build an unforced cache of the `T`-valued thunk `f`.

# Arguments
- `::Type{T}`: element type the thunk returns.
- `f`: zero-arg thunk, run once on the first [`force`](@ref).
"""
make(::Type{T}, f) where {T} = Lazy{T}(f)

"""
    force(c::Lazy{T}) -> T

Run `c`'s thunk once, under its lock, then return the stored value on every call.

# Arguments
- `c`: the cache to force.
"""
function force(c::Lazy{T})::T where {T}
    @lock c.lock begin
        if !c.done
            c.value = c.f()::T
            c.done = true
        end
    end
    return c.value
end

end # module
