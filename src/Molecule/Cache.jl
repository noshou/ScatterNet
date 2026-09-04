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

"Build an unforced cache of a `T`-valued thunk."
make(::Type{T}, f) where {T} = Lazy{T}(f)

"Force a cache: run the thunk once, then return the stored value on every call."
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
