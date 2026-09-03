"""
Atomic/ionic radii: parse an ion string, then resolve its radius through a
fallback chain over the bundled `atomic_radii.sqlite3` (loaded once into `Dict`s).
"""
module AtomicRadii

using SQLite: SQLite
using DBInterface: DBInterface
import ..Interfaces
using ..Interfaces: RadiiSource

export Ion, tryparse_ion, ion_key, resolve_one, lookup, AtomicRadiiSource

# ---- ion string -> (element, signed charge) -------------------------------

"A parsed ion, e.g. `\"fe3+\"` -> `Ion(\"fe\", 3)`."
struct Ion
    element::String
    charge::Int
end

# element = 1-2 lowercase; magnitude has a non-zero leading digit (no charge 0);
# a bare sign means +/-1.
const _ION_RE = r"^\s*([a-z]{1,2})\s*(?:([1-9][0-9]*)\s*([+-])|([+-])\s*([1-9][0-9]*)?)?\s*$"

"Parse an ion string; `nothing` for a bare element like `\"fe\"`."
function tryparse_ion(s::AbstractString)::Union{Ion,Nothing}
    m = match(_ION_RE, s)
    m === nothing && return nothing
    elem = String(m.captures[1]::AbstractString)
    sign = m.captures[3] === nothing ? m.captures[4] : m.captures[3]
    sign === nothing && return nothing
    magcap = m.captures[2] === nothing ? m.captures[5] : m.captures[2]
    mag = magcap === nothing ? 1 : parse(Int, magcap)
    return Ion(elem, sign == "-" ? -mag : mag)
end

"`Ion` -> `ionic_radii` key, digits-then-sign (`Ion(\"fe\", 3)` -> `\"fe3+\"`)."
function ion_key(ion::Ion)::String
    ion.charge == 0 && throw(ArgumentError("charge 0 has no ion-string form"))
    string(ion.element, abs(ion.charge), ion.charge > 0 ? '+' : '-')
end

# ---- static tables -------------------------------------------------------

const _IONIC   = Dict{String,Float64}()               # "fe3+" => radius Å
const _ATOMIC  = Dict{String,Tuple{Float64,String}}() # "fe"   => (radius Å, type)
const _CHARGES = Dict{String,Vector{Int}}()           # "fe"   => sorted charges

_dbpath()::String = joinpath(pkgdir(@__MODULE__)::String, "data", "atomic_radii.sqlite3")

function _load!(path::String = _dbpath())
    empty!(_IONIC); empty!(_ATOMIC); empty!(_CHARGES)
    isfile(path) || error("AtomicRadii: missing data file $path")
    db = SQLite.DB(path)
    try
        DBInterface.execute(db, "PRAGMA query_only = ON;")
        for row in DBInterface.execute(db, "SELECT ion, radius FROM ionic_radii")
            _IONIC[String(row.ion)] = Float64(row.radius) / 100.0   # pm -> Å
        end
        for row in DBInterface.execute(db, "SELECT element, radius, radius_type FROM atomic_radii")
            _ATOMIC[String(row.element)] = (Float64(row.radius), String(row.radius_type))
        end
        for row in DBInterface.execute(db, "SELECT element, charge FROM element_charges")
            push!(get!(() -> Int[], _CHARGES, String(row.element)), Int(row.charge))
        end
        for v in values(_CHARGES); sort!(v); end
    finally
        DBInterface.close!(db)
    end
    return nothing
end

__init__() = _load!()

ion_radius(key::AbstractString)::Union{Float64,Nothing} = get(_IONIC, key, nothing)
element_radius(el::AbstractString)::Union{Tuple{Float64,String},Nothing} = get(_ATOMIC, el, nothing)

"Closest charge state on file for `element` to `charge`; `nothing` if none."
function nearest_ion(element::AbstractString, charge::Int)::Union{String,Nothing}
    cs = get(_CHARGES, element, nothing)
    cs === nothing && return nothing
    best = cs[1]
    for c in cs
        abs(c - charge) < abs(best - charge) && (best = c)
    end
    return string(element, abs(best), best > 0 ? '+' : '-')
end

"""
Resolve one ion/element to a radius (Å): exact charge match, else nearest
charge state, else bare element, else `nothing`. Unparseable strings go
straight to the bare-element table.
"""
function resolve_one(ion::AbstractString)::Union{Float64,Nothing}
    parsed = tryparse_ion(ion)
    if parsed === nothing
        er = element_radius(ion)
        return er === nothing ? nothing : er[1]
    end
    exact = ion_radius(ion_key(parsed))
    exact === nothing || return exact
    near = nearest_ion(parsed.element, parsed.charge)
    if near !== nothing
        nr = ion_radius(near)
        nr === nothing || return nr
    end
    er = element_radius(parsed.element)
    return er === nothing ? nothing : er[1]
end

struct _Miss end
const _MISS = _Miss()

"Batch resolve; input order and count preserved, repeats deduped per call."
function lookup(ions::AbstractVector{<:AbstractString})
    cache = Dict{String,Union{Float64,Nothing}}()
    out = Vector{Tuple{String,Union{Float64,Nothing}}}(undef, length(ions))
    @inbounds for i in eachindex(ions)
        s = String(ions[i])
        hit = get(cache, s, _MISS)
        r::Union{Float64,Nothing} = hit isa _Miss ? resolve_one(s) : hit
        hit isa _Miss && (cache[s] = r)
        out[i] = (s, r)
    end
    return out
end

"Concrete [`Interfaces.RadiiSource`](@ref) backed by [`lookup`](@ref)."
struct AtomicRadiiSource <: RadiiSource end
Interfaces.lookup(::AtomicRadiiSource, ions::AbstractVector{<:AbstractString}) = lookup(ions)

end # module
