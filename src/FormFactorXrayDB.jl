"""
X-ray form factors `f(q,E) = f0(s) + f1(E) + i f2(E)` from the Python `xraydb`
package (via PythonCall). `py/FormFact_py.py` holds the tier logic.
"""
module FormFactorXrayDB

using PythonCall: pyimport, pylist, pyconvert, Py
using ..Interfaces: FormFactorSource

export FF, FormFactorError, FormFactorSourceXrayDB, compute_form_factors

struct FormFactorError <: Exception; msg::String end
Base.showerror(io::IO, e::FormFactorError) = print(io, "FormFactorError: ", e.msg)

"Per-batch form factors: `tbl` ion => row aligned to `qmp` (qval => index), plus a build `log`."
struct FF
    tbl::Dict{String,Vector{ComplexF64}}
    qmp::Dict{Float64,Int}
    log::Vector{String}
end

"Marker for the xraydb backend (single implementation of [`Interfaces.FormFactorSource`](@ref))."
struct FormFactorSourceXrayDB <: FormFactorSource end

const _PYMOD = Ref{Py}()

function _pymod()::Py
    if !isassigned(_PYMOD)
        sys = pyimport("sys")
        d = pkgdir(@__MODULE__, "py")
        d ∉ sys.path && sys.path.insert(0, d)
        try
            _PYMOD[] = pyimport("FormFact_py")
        catch e
            throw(FormFactorError("cannot import FormFact_py ($e); build the CondaPkg env"))
        end
    end
    return _PYMOD[]
end

"Form factors for a batch of ions at `energy` (eV) over `qvals` (Å⁻¹); one row per unique ion."
function compute_form_factors(
    ions::AbstractVector{<:AbstractString}, energy::Real,
    qvals::AbstractVector{<:Real}
)::FF
    qs = collect(Float64, qvals)
    np = pyimport("numpy")
    try
        res = _pymod().compute_form_factors(pylist(collect(String, ions)), Float64(energy),
                                            np.asarray(qs; dtype = np.float64))
        tbl = Dict{String,Vector{ComplexF64}}()
        for kv in res[0]
            tbl[pyconvert(String, kv[0])] = pyconvert(Vector{ComplexF64}, kv[1])
        end
        log_lines = pyconvert(Vector{String}, res[1])
        qmp = Dict{Float64,Int}()
        @inbounds for i in eachindex(qs)
            qmp[qs[i]] = i
        end
        return FF(tbl, qmp, log_lines)
    catch e
        e isa FormFactorError && rethrow()
        throw(FormFactorError("compute_form_factors failed ($e)"))
    end
end

"Build a container for `ions` at one `energy` over a q grid."
create(energy::Real, ions, qvals)::FF =
    compute_form_factors(collect(String, ions), energy, collect(Float64, qvals))

"Construction-time log (one line per non-full ion)."
log(t::FF)::Vector{String} = t.log

"""
    lookup(t, ions, qvals) -> Vector{Tuple{String,Vector{ComplexF64}}}

Rows for `ions` at `qvals`; ions with no data are dropped. Errors if a queried q
is not one of `t`'s grid points (exact match).
"""
function lookup(t::FF, ions::AbstractVector{<:AbstractString}, qvals::AbstractVector{<:Real})
    idx = Vector{Int}(undef, length(qvals))
    @inbounds for i in eachindex(qvals)
        q = Float64(qvals[i])
        j = get(t.qmp, q, 0)
        j == 0 && throw(FormFactorError("q=$q is not one of this container's q-points"))
        idx[i] = j
    end
    out = Tuple{String,Vector{ComplexF64}}[]
    for ion in ions
        key = String(ion)
        row = get(t.tbl, key, nothing)
        row === nothing && continue
        push!(out, (key, ComplexF64[row[k] for k in idx]))
    end
    return out
end

end # module
