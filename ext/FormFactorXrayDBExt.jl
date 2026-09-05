"""
Package extension supplying the one part of `FormFactorXrayDB` that crosses the
language boundary: [`compute_form_factors`](@ref), which calls `py/FormFact_py.py`
through PythonCall. Loads automatically once `PythonCall` is in the session.

Keeping this out of the package proper means an environment that never touches
form factors.
"""
module FormFactorXrayDBExt

using PythonCall: pyimport, pylist, pyconvert, Py
using ScatterNet.Scattering.FormFactorXrayDB: FormFactorXrayDB, FF, FormFactorError

const _PYMOD = Ref{Py}()

"""
    _pymod() -> Py

Import and memoize the `FormFact_py` Python module, prepending the package `py/`
dir to `sys.path` on first call. Throws [`FormFactorError`](@ref) if the import
fails.
"""
function _pymod()::Py
    if !isassigned(_PYMOD)
        sys = pyimport("sys")
        d = pkgdir(FormFactorXrayDB, "py")
        d ∉ sys.path && sys.path.insert(0, d)
        try
            _PYMOD[] = pyimport("FormFact_py")
        catch e
            throw(FormFactorError("cannot import FormFact_py ($e); build the CondaPkg env"))
        end
    end
    return _PYMOD[]
end

"""
    compute_form_factors(ions, energy::Real, qvals) -> FF

The real xraydb-backed implementation. More specific than the package's
catch-all `Vararg` stub, so loading this extension makes it the one that
dispatch selects.
"""
function FormFactorXrayDB.compute_form_factors(
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

end # module
