"""
X-ray form factors `f(q,E) = f0(s) + f1(E) + i f2(E)` from the Python `xraydb`
package. `py/FormFact_py.py` holds the tier logic.

Everything here is pure Julia; the one function that actually crosses into
Python, [`compute_form_factors`](@ref), is supplied by the `FormFactorXrayDBExt`
package extension and only exists once `PythonCall` is loaded. Consumers that
just want the scattering geometry therefore never pull in a Python stack.
"""
module FormFactorXrayDB

using ...Interfaces: FormFactorSource

export FF, FormFactorError, FormFactorSourceXrayDB, compute_form_factors

"Raised on any failure building or querying form factors."
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

"""
    compute_form_factors(ions, energy::Real, qvals) -> FF

Form factors for a batch of ions at one `energy` over a q grid; one row per
unique ion, aligned to the returned container's q index.

Implemented by the `FormFactorXrayDBExt` extension, which loads with
`PythonCall`. Without it this catch-all method is the only one defined and it
raises [`FormFactorError`](@ref) -- it is deliberately less specific than the
extension's typed method, so the real implementation wins on dispatch without
redefining anything.

# Arguments
- `ions`: vector of ion strings.
- `energy`: photon energy in eV.
- `qvals`: vector of q values in Å⁻¹.
"""
compute_form_factors(args...) = throw(FormFactorError(
    "the xraydb backend is not loaded; run `using PythonCall` to activate the " *
    "FormFactorXrayDBExt extension (and make sure PythonCall is in your environment)"))

"""
    create(energy::Real, ions, qvals) -> FF

Build an [`FF`](@ref) container for `ions` at one `energy` (eV) over the `qvals`
(Å⁻¹) grid. Thin wrapper over [`compute_form_factors`](@ref).

# Arguments
- `energy`: photon energy in eV.
- `ions`: vector of ion strings.
- `qvals`: vector of q values in Å⁻¹.
"""
create(energy::Real, ions, qvals)::FF =
    compute_form_factors(collect(String, ions), energy, collect(Float64, qvals))

"Construction-time log (one line per non-full ion)."
log(t::FF)::Vector{String} = t.log

"""
    lookup(t::FF, ions, qvals) -> Vector{Tuple{String,Vector{ComplexF64}}}

Rows for `ions` at `qvals`, drawn from container `t`; ions with no data are dropped.
Throws [`FormFactorError`](@ref) if a queried q is not one of `t`'s grid points (exact match).

# Arguments
- `t`: form-factor container from [`create`](@ref) / [`compute_form_factors`](@ref).
- `ions`: vector of ion strings to fetch.
- `qvals`: vector of q values; each must match a grid point of `t` exactly.
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
