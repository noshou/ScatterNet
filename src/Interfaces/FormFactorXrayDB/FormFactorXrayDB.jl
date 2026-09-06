"""
X-ray form factors `f(q,E) = f0(s) + f1(E) + i f2(E)` from the Python `xraydb`
package. `py/FormFact_py.py` (bundled alongside this file) holds the tier logic.

Everything here is pure Julia; the one function that actually crosses into
Python, [`compute_form_factors`](@ref), is supplied by the `FormFactorXrayDBExt`
package extension and only exists once `PythonCall` is loaded. Consumers that
just want the scattering geometry therefore never pull in a Python stack.

The public trio ([`Interfaces.form_factor_table`](@ref) /
[`Interfaces.form_factors`](@ref) / [`Interfaces.form_factor_log`](@ref)) are
methods on the `Interfaces` generics; nothing here is meant to be reached as
`FormFactorXrayDB.<name>` from outside `Interfaces` except through the parent
facade.
"""
module FormFactorXrayDB

import ..Interfaces
using  ..Interfaces: FormFactorSource

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

# Arguments
- `ions`: vector of ion strings.
- `energy`: photon energy in eV.
- `qvals`: vector of q values in Å⁻¹.
"""
compute_form_factors(args...) = throw(FormFactorError(
    "the xraydb backend is not loaded; run `using PythonCall` to activate the " *
    "FormFactorXrayDBExt extension (and make sure PythonCall is in your environment)"))

"""
    Interfaces.form_factor_table([src::FormFactorSourceXrayDB,] energy::Real, ions, qvals) -> FF

Build an [`FF`](@ref) container for `ions` at one `energy` (eV) over the `qvals`
(Å⁻¹) grid. Thin wrapper over [`compute_form_factors`](@ref). The `src`-less
form defaults the backend to `FormFactorSourceXrayDB()`.

# Arguments
- `src`: the xraydb backend marker (optional).
- `energy`: photon energy in eV.
- `ions`: vector of ion strings.
- `qvals`: vector of q values in Å⁻¹.
"""
Interfaces.form_factor_table(::FormFactorSourceXrayDB, energy::Real, ions, qvals)::FF =
    compute_form_factors(collect(String, ions), energy, collect(Float64, qvals))

Interfaces.form_factor_table(energy::Real, ions, qvals)::FF =
    Interfaces.form_factor_table(FormFactorSourceXrayDB(), energy, ions, qvals)

"""
    Interfaces.form_factor_log(t::FF) -> Vector{String}

Construction-time log (one line per non-full ion).
"""
Interfaces.form_factor_log(t::FF)::Vector{String} = t.log

"""
    Interfaces.form_factors(t::FF, ions, qvals) -> Matrix{ComplexF64}

Per-ion form-factor rows from container `t` as a `(length(ions), length(qvals))`
matrix: row `i` is the form factor for `ions[i]`, its columns aligned to `qvals`
in the order given. Pass the per-atom ion vector and the result drops straight
into [`compute_B_lm`](@ref) as `f_atoms` — no intermediate mapping step.

Throws [`FormFactorError`](@ref) if a queried q is not one of `t`'s grid points
(exact match), or if any `ions[i]` is absent from `t` (a dummy site or an ion
the backend could not resolve — see [`Interfaces.form_factor_log`](@ref)).

# Arguments
- `t`: form-factor container from [`Interfaces.form_factor_table`](@ref) / [`compute_form_factors`](@ref).
- `ions`: ion strings to fetch, one per output row (typically one per atom).
- `qvals`: vector of q values; each must match a grid point of `t` exactly.
"""
function Interfaces.form_factors(t::FF, ions::AbstractVector{<:AbstractString},
                                 qvals::AbstractVector{<:Real})::Matrix{ComplexF64}
    idx = Vector{Int}(undef, length(qvals))
    @inbounds for i in eachindex(qvals)
        q = Float64(qvals[i])
        j = get(t.qmp, q, 0)
        j == 0 && throw(FormFactorError("q=$q is not one of this container's q-points"))
        idx[i] = j
    end
    f = Matrix{ComplexF64}(undef, length(ions), length(qvals))
    @inbounds for (a, ion) in enumerate(ions)
        row = get(t.tbl, String(ion), nothing)
        row === nothing && throw(FormFactorError("ion \"$ion\" is not in this container"))
        for (c, k) in enumerate(idx)
            f[a, c] = row[k]
        end
    end
    return f
end

end # module
