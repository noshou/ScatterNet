# Exercises the xraydb form-factor backend: src/Interfaces/FormFactorXrayDB/ plus
# the FormFactorXrayDBExt package extension (the Python-facing half, loaded with
# PythonCall). It needs the CondaPkg env (numpy + xraydb).
#
# PYTHON AIRLOCK -- this file is opt-in.
# It is the only part of the suite that can trip CondaPkg into provisioning a
# Python env. The value tests therefore run ONLY when the environment variable
# SCATTERNET_TEST_XRAYDB is set to a truthy value ("1" / "true" / "yes",
# case-insensitive). Without it: PythonCall is never imported here,
# `compute_form_factors` is never called, the conda env is never built, and the
# Python-facing testsets are @test_skip'd. The pure-Julia marker assertions
# (FormFactorSourceXrayDB, FormFactorError) always run.
#
# PythonCall + CondaPkg are not in `test/Project.toml` (loading PythonCall at all
# -- even in a `Pkg.test` precompile worker -- provisions the conda env). They
# live in `test/xraydb/Project.toml`; the opt-in branch below puts that env on
# LOAD_PATH and only then imports PythonCall.

const IFACE = ScatterNet.Interfaces
using ScatterNet.Interfaces.FormFactorXrayDB: compute_form_factors, FF, FormFactorError, FormFactorSourceXrayDB

# tiny truthy-env check
_ff_truthy(v) = lowercase(strip(String(v))) in ("1", "true", "yes")
const _optin = _ff_truthy(get(ENV, "SCATTERNET_TEST_XRAYDB", ""))

check_c(a, b) = abs(a - b) < 1e3 * DEFAULT_ATOL   # looser: against tabulated xraydb reference values, not a closed-form identity
qvals = [0.1, 0.2]
qgrid = [0.0, 0.1, 0.5, 1.0]

if !_optin
    _available = false
    @info "FormFactorXrayDB xraydb tests SKIPPED: opt in with SCATTERNET_TEST_XRAYDB=1 " *
          "(that path loads PythonCall and lets CondaPkg provision the numpy + xraydb env)."
else
    # Reach the side environment that holds PythonCall + CondaPkg, then load it.
    # This is the point where the CondaPkg (numpy + xraydb) env gets provisioned.
    let xrenv = joinpath(@__DIR__, "xraydb")
        xrenv in LOAD_PATH || push!(LOAD_PATH, xrenv)
    end
    try
        @eval import PythonCall
    catch e
        @warn "SCATTERNET_TEST_XRAYDB is set but PythonCall could not be loaded from test/xraydb; xraydb tests will be skipped" exception = e
    end
    global _available = try
        compute_form_factors(["fe3+"], 8000.0, [0.1, 0.2]); true
    catch e
        @warn "FormFactorXrayDB backend unavailable; its tests will be skipped" exception = e
        false
    end
    _available || @info "FormFactorXrayDB tests skipped (no numpy/xraydb env)"
end

@testset "FormFactorXrayDB" begin

    # The marker type is pure Julia and testable with or without the Python env.
    @testset "backend marker" begin
        @test FormFactorSourceXrayDB() isa IFACE.FormFactorSource
        @test sprint(showerror, FormFactorError("boom")) == "FormFactorError: boom"
        @test FormFactorError("boom") isa Exception
    end

    if !_available
        @testset "xraydb form-factor tests (opt-in: SCATTERNET_TEST_XRAYDB=1)" begin
            @test_skip "the xraydb extension is actually loaded"
            @test_skip "known fe3+ values at 8000 eV"
            @test_skip "FF container shape / q-dependence / monotonicity"
            @test_skip "lookup ordering, dedup, column selection, off-grid raises"
            @test_skip "compute_form_factors input guards"
        end
    else
        @testset "the xraydb extension is actually loaded" begin
            # Guards against this whole file going vacuous: if PythonCall were not
            # loaded, `compute_form_factors` would still resolve -- to the package's
            # catch-all stub, which throws -- and `_available` would quietly turn
            # every assertion below into a skip.
            @test Base.get_extension(ScatterNet, :FormFactorXrayDBExt) !== nothing
            m = only(methods(compute_form_factors, (Vector{String}, Float64, Vector{Float64})))
            @test parentmodule(m) === Base.get_extension(ScatterNet, :FormFactorXrayDBExt)
            @test _available          # the backend really answered
        end

        @testset "known fe3+ values at 8000 eV" begin
            t = IFACE.form_factor_table(8000.0, ["fe3+"], qvals)
            res = IFACE.form_factors(t, ["fe3+"], qvals)
            @test res isa Matrix{ComplexF64} && size(res) == (1, 2)
            @test check_c(res[1, 1], 21.73320334 + 3.20285267im)
            @test check_c(res[1, 2], 21.71503439 + 3.20285267im)
        end

        @testset "FF container shape" begin
            t = IFACE.form_factor_table(8000.0, ["fe3+", "o2-"], qgrid)
            @test t isa FF
            @test sort(collect(keys(t.tbl))) == ["fe3+", "o2-"]
            @test all(v -> length(v) == length(qgrid), values(t.tbl))
            @test all(v -> v isa Vector{ComplexF64}, values(t.tbl))
            @test t.qmp == Dict(q => i for (i, q) in enumerate(qgrid))  # q => column index
            @test IFACE.form_factor_log(t) isa Vector{String}
            @test IFACE.form_factor_table(8000.0, ["fe3+"], qvals).tbl == compute_form_factors(["fe3+"], 8000.0, qvals).tbl   # form_factor_table is a thin wrapper
        end

        @testset "f0 carries all the q dependence; f1/f2 carry none" begin
            # f(q, E) = f0(s) + f1(E) + i f2(E), so Im(f) must be flat in q...
            t = IFACE.form_factor_table(8000.0, ["fe3+", "o2-"], qgrid)
            for row in values(t.tbl)
                @test all(v -> check_float(imag(v), imag(row[1])), row)
            end
            # ...and changing only E must shift Re(f) by the same amount at every q
            a = IFACE.form_factor_table(8000.0, ["fe3+"], qgrid).tbl["fe3+"]
            b = IFACE.form_factor_table(12000.0, ["fe3+"], qgrid).tbl["fe3+"]
            d = real.(a) .- real.(b)
            @test all(v -> check_float(v, d[1]), d)
            @test !check_float(d[1], 0.0)             # the two energies really do differ
        end

        @testset "f0 decreases monotonically with q" begin
            t = IFACE.form_factor_table(8000.0, ["fe3+", "o2-", "h"], qgrid)
            for row in values(t.tbl)
                re = real.(row)
                @test all(i -> re[i] > re[i + 1], 1:(length(re) - 1))
            end
        end

        @testset "the q -> 0 limit recovers the electron count" begin
            # f0(0) = Z - charge for an ion, Z for a neutral atom; at 1 eV the
            # backend drops to the f0-only tier, so there is no f1/f2 offset.
            t = IFACE.form_factor_table(1.0, ["fe3+", "o2-", "h"], [0.0])
            # atol covers the Cromer-Mann parameterization's own residual at s = 0
            @test isapprox(real(t.tbl["fe3+"][1]), 23.0; atol = 5e-3)   # Fe: Z = 26
            @test isapprox(real(t.tbl["o2-"][1]), 10.0;  atol = 5e-3)   # O:  Z = 8
            @test isapprox(real(t.tbl["h"][1]),    1.0;  atol = 5e-3)
        end

        @testset "f0-only tier is logged and has no imaginary part" begin
            # 1 eV is below the Chantler tabulation range for Fe
            t = IFACE.form_factor_table(1.0, ["fe3+"], qgrid)
            @test any(==("F0-ONLY fe3+"), IFACE.form_factor_log(t))
            @test all(v -> imag(v) == 0.0, t.tbl["fe3+"])
            # at 8000 eV the same ion is full-tier: logged nowhere, f2 non-zero
            full = IFACE.form_factor_table(8000.0, ["fe3+"], qgrid)
            @test isempty(IFACE.form_factor_log(full))
            @test all(v -> imag(v) != 0.0, full.tbl["fe3+"])
        end

        @testset "dummy ion is logged, dropped from the table, and rejected by form_factors" begin
            t = IFACE.form_factor_table(8000.0, ["fe3+", "xx"], qvals)
            @test any(==("DUMMY   xx"), IFACE.form_factor_log(t))
            @test !haskey(t.tbl, "xx")
            # querying an ion the table doesn't hold is a hard error, not a silent drop
            @test_throws FormFactorError IFACE.form_factors(t, ["fe3+", "xx"], qvals)
        end

        @testset "the ion batch is deduped; the query is one row per requested ion" begin
            t = IFACE.form_factor_table(8000.0, ["fe3+", "fe3+", "o2-", "fe3+"], qvals)
            @test length(t.tbl) == 2                      # one row per *unique* ion
            res = IFACE.form_factors(t, ["fe3+", "fe3+"], qvals)
            @test size(res) == (2, length(qvals))        # one row per *requested* ion
            @test res[1, :] == res[2, :]
        end

        @testset "form_factors rows follow the requested ion order" begin
            t = IFACE.form_factor_table(8000.0, ["fe3+", "o2-", "h"], qvals)
            ref = Dict(i => IFACE.form_factors(t, [i], qvals)[1, :] for i in ("fe3+", "o2-", "h"))
            got = IFACE.form_factors(t, ["h", "fe3+", "o2-"], qvals)
            @test got[1, :] == ref["h"] && got[2, :] == ref["fe3+"] && got[3, :] == ref["o2-"]
            got2 = IFACE.form_factors(t, ["o2-", "h"], qvals)
            @test got2[1, :] == ref["o2-"] && got2[2, :] == ref["h"]
        end

        @testset "form_factors selects columns by q, in the requested order" begin
            t = IFACE.form_factor_table(8000.0, ["fe3+"], qgrid)
            row = t.tbl["fe3+"]
            got = IFACE.form_factors(t, ["fe3+"], [0.5, 0.0, 0.5, 1.0])
            @test size(got) == (1, 4)
            @test got[1, :] == row[[3, 1, 3, 4]]
            @test size(IFACE.form_factors(t, ["fe3+"], qgrid)) == (1, length(qgrid))
        end

        @testset "form_factors edge cases: empty ions, empty q" begin
            t = IFACE.form_factor_table(8000.0, ["fe3+"], qvals)
            @test size(IFACE.form_factors(t, String[], qvals)) == (0, length(qvals))
            @test size(IFACE.form_factors(t, ["fe3+"], Float64[])) == (1, 0)
        end

        @testset "form_factors rejects an ion never in the container" begin
            t = IFACE.form_factor_table(8000.0, ["fe3+"], qvals)
            @test_throws FormFactorError IFACE.form_factors(t, ["not_built"], qvals)
            @test_throws FormFactorError IFACE.form_factors(t, ["fe2+"], qvals)   # a different charge state is a miss
        end

        @testset "lookup at an off-grid q raises" begin
            t = IFACE.form_factor_table(8000.0, ["fe3+"], qvals)
            @test_throws FormFactorError IFACE.form_factors(t, ["fe3+"], [0.15])
            @test_throws FormFactorError IFACE.form_factors(t, ["fe3+"], [0.1, 0.15])
            @test_throws FormFactorError IFACE.form_factors(t, ["fe3+"], [nextfloat(0.1)])  # exact match only
            @test_throws FormFactorError IFACE.form_factors(t, String[], [0.15])            # q checked first
        end

        @testset "compute_form_factors input guards raise" begin
            @test_throws FormFactorError compute_form_factors(String[], 8000.0, qvals)
            @test_throws FormFactorError compute_form_factors(["fe3+"], 8000.0, Float64[])
            @test_throws FormFactorError compute_form_factors(["fe3+"], 0.0, qvals)
            @test_throws FormFactorError compute_form_factors(["fe3+"], -1.0, qvals)
            @test_throws FormFactorError compute_form_factors(["fe3+"], 8000.0, [-0.1])
            @test_throws FormFactorError compute_form_factors(["fe3+"], 8000.0, [0.1, -0.1])
        end

        @testset "integer-typed energy and q are accepted" begin
            t = compute_form_factors(["fe3+"], 8000, [0, 1])
            @test sort(collect(keys(t.qmp))) == [0.0, 1.0]
            @test t.tbl["fe3+"] == IFACE.form_factor_table(8000.0, ["fe3+"], [0.0, 1.0]).tbl["fe3+"]
        end

        @testset "log is empty when every ion resolves fully" begin
            @test isempty(IFACE.form_factor_log(IFACE.form_factor_table(8000.0, ["fe3+"], qvals)))
            @test isempty(IFACE.form_factor_log(IFACE.form_factor_table(8000.0, ["fe3+", "o2-", "h"], qvals)))
        end
    end
end
