# Exercises src/Scattering/FormFactorXrayDB.jl. The Python-facing half lives in
# the FormFactorXrayDBExt package extension, which loads with PythonCall
# (`runtests.jl` does `using PythonCall` for exactly that reason). Needs the
# CondaPkg env (numpy + xraydb); the set is skipped, loudly, if unavailable.
using .FormFactorXrayDB: compute_form_factors, FF, FormFactorError,
                         FormFactorSourceXrayDB
const F = FormFactorXrayDB

const _available = try
    compute_form_factors(["fe3+"], 8000.0, [0.1, 0.2]); true
catch e
    @warn "FormFactorXrayDB backend unavailable; its tests will be skipped" exception = e
    false
end

_available || @info "FormFactorXrayDB tests skipped (no numpy/xraydb env)"

check_c(a, b) = abs(a - b) < 1e3 * DEFAULT_ATOL   # looser: against tabulated xraydb reference values, not a closed-form identity
qvals = [0.1, 0.2]
qgrid = [0.0, 0.1, 0.5, 1.0]

@testset "FormFactorXrayDB" begin

    @testset "the xraydb extension is actually loaded" begin
        # Guards against this whole file going vacuous: if PythonCall were not
        # loaded, `compute_form_factors` would still resolve -- to the package's
        # catch-all stub, which throws -- and `_available` would quietly turn
        # every assertion below into a skip.
        @test Base.get_extension(ScatterNet, :FormFactorXrayDBExt) !== nothing
        m = only(methods(compute_form_factors,
                         (Vector{String}, Float64, Vector{Float64})))
        @test parentmodule(m) === Base.get_extension(ScatterNet, :FormFactorXrayDBExt)
        @test _available          # the backend really answered
    end

    # The marker type is pure Julia and testable with or without the Python env.
    @testset "backend marker" begin
        @test FormFactorSourceXrayDB() isa ScatterNet.Interfaces.FormFactorSource
        @test sprint(showerror, FormFactorError("boom")) == "FormFactorError: boom"
        @test FormFactorError("boom") isa Exception
    end

    if !_available
        @test_skip false
    else
        @testset "known fe3+ values at 8000 eV" begin
            t = F.create(8000.0, ["fe3+"], qvals)
            res = F.lookup(t, ["fe3+"], qvals)
            @test length(res) == 1
            ion, arr = res[1]
            @test ion == "fe3+" && length(arr) == 2
            @test check_c(arr[1], 21.73320334 + 3.20285267im)
            @test check_c(arr[2], 21.71503439 + 3.20285267im)
        end

        @testset "FF container shape" begin
            t = F.create(8000.0, ["fe3+", "o2-"], qgrid)
            @test t isa FF
            @test sort(collect(keys(t.tbl))) == ["fe3+", "o2-"]
            @test all(v -> length(v) == length(qgrid), values(t.tbl))
            @test all(v -> v isa Vector{ComplexF64}, values(t.tbl))
            @test t.qmp == Dict(q => i for (i, q) in enumerate(qgrid))  # q => column index
            @test F.log(t) isa Vector{String}
            @test F.create(8000.0, ["fe3+"], qvals).tbl ==
                  compute_form_factors(["fe3+"], 8000.0, qvals).tbl   # create is a thin wrapper
        end

        @testset "f0 carries all the q dependence; f1/f2 carry none" begin
            # f(q, E) = f0(s) + f1(E) + i f2(E), so Im(f) must be flat in q...
            t = F.create(8000.0, ["fe3+", "o2-"], qgrid)
            for row in values(t.tbl)
                @test all(v -> check_float(imag(v), imag(row[1])), row)
            end
            # ...and changing only E must shift Re(f) by the same amount at every q
            a = F.create(8000.0, ["fe3+"], qgrid).tbl["fe3+"]
            b = F.create(12000.0, ["fe3+"], qgrid).tbl["fe3+"]
            d = real.(a) .- real.(b)
            @test all(v -> check_float(v, d[1]), d)
            @test !check_float(d[1], 0.0)             # the two energies really do differ
        end

        @testset "f0 decreases monotonically with q" begin
            t = F.create(8000.0, ["fe3+", "o2-", "h"], qgrid)
            for row in values(t.tbl)
                re = real.(row)
                @test all(i -> re[i] > re[i + 1], 1:(length(re) - 1))
            end
        end

        @testset "the q -> 0 limit recovers the electron count" begin
            # f0(0) = Z - charge for an ion, Z for a neutral atom; at 1 eV the
            # backend drops to the f0-only tier, so there is no f1/f2 offset.
            t = F.create(1.0, ["fe3+", "o2-", "h"], [0.0])
            # atol covers the Cromer-Mann parameterization's own residual at s = 0
            @test isapprox(real(t.tbl["fe3+"][1]), 23.0; atol = 5e-3)   # Fe: Z = 26
            @test isapprox(real(t.tbl["o2-"][1]), 10.0;  atol = 5e-3)   # O:  Z = 8
            @test isapprox(real(t.tbl["h"][1]),    1.0;  atol = 5e-3)
        end

        @testset "f0-only tier is logged and has no imaginary part" begin
            # 1 eV is below the Chantler tabulation range for Fe
            t = F.create(1.0, ["fe3+"], qgrid)
            @test any(==("F0-ONLY fe3+"), F.log(t))
            @test all(v -> imag(v) == 0.0, t.tbl["fe3+"])
            # at 8000 eV the same ion is full-tier: logged nowhere, f2 non-zero
            full = F.create(8000.0, ["fe3+"], qgrid)
            @test isempty(F.log(full))
            @test all(v -> imag(v) != 0.0, full.tbl["fe3+"])
        end

        @testset "dummy ion is logged and dropped from lookup" begin
            t = F.create(8000.0, ["fe3+", "xx"], qvals)
            @test any(==("DUMMY   xx"), F.log(t))
            @test !haskey(t.tbl, "xx")
            @test length(F.lookup(t, ["fe3+", "xx"], qvals)) == 1
            @test first(F.lookup(t, ["fe3+", "xx"], qvals))[1] == "fe3+"
        end

        @testset "the ion batch is deduped; lookup is not" begin
            t = F.create(8000.0, ["fe3+", "fe3+", "o2-", "fe3+"], qvals)
            @test length(t.tbl) == 2                      # one row per *unique* ion
            res = F.lookup(t, ["fe3+", "fe3+"], qvals)
            @test length(res) == 2                        # one result per *queried* ion
            @test res[1] == res[2]
        end

        @testset "lookup preserves the query's ion order" begin
            t = F.create(8000.0, ["fe3+", "o2-", "h"], qvals)
            @test [i for (i, _) in F.lookup(t, ["h", "fe3+", "o2-"], qvals)] ==
                  ["h", "fe3+", "o2-"]
            @test [i for (i, _) in F.lookup(t, ["o2-", "h"], qvals)] == ["o2-", "h"]
        end

        @testset "lookup selects columns by q, in the query's order" begin
            t = F.create(8000.0, ["fe3+"], qgrid)
            row = t.tbl["fe3+"]
            _, got = F.lookup(t, ["fe3+"], [0.5, 0.0, 0.5, 1.0])[1]
            @test got == row[[3, 1, 3, 4]]
            @test length(F.lookup(t, ["fe3+"], qgrid)[1][2]) == length(qgrid)
        end

        @testset "lookup edge cases: empty ions, empty q" begin
            t = F.create(8000.0, ["fe3+"], qvals)
            @test isempty(F.lookup(t, String[], qvals))
            res = F.lookup(t, ["fe3+"], Float64[])
            @test length(res) == 1 && isempty(res[1][2])
        end

        @testset "lookup drops an ion never in the container" begin
            t = F.create(8000.0, ["fe3+"], qvals)
            @test F.lookup(t, ["not_built"], qvals) == []
            @test F.lookup(t, ["fe2+"], qvals) == []      # a different charge state is a miss
        end

        @testset "lookup at an off-grid q raises" begin
            t = F.create(8000.0, ["fe3+"], qvals)
            @test_throws FormFactorError F.lookup(t, ["fe3+"], [0.15])
            @test_throws FormFactorError F.lookup(t, ["fe3+"], [0.1, 0.15])
            @test_throws FormFactorError F.lookup(t, ["fe3+"], [nextfloat(0.1)])  # exact match only
            @test_throws FormFactorError F.lookup(t, String[], [0.15])            # q checked first
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
            @test t.tbl["fe3+"] == F.create(8000.0, ["fe3+"], [0.0, 1.0]).tbl["fe3+"]
        end

        @testset "log is empty when every ion resolves fully" begin
            @test isempty(F.log(F.create(8000.0, ["fe3+"], qvals)))
            @test isempty(F.log(F.create(8000.0, ["fe3+", "o2-", "h"], qvals)))
        end
    end
end
