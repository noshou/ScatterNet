# Exercises the pure-Julia form-factor backend: src/Interfaces/FormFactor/.
# f(q,E) = f0(s) + f1(E) + i*f2(E), s = q/(4pi), from the bundled
# form_factors.sqlite3 (Waasmaier-Kirfel f0, Chantler FFAST anomalous terms).

const IFACE = ScatterNet.Interfaces
using ScatterNet.Interfaces.FormFactor: compute_form_factors, FF, FormFactorError,
                                        FormFactorSourceTables, f0, f1f2, S_MAX

check_c(a, b) = abs(a - b) < 1e3 * DEFAULT_ATOL
qvals = [0.1, 0.2]
qgrid = [0.0, 0.1, 0.5, 1.0]

@testset "FormFactor" begin

    @testset "backend marker" begin
        @test FormFactorSourceTables() isa IFACE.FormFactorSource
        @test sprint(showerror, FormFactorError("boom")) == "FormFactorError: boom"
        @test FormFactorError("boom") isa Exception
    end

    # -----------------------------------------------------------------------
    # Anti-vacuity guard. This is what stops the rest of the file being a set of
    # assertions about nothing: every value below is checked against a reference
    # dumped from the predecessor xraydb/PythonCall implementation, at full
    # precision, over species and energies this file otherwise never touches.
    # Provenance and tolerances: test/fixtures/README.md.
    #
    # It replaces a check that the xraydb package extension was loaded and owned
    # the answering method. That guarded provenance -- *who* answered. With the
    # backend now in-package there is no stub to fall through to, so the useful
    # guard is *what* it answers, pinned to the implementation being replaced.
    # -----------------------------------------------------------------------

    "Rows of a fixture CSV, minus its header."
    _fx(name) = Iterators.drop(eachline(joinpath(@__DIR__, "fixtures", name)), 1)

    "Agreement to within `k` units in the last place at `ref`'s own magnitude."
    _ulp(got, ref, k = 2) = abs(got - ref) <= k * eps(abs(ref))

    @testset "f0 matches the reference to <= 1 ulp over 348 species/s points" begin
        # Exact for all but one point (u6+ at s = 1.989, 1 ulp). The residual is
        # summation order inside `c + Σ a_i exp(...)`, not a different formula --
        # the NumPy exponent association is already reproduced. Chasing the last
        # ulp would pin us to a NumPy implementation detail for no physical gain:
        # independent form-factor tabulations disagree at the 0.4% level.
        n = 0; exact = 0; worst = 0.0
        for ln in _fx("fx_f0.csv")
            ion, s_, ref = split(ln, ',')
            got = f0(String(ion), parse(Float64, s_)); r = parse(Float64, ref)
            n += 1; got === r && (exact += 1)
            @test _ulp(got, r)
            worst = max(worst, abs(got - r) / abs(r))
        end
        @test n == 348                      # the fixture is actually being read
        @test exact >= 347                  # essentially all of it is bit-for-bit
        @test worst < 1e-15
    end

    @testset "f1/f2 match the reference, incl. 500 points across the Fe K edge" begin
        # f2 is exact on ~99% of points; the rest differ by 1 ulp because NumPy's
        # and Julia's `log`/`exp` differ by that much on some arguments. That is a
        # libm difference, not an algorithmic one, and it is not matchable.
        # f1 additionally carries a 7x7 dense solve for the spline coefficients,
        # hence the looser (but still ~5 orders inside the suite's 1e-6) bound.
        n = 0; exact2 = 0; w1 = 0.0; w2 = 0.0
        for ln in _fx("fx_f1f2.csv")
            el, E, r1, r2 = split(ln, ',')
            g1, g2 = f1f2(String(el), parse(Float64, E))
            a = parse(Float64, r1); b = parse(Float64, r2)
            n += 1; g2 === b && (exact2 += 1)
            @test _ulp(g2, b)
            w1 = max(w1, abs(g1 - a)); w2 = max(w2, abs(g2 - b) / abs(b))
        end
        @test n == 852
        @test w1 < 1e-11                    # measured ~6e-14
        @test exact2 / n > 0.98
        @test w2 < 1e-15
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
        @test IFACE.form_factor_table(
            8000.0, 
            ["fe3+"], qvals).tbl == compute_form_factors(["fe3+"], 
            8000.0, 
            qvals
        ).tbl #form_factor_table is a thin wrapper
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

    @testset "f0 is guarded against out-of-range s rather than extrapolating" begin
        # The ionic fits carry large negative constant terms (fe3+: c = -61.93)
        # and go negative well past their fit range, so this is a hard error.
        @test f0("fe3+", S_MAX) isa Float64
        @test_throws FormFactorError f0("fe3+", S_MAX + 1e-9)
        @test_throws FormFactorError f0("fe3+", -1e-9)
        @test_throws FormFactorError f0("not_an_ion", 0.1)
    end

    @testset "an unknown charge state falls back to the neutral atom, and says so" begin
        # The predecessor made this substitution silently; fe4+ scattering as
        # 26 electrons instead of 22 is a real approximation, so it is logged.
        t = IFACE.form_factor_table(8000.0, ["fe4+"], qvals)
        @test any(==("NEUTRAL fe4+"), IFACE.form_factor_log(t))
        @test check_c(real(t.tbl["fe4+"][1]) - real(IFACE.form_factor_table(8000.0, ["fe"], qvals).tbl["fe"][1]), 0.0)
    end

    @testset "f1f2 rejects an element with no Chantler data or an out-of-range energy" begin
        # Chantler covers Z = 1..92; Waasmaier-Kirfel reaches Z = 98, so the
        # actinides past U are f0-only rather than an error at the table level.
        @test_throws FormFactorError f1f2("pu", 8000.0)
        @test_throws FormFactorError f1f2("fe", 0.5)          # below 1.01 eV
        @test_throws FormFactorError f1f2("fe", 1.0e9)
        t = IFACE.form_factor_table(8000.0, ["pu"], qvals)
        @test any(==("F0-ONLY pu"), IFACE.form_factor_log(t))
        @test all(v -> imag(v) == 0.0, t.tbl["pu"])
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
