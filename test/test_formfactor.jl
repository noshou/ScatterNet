# Needs the CondaPkg env (numpy + xraydb); skipped if unavailable.
using .FormFactorXrayDB: compute_form_factors, FormFactorError
const F = FormFactorXrayDB

const _available = try
    compute_form_factors(["fe3+"], 8000.0, [0.1, 0.2]); true
catch; false end

_available || @info "FormFactorXrayDB tests skipped (no numpy/xraydb env)"

check_c(a, b) = abs(a - b) < 1e-6
qvals = [0.1, 0.2]

@testset "FormFactorXrayDB" begin
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

        @testset "dummy ion is logged and dropped from lookup" begin
            t = F.create(8000.0, ["fe3+", "xx"], qvals)
            @test any(==("DUMMY   xx"), F.log(t))
            @test length(F.lookup(t, ["fe3+", "xx"], qvals)) == 1
        end

        @testset "lookup drops an ion never in the container" begin
            t = F.create(8000.0, ["fe3+"], qvals)
            @test F.lookup(t, ["not_built"], qvals) == []
        end

        @testset "lookup at an off-grid q raises" begin
            t = F.create(8000.0, ["fe3+"], qvals)
            @test_throws FormFactorError F.lookup(t, ["fe3+"], [0.15])
        end

        @testset "compute_form_factors input guards raise" begin
            @test_throws FormFactorError compute_form_factors(String[], 8000.0, qvals)
            @test_throws FormFactorError compute_form_factors(["fe3+"], 8000.0, Float64[])
            @test_throws FormFactorError compute_form_factors(["fe3+"], 0.0, qvals)
            @test_throws FormFactorError compute_form_factors(["fe3+"], 8000.0, [-0.1])
        end

        @testset "log is empty when every ion resolves fully" begin
            @test isempty(F.log(F.create(8000.0, ["fe3+"], qvals)))
        end
    end
end
