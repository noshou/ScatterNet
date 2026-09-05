# Exercises src/Molecule/AtomicRadii.jl: ion-string parsing, the ion_key
# round-trip, the raw table lookups, and the resolve_one fallback chain.
using .Molecules: Ion, tryparse_ion, ion_key, ion_radius, element_radius,
                  nearest_ion, resolve_one, _resolve_all, AtomicRadiiSource

lookup_one(ion) = _resolve_all([ion])[1][2]

@testset "AtomicRadii" begin

    @testset "tryparse_ion: magnitude-then-sign and sign-then-magnitude" begin
        @test tryparse_ion("fe3+") == Ion("fe", 3)
        @test tryparse_ion("fe+3") == Ion("fe", 3)
        @test tryparse_ion("fe3-") == Ion("fe", -3)
        @test tryparse_ion("fe-3") == Ion("fe", -3)
        @test tryparse_ion("fe12+") == Ion("fe", 12)     # multi-digit magnitude
        @test tryparse_ion("h1+")  == Ion("h", 1)
        @test tryparse_ion("x+")   == Ion("x", 1)        # 1-letter element
    end

    @testset "tryparse_ion: a bare sign means +/-1" begin
        @test tryparse_ion("fe+") == Ion("fe", 1)
        @test tryparse_ion("fe-") == Ion("fe", -1)
        @test tryparse_ion("fe+") == tryparse_ion("fe1+") == tryparse_ion("fe+1")
        @test tryparse_ion("fe-") == tryparse_ion("fe1-") == tryparse_ion("fe-1")
    end

    @testset "tryparse_ion: surrounding and interior whitespace is tolerated" begin
        @test tryparse_ion(" fe3+ ")  == Ion("fe", 3)
        @test tryparse_ion("fe 3 +")  == Ion("fe", 3)
        @test tryparse_ion("\tfe+3\n") == Ion("fe", 3)
    end

    @testset "tryparse_ion: bare elements are not ions" begin
        # `nothing` here is the signal that sends resolve_one straight to the
        # bare-element table, so it is a contract, not an error.
        @test tryparse_ion("fe") === nothing
        @test tryparse_ion("h")  === nothing
        @test tryparse_ion("rn") === nothing
    end

    @testset "tryparse_ion: charge 0 is rejected outright" begin
        # the magnitude alternation starts at [1-9], so "0" cannot be a magnitude
        @test tryparse_ion("fe0+") === nothing
        @test tryparse_ion("fe+0") === nothing
        @test tryparse_ion("fe0-") === nothing
    end

    @testset "tryparse_ion: unparseable junk" begin
        @test tryparse_ion("")      === nothing
        @test tryparse_ion("Fe3+")  === nothing   # uppercase: table keys are lowercase
        @test tryparse_ion("fe3")   === nothing   # magnitude with no sign
        @test tryparse_ion("abc3+") === nothing   # element is 1-2 letters
        @test tryparse_ion("3+")    === nothing   # no element
        @test tryparse_ion("fe!!")  === nothing
        @test tryparse_ion("fe3++") === nothing
        @test tryparse_ion("fe3.5+") === nothing
    end

    @testset "ion_key: canonical digits-then-sign form, charge 0 rejected" begin
        @test ion_key(Ion("fe", 3))  == "fe3+"
        @test ion_key(Ion("fe", -3)) == "fe3-"
        @test ion_key(Ion("fe", 1))  == "fe1+"
        @test ion_key(Ion("au", -1)) == "au1-"
        @test_throws ArgumentError ion_key(Ion("fe", 0))
    end

    @testset "ion_key normalizes either token order to one key" begin
        @test ion_key(tryparse_ion("fe+3")) == ion_key(tryparse_ion("fe3+")) == "fe3+"
        @test ion_key(tryparse_ion("fe-"))  == "fe1-"
    end

    @testset "raw table lookups" begin
        @test check_float(ion_radius("fe3+"), 0.49)      # pm -> Å conversion applied
        @test check_float(ion_radius("o2-"), 1.35)
        @test ion_radius("qq9+") === nothing
        er = element_radius("fe")
        @test er isa Tuple{Float64,String}
        @test check_float(er[1], 1.274) && er[2] == "metallic"
        @test element_radius("o")[2] == "vdw"
        @test element_radius("qq") === nothing
    end

    @testset "nearest_ion picks the closest on-file charge state" begin
        # fe has 2+, 3+, 4+, 6+ on file (no 5+)
        @test nearest_ion("fe", 3) == "fe3+"          # exact state is also the nearest
        @test nearest_ion("fe", 5) == "fe4+"          # tie 4+/6+ broken toward the first
        @test nearest_ion("fe", 9) == "fe6+"
        @test nearest_ion("fe", -3) == "fe2+"         # distance is on signed charge
        @test nearest_ion("h", -1) == "h1-"           # h has both signs on file
        @test nearest_ion("h", 4) == "h1+"
        @test nearest_ion("qq", 1) === nothing        # element has no charge states
        @test nearest_ion("rn", 1) === nothing        # noble gas: no charge states on file
    end

    @testset "resolve_one, step 1: exact ionic match, both token orderings" begin
        @test check_float(0.49, resolve_one("fe3+"))
        @test check_float(0.49, resolve_one("fe+3"))
        @test check_float(2.2,  resolve_one("au1-"))
        @test check_float(1.35, resolve_one("o2-"))
        @test resolve_one("fe3+") === ion_radius("fe3+")   # straight from the ionic table
    end

    @testset "resolve_one, step 2: nearest-charge fallback" begin
        # fe5+ is not on file; the chain must land on fe4+, not on bare fe
        v5, v4 = resolve_one("fe5+"), resolve_one("fe4+")
        @test v5 !== nothing && v4 !== nothing && check_float(v5, v4)
        @test v5 !== element_radius("fe")[1]              # it did NOT fall to bare fe
        @test check_float(resolve_one("fe5+"), ion_radius(nearest_ion("fe", 5)))
        @test check_float(resolve_one("fe9+"), ion_radius("fe6+"))
    end

    @testset "resolve_one, step 3: bare-element fallback for a parsed ion" begin
        # rn parses as an ion but has no ionic entry and no charge states at all,
        # so the only rung left is atomic_radii
        @test check_float(resolve_one("rn3+"), 2.24)
        @test check_float(resolve_one("rn3+"), element_radius("rn")[1])
        @test check_float(resolve_one("ne1+"), element_radius("ne")[1])
    end

    @testset "resolve_one: unparseable strings go straight to the bare table" begin
        @test check_float(1.274, resolve_one("fe"))
        @test check_float(2.24,  resolve_one("rn"))
        @test check_float(1.1,   resolve_one("h"))
        @test resolve_one("Fe")   === nothing   # unparseable AND not a table key
        @test resolve_one("fe!!") === nothing
    end

    @testset "resolve_one, step 4: total miss returns nothing" begin
        @test resolve_one("qq3+") === nothing   # parses, but no ionic/charge/atomic data
        @test resolve_one("zzzz9+") === nothing # does not even parse
        @test resolve_one("") === nothing
    end

    @testset "resolve_one is exhaustive over the chain's return type" begin
        for s in ("fe3+", "fe5+", "rn3+", "fe", "qq3+", "")
            v = resolve_one(s)
            @test v === nothing || v isa Float64
        end
    end

    @testset "_resolve_all: order/count preserved, repeats deduped" begin
        input = ["fe3+", "zzzz9+", "fe3+", "rn"]
        res = _resolve_all(input)
        @test res isa Vector{Tuple{String,Union{Float64,Nothing}}}
        @test length(res) == length(input)
        @test [k for (k, _) in res] == input
        r1, m, r2, r3 = (v for (_, v) in res)
        @test r1 !== nothing && m === nothing && r3 !== nothing && check_float(r1, r2)
    end

    @testset "_resolve_all: repeated misses are deduped too" begin
        res = _resolve_all(["qq3+", "qq3+", "fe", "qq3+"])
        @test [v for (_, v) in res] == [nothing, nothing, resolve_one("fe"), nothing]
    end

    @testset "_resolve_all agrees with resolve_one entry by entry" begin
        input = ["fe3+", "fe+3", "fe5+", "rn3+", "o2-", "h", "qq3+", "fe!!"]
        @test [v for (_, v) in _resolve_all(input)] == [resolve_one(s) for s in input]
    end

    @testset "_resolve_all: empty input" begin
        res = _resolve_all(String[])
        @test isempty(res)
        @test res isa Vector{Tuple{String,Union{Float64,Nothing}}}
    end

    @testset "_resolve_all: non-String element types are normalized to String" begin
        res = _resolve_all(SubString.(["fe3+", "rn"]))
        @test [k for (k, _) in res] == ["fe3+", "rn"]
        @test all(k -> k isa String, first.(res))
    end

    @testset "AtomicRadiiSource satisfies the Interfaces.lookup contract" begin
        src = AtomicRadiiSource()
        @test src isa ScatterNet.Interfaces.RadiiSource
        ions = ["fe3+", "qq3+", "rn"]
        @test ScatterNet.Interfaces.lookup(src, ions) == _resolve_all(ions)
    end
end
