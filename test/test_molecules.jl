# Exercises src/Molecule/Molecules.jl: construction/centering, the two
# coordinate frames, the lazy radii/vols/r_max accessors, and the error contract.
using .Molecules: Molecule, create, coords_cartesian, coords_spherical,
                  radii, vols, r_max, elms, name, sphere_volume,
                  MoleculeError, _to_tuples
using ScatterNet.Molecule: SASA

# row 1 = r, row 2 = theta, row 3 = phi
r_(m)     = coords_spherical(m)[1, :]
theta_(m) = coords_spherical(m)[2, :]
phi_(m)   = coords_spherical(m)[3, :]

# `done` on the private Lazy fields: the only way to assert that an accessor is
# still unforced without forcing it.
_forced_radii(m) = getfield(getfield(m, :_radii), :done)
_forced_vols(m)  = getfield(getfield(m, :_vols),  :done)
_forced_rmax(m)  = getfield(getfield(m, :_r_max), :done)

# A stand-in RadiiSource, to prove `radii_source` is actually consulted rather
# than the default AtomicRadiiSource being hard-wired in.
struct ConstantRadii <: ScatterNet.Interfaces.RadiiSource
    value::Float64
end
ScatterNet.Interfaces.lookup(s::ConstantRadii, ions::AbstractVector{<:AbstractString}) =
    Tuple{String,Union{Float64,Nothing}}[(String(i), s.value) for i in ions]

struct NeverResolves <: ScatterNet.Interfaces.RadiiSource end
ScatterNet.Interfaces.lookup(::NeverResolves, ions::AbstractVector{<:AbstractString}) =
    Tuple{String,Union{Float64,Nothing}}[(String(i), nothing) for i in ions]

@testset "Molecules" begin

    @testset "two atoms on x axis" begin
        m = create("test", ["h", "h"], [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)])
        @test length(r_(m)) == 2
        @test check_float(r_(m)[1], 1.0) && check_float(r_(m)[2], 1.0)
        @test length(theta_(m)) == 2 && length(phi_(m)) == 2
        @test check_float(theta_(m)[1], π / 2) && check_float(phi_(m)[1], 0.0)
        @test check_float(theta_(m)[2], π / 2) && check_float(phi_(m)[2], π)
    end

    @testset "centering shifts to centroid" begin
        m = create("test", ["h", "h"], [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)])
        @test check_float(r_(m)[1], 1.0) && check_float(r_(m)[2], 1.0)
    end

    @testset "coords_cartesian: shape, centroid, and translation invariance" begin
        pts = [(1.0, 2.0, 3.0), (-4.0, 0.5, 2.0), (7.0, -1.0, 0.0), (0.0, 0.0, 0.0)]
        m = create("test", ["h", "h", "h", "h"], pts)
        c = coords_cartesian(m)
        @test c isa Matrix{Float64}
        @test size(c) == (3, 4)                      # (3, n): rows x/y/z, cols atoms
        for row in 1:3
            @test check_float(sum(c[row, :]) / 4, 0.0)   # centroid is the origin
        end
        # the shape is preserved: only the origin moved
        for i in 1:4, j in 1:4
            d0 = sqrt(sum((pts[i][k] - pts[j][k])^2 for k in 1:3))
            d1 = sqrt(sum((c[k, i] - c[k, j])^2 for k in 1:3))
            @test check_float(d0, d1)
        end
        # translating the input does not change the centered output at all
        shifted = [(p[1] + 10.0, p[2] - 3.0, p[3] + 0.5) for p in pts]
        @test coords_cartesian(create("test", ["h", "h", "h", "h"], shifted)) ≈ c
    end

    @testset "coords_spherical is the polar form of coords_cartesian" begin
        m = create("test", ["h", "h", "h"],
                   [(1.0, 2.0, 3.0), (-2.0, 1.0, -4.0), (0.5, -0.5, 2.0)])
        c, s = coords_cartesian(m), coords_spherical(m)
        @test size(s) == size(c) == (3, 3)
        for j in 1:3
            x, y, z = c[1, j], c[2, j], c[3, j]
            r, θ, φ = s[1, j], s[2, j], s[3, j]
            @test check_float(r, sqrt(x^2 + y^2 + z^2))
            @test 0.0 <= θ <= π
            @test -π <= φ <= π
            # round-trip back to cartesian
            @test check_float(r * sin(θ) * cos(φ), x)
            @test check_float(r * sin(θ) * sin(φ), y)
            @test check_float(r * cos(θ), z)
        end
    end

    @testset "single atom r = 0, theta not NaN" begin
        m = create("test", ["h"], [(5.0, 5.0, 5.0)])
        @test check_float(r_(m)[1], 0.0)
        @test !isnan(theta_(m)[1]) && !isnan(phi_(m)[1])
        # the rsafe = 1.0 clamp makes theta = acos(0) and phi = atan(0, 0)
        @test check_float(theta_(m)[1], π / 2)
        @test check_float(phi_(m)[1], 0.0)
        @test coords_cartesian(m) == zeros(3, 1)
    end

    @testset "an atom sitting exactly on the centroid is r = 0, not NaN" begin
        # 0/0 would also bite an interior atom of a larger molecule
        m = create("test", ["h", "h", "h"],
                   [(-1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
        @test check_float(r_(m)[2], 0.0)
        @test !any(isnan, coords_spherical(m))
    end

    @testset "poles: theta = 0 and theta = pi have no NaN phi" begin
        m = create("test", ["h", "h"], [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0)])
        @test check_float(theta_(m)[1], 0.0)
        @test check_float(theta_(m)[2], π)
        @test !any(isnan, phi_(m))
    end

    @testset "name and elms accessors round-trip the inputs" begin
        m = create("water", ["o", "h", "h"],
                   [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)])
        @test name(m) == "water"
        @test elms(m) == ["o", "h", "h"]
        @test elms(m) isa Vector{String}
        @test name(create(SubString("abc", 1, 2), ["h"], [(0.0, 0.0, 0.0)])) == "ab"
    end

    @testset "vols for known elements" begin
        m = create(
            "test", ["fe", "o", "rn"],
            [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
        )
        @test length(vols(m)) == 3
        ev(x) = (4.0 / 3.0) * π * x^3
        @test check_float(vols(m)[1], ev(1.274))
        @test check_float(vols(m)[3], ev(2.24))
        @test vols(m) == sphere_volume.(radii(m))     # vols is exactly (4/3)πr³ of radii
    end

    @testset "sphere_volume closed form" begin
        @test check_float(sphere_volume(0.0), 0.0)
        @test check_float(sphere_volume(1.0), 4π / 3)
        @test check_float(sphere_volume(2.0), 8 * 4π / 3)   # scales as r³
    end

    @testset "r_max is the largest per-atom radius" begin
        m = create("test", ["fe", "o", "rn"],
                   [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)])
        @test r_max(m) isa Float64
        @test check_float(r_max(m), maximum(radii(m)))
        @test check_float(r_max(m), 2.24)               # rn is the largest of the three
        @test r_max(m) === r_max(m)
        # a single-atom molecule's r_max is just that atom's radius
        @test check_float(r_max(create("t", ["fe"], [(0.0, 0.0, 0.0)])), 1.274)
        # order of the atoms does not matter
        m2 = create("test", ["rn", "fe", "o"],
                    [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)])
        @test check_float(r_max(m), r_max(m2))
    end

    @testset "radii/vols/r_max are lazy and memoized" begin
        m = create("test", ["fe", "o"], [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
        @test !_forced_radii(m) && !_forced_vols(m) && !_forced_rmax(m)
        @test radii(m) === radii(m)          # same object on repeat, not a recompute
        @test _forced_radii(m)
        @test !_forced_vols(m) && !_forced_rmax(m)   # forcing one does not force the others
        @test vols(m) === vols(m)
        @test _forced_vols(m) && !_forced_rmax(m)
        r_max(m)
        @test _forced_rmax(m)
    end

    @testset "vols and r_max force radii as a side effect" begin
        m1 = create("test", ["fe", "o"], [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
        @test !_forced_radii(m1)
        vols(m1)
        @test _forced_radii(m1)                       # vols is defined over force(rad)
        m2 = create("test", ["fe", "o"], [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
        r_max(m2)
        @test _forced_radii(m2) && !_forced_vols(m2)  # r_max needs radii but not vols
    end

    @testset "coordinate frames are computed eagerly, not lazily" begin
        m = create("test", ["fe"], [(1.0, 1.0, 1.0)])
        @test coords_cartesian(m) === coords_cartesian(m)   # a stored field, not a thunk
        @test coords_spherical(m) === coords_spherical(m)
    end

    @testset "unknown element raises only when a radius accessor is forced" begin
        # construction must stay clean: the radii backend is not consulted here
        m = create("test", ["zzzz"], [(0.0, 0.0, 0.0)])
        @test m isa Molecule
        @test size(coords_cartesian(m)) == (3, 1)     # geometry still usable
        @test_throws MoleculeError vols(m)
        @test_throws MoleculeError radii(create("test", ["zzzz"], [(0.0, 0.0, 0.0)]))
        @test_throws MoleculeError r_max(create("test", ["zzzz"], [(0.0, 0.0, 0.0)]))
        # the message names the offending element
        err = try; radii(m); catch e; e; end
        @test occursin("zzzz", sprint(showerror, err))
    end

    @testset "one unresolvable element out of many still raises" begin
        m = create("test", ["fe", "zzzz", "o"],
                   [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)])
        @test_throws MoleculeError radii(m)
    end

    @testset "create computes coords_spherical/vols; repeat access stable" begin
        m = create(
            "test", ["o", "h", "h"],
            [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
        )
        @test length(r_(m)) == 3
        @test check_float(r_(m)[1], 0.4714045208)
        @test length(theta_(m)) == 3 && length(phi_(m)) == 3 && length(vols(m)) == 3
        @test radii(m) === radii(m)          # lazy accessor: memoized, same object on repeat
    end

    @testset "a custom radii_source is honoured" begin
        pts = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
        m = create("test", ["o", "h", "h"], pts; radii_source = ConstantRadii(2.5))
        @test radii(m) == [2.5, 2.5, 2.5]
        @test check_float(r_max(m), 2.5)
        @test all(v -> check_float(v, sphere_volume(2.5)), vols(m))
        # the same elements through the default source give something else entirely
        @test radii(create("test", ["o", "h", "h"], pts)) != radii(m)
        # a source that resolves nothing raises through the same MoleculeError path
        @test_throws MoleculeError radii(create("t", ["o"], [(0.0, 0.0, 0.0)];
                                                radii_source = NeverResolves()))
        # a custom source also lets otherwise-unknown element labels work
        mystery = create("t", ["zzzz"], [(0.0, 0.0, 0.0)]; radii_source = ConstantRadii(1.0))
        @test radii(mystery) == [1.0]
    end

    @testset "_to_tuples accepts any iterable of 3 components" begin
        already = NTuple{3,Float64}[(1.0, 2.0, 3.0)]
        @test _to_tuples(already) === already          # identity fast path, no copy
        @test _to_tuples([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]) ==
              [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0)]
        @test _to_tuples([(1, 2, 3)]) == [(1.0, 2.0, 3.0)]          # Int -> Float64
        @test _to_tuples([(1, 2, 3)]) isa Vector{NTuple{3,Float64}}
        @test _to_tuples([1:3]) == [(1.0, 2.0, 3.0)]                 # a range works too
        @test _to_tuples(((1.0, 2.0, 3.0), (4.0, 5.0, 6.0))) ==
              [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0)]                     # tuple of tuples
    end

    @testset "_to_tuples rejects wrong-length entries" begin
        @test_throws MoleculeError _to_tuples([(1.0, 2.0)])
        @test_throws MoleculeError _to_tuples([(1.0, 2.0, 3.0, 4.0)])
        @test_throws MoleculeError _to_tuples([[1.0]])
        @test_throws MoleculeError _to_tuples([(1.0, 2.0, 3.0), (1.0, 2.0)])  # only one bad
    end

    @testset "create accepts non-tuple coordinate containers" begin
        pts_v = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]
        pts_i = [(1, 0, 0), (-1, 0, 0)]
        ref = coords_cartesian(create("t", ["h", "h"], [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)]))
        @test coords_cartesian(create("t", ["h", "h"], pts_v)) == ref
        @test coords_cartesian(create("t", ["h", "h"], pts_i)) == ref
        @test_throws MoleculeError create("t", ["h"], [(1.0, 0.0)])
    end

    @testset "empty coords raises" begin
        @test_throws MoleculeError create("empty", String[], NTuple{3,Float64}[])
        @test_throws MoleculeError create("empty", String[], Vector{Float64}[])
    end

    @testset "length mismatch raises, in both directions" begin
        @test_throws MoleculeError create("bad", ["o", "h"], [(0.0, 0.0, 0.0)])
        @test_throws MoleculeError create("bad", ["o"],
                                          [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
        @test_throws MoleculeError create("bad", String[], [(0.0, 0.0, 0.0)])
    end

    @testset "MoleculeError prints its message" begin
        @test sprint(showerror, MoleculeError("boom")) == "MoleculeError: boom"
        @test MoleculeError("boom") isa Exception
    end

    @testset "negative ionic radii are clamped to zero" begin
        # Shannon's table stores h1+/c4+/n5+ with negative radii -- extrapolation
        # artifacts, not physical sizes. `_compute_radii` clamps them to 0.0 so
        # volumes stay non-negative and SASA's `r + probe` never inverts.
        for ion in ("h1+", "c4+", "n5+")
            m = create("artifact", [ion], [(0.0, 0.0, 0.0)])
            @test radii(m) == [0.0]
            @test vols(m)  == [0.0]
            @test r_max(m) == 0.0
        end

        # the clamp is floor-only: it must not disturb ordinary positive radii
        m = create("normal", ["fe", "o", "rn"],
                   [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)])
        @test all(>(0.0), radii(m))
        @test all(>(0.0), vols(m))

        # a clamped ion alongside normal atoms leaves the others untouched
        m = create("mixed", ["h1+", "fe"], [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0)])
        @test radii(m)[1] == 0.0
        @test check_float(radii(m)[2], 1.274)
        @test r_max(m) == radii(m)[2]     # the clamped atom never wins r_max
    end

    @testset "a zero-radius atom still has a well-defined SASA" begin
        # clamped to r = 0, the atom is a bare probe-radius sphere rather than
        # an inverted one -- the invariant the clamp exists to protect.
        m = create("proton", ["h1+"], [(0.0, 0.0, 0.0)])
        a = SASA.sasa_atoms(m; n_occ = 64, n_exp = 256, probe = 1.4)
        @test length(a) == 1
        @test check_float(a[1], 4π * 1.4^2)
    end
end
