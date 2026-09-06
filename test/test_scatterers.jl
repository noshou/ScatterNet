# Exercises src/Scattering/Scatterers.jl, the per-species partial-wave terms:
# `_gaussian_dummy` (the Fraser/MacRae/Suzuki dummy amplitude), `excluded`,
# `hydration`, `vacuo` and the `SHELL_THICKNESS` constant. 
# The backbone of the file is a comparison against the Debye double sum
#
#     I(q) = Re Σ_i Σ_j f_i(q) conj(f_j(q)) j_0(q r_ij),   r_ij = |r_i - r_j|
#
# which is the exact orientational average and shares nothing with the
# spherical-harmonic code. The partial-wave sum converges to it as lMax
# grows, so each comparison is run as a convergence ladder (error falling with
# lMax) plus a machine-precision check at the top and a deliberately truncated
# case that must disagree loudly.
#
# Two known gaps:
#   * The Debye sum depends only on pairwise distances, so it is invariant under
#     any rigid motion and cannot see a bug in `_center` or in the cartesian ->
#     spherical conversion. "molecule frame conventions" tests those directly.
#   * `hydration`'s dummy positions come from `SASA.shell_points`, which is
#     pipeline and cannot be meaningfully hardcoded (the cloud depends on the
#     sampling). There the oracle necessarily consumes `shell_points`' output as
#     given, so those tests check only the `shell_points -> B_lm -> I(q)` 
using   .Scattering: _gaussian_dummy, vacuo, excluded, hydration, SHELL_THICKNESS,
        compute_B_lm, partial_wave_weights, self_scatter, cross_scatter
using   .Molecules: create, coords_cartesian, coords_spherical, to_spherical,
        radii, vols, elms, Molecule

using ScatterNet.Molecule: SASA

"""
Per-element van der Waals radii in Å, EXACTLY as the live `AtomicRadii` backend
returns them.

Provenance: dumped at full `Float64` precision from
`Interfaces.lookup(AtomicRadii.AtomicRadiiSource(), [e])` on 2026-09-06, against
`data/atomic_radii.sqlite3` as of commit 071e4fd.
"""
const SCAT_RADII = Dict("c" => 1.77, "o" => 1.5, "h" => 1.2, "fe3+" => 0.49, "o2-" => 1.35)

"Sphere volume from a hardcoded radius, written out rather than taken from `Molecules`."
scat_vol(e) = (4.0 / 3.0) * π * SCAT_RADII[e]^3

"""
Dummy amplitude `f = v exp(-q² v^(2/3) / 4π)`, transcribed from the
Fraser/MacRae/Suzuki form in `_gaussian_dummy`'s docstring rather than called
from it.
"""
scat_gauss(v, q) = v * exp(-q^2 * v^(2 / 3) / (4π))

"Literal momentum-transfer grid in Å⁻¹. q = 0 is present: several closed forms need it."
const SCAT_Q = [0.0, 0.05, 0.15, 0.3, 0.5]

# --- literal fixture geometries -------------------------------------------
# Each is a `(elements, (3, N) cartesian literal)` pair. The molecule under test
# is built from these with `create`, which re-centres them; the oracle reads the
# matrix directly and never asks the molecule for coordinates. That is sound
# because the Debye sum uses only `r_ij`, which no translation can change.

"Two atoms 4 Å apart: the one case whose Debye sum is hand-checkable in closed form."
const SCAT_DIMER_E = ["c", "c"]
const SCAT_DIMER_X = [0.0 4.0; 0.0 0.0; 0.0 0.0]

"""
Cube of 8 identical atoms on ±1.5 Å. Highly symmetric on purpose: whole blocks
of `(l, m)` cancel, so an index/packing bug that a generic molecule would smear
into a small error can show up here as a gross one.
"""
const SCAT_CUBE_E = fill("c", 8)
const SCAT_CUBE_X = reduce(
    hcat, 
    [Float64[1.5x, 1.5y, 1.5z]
    for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
)

"""
20 atoms of three different volumes at irregular positions in a ±4.5 Å box: no
symmetry at all, so the `m != 0` rows carry real weight instead of cancelling.
This is the fixture the convergence ladders lean on. The coordinates were drawn
once from a seeded generator and then frozen as literals, so no test in this
file depends on any RNG state.
"""
const SCAT_BLOB_E = ["o", "h", "c", "o", "h", "c", "o", "h", "c", "o",
                    "h", "c", "o", "h", "c", "o", "h", "c", "o", "h"]
# Each axis is ONE line of 20 numbers. Do not re-wrap: a newline inside a
# matrix literal starts a new row, so splitting these makes it (4, 10).
const SCAT_BLOB_X = [
    -2.660 -0.142  1.295  1.009  4.166  1.157 -4.056  1.445 -1.888 -3.607  2.178 -1.018 -0.333  0.621  3.352  0.320  1.969  0.482 -0.998  4.099
    -0.119 -3.789  3.422 -0.603  2.669  3.933 -3.162  1.902 -0.207  2.393 -3.760 -4.144 -2.267  3.580 -2.385 -2.650  3.236 -3.208 -0.561  0.987
    -1.937  0.693 -0.720  0.261  3.099  0.801 -4.123  1.298  1.323  2.741  2.887 -1.076  4.101  3.973 -2.176 -0.575  1.344 -0.462 -3.859  3.967
]

"""
The 5-atom molecule the opt-in `vacuo` tests use, and its form factors at
8000 eV on `SCAT_Q`, hardcoded.

Provenance: `Interfaces.form_factor_table(8000.0, ["fe3+", "o2-", "h"], SCAT_Q)`
followed by `Interfaces.form_factors`, dumped at full `Float64` precision on
2026-09-06 from the xraydb backend. Only three distinct ions are tabulated; the
per-atom matrix is assembled from them by lookup. `fe3+`'s K edge is at ~7.1
keV, so `f''` is large here and `compute_B_lm` must take its `C = 2` branch --
`f''` is q-independent, which is why each ion's imaginary part is constant
across the row.
"""
const SCAT_FF_E = ["fe3+", "o2-", "h", "o2-", "h"]
const SCAT_FF_X = [ 0.0  2.1  0.0 -1.4  0.9
                    0.0  0.0  1.6  1.9 -2.3
                    0.0  0.0  0.6 -0.7  1.1 ]
const SCAT_FF_ENERGY = 8000.0
const SCAT_FF_TABLE = Dict(
    "fe3+" => [ 21.739265778665303 + 3.2028526697246287im,
                21.737749881791149 + 3.2028526697246287im,
                21.725629596540660 + 3.2028526697246287im,
                21.684813884461967 + 3.2028526697246287im,
                21.588618751767171 + 3.2028526697246287im],
    "o2-"  => [ 10.044358726501459 + 0.03235399864710075im,
                10.041666212623870 + 0.03235399864710075im,
                10.020175140088888 + 0.03235399864710075im,
                9.9482814754245616 + 0.03235399864710075im,
                9.7817211659042105 + 0.03235399864710075im],
    "h"    => [ 0.99944620932546357 + 1.08869064958089e-06im,
                0.99909667497428978 + 1.08869064958089e-06im,
                0.99630692714662494 + 1.08869064958089e-06im,
                0.98697650916313695 + 1.08869064958089e-06im,
                0.96537230148199171 + 1.08869064958089e-06im],
)

# ===========================================================================
# ORACLE
# ===========================================================================

"""
The Debye formula: `I(q) = Re Σ_i Σ_j f_i(q) conj(f_j(q)) j_0(q r_ij)`, the
exact orientational average of `|Σ_i f_i exp(i q·r_i)|²`. Written directly over
cartesian coordinates with `j_0(x) = sin(x)/x`, `j_0(0) = 1`; it borrows nothing
from `SphFuncs`/`PartialWave`. Given two different point sets it is the exact
cross term `Σ_{i∈a} Σ_{j∈b}`, which is what `cross_scatter` converges to.
"""
function scat_debye(qvals, crd_a, f_a, crd_b = crd_a, f_b = f_a)
    out = zeros(Float64, length(qvals))
    for k in eachindex(qvals)
        q = qvals[k]
        acc = zero(ComplexF64)
        for i in axes(crd_a, 2), j in axes(crd_b, 2)
            d = hypot(
                crd_a[1, i] - crd_b[1, j],
                crd_a[2, i] - crd_b[2, j],
                crd_a[3, i] - crd_b[3, j]
            )
            x = q * d
            j0 = iszero(x) ? 1.0 : sin(x) / x   # j_0(0) = 1 by continuity
            acc += f_a[i, k] * conj(f_b[j, k]) * j0
        end
        out[k] = real(acc)
    end
    return out
end

"`(N, Q)` dummy amplitudes for `es` on `qvals`, built from the hardcoded radii alone."
scat_dummy_amp(es, qvals) = [scat_gauss(scat_vol(es[i]), qvals[k]) for i in eachindex(es), k in eachindex(qvals)]

"Largest relative deviation of `got` from `ref`, elementwise."
scat_relerr(got, ref) = maximum(abs.(got .- ref) ./ abs.(ref))

# ===========================================================================
# PIPELINE-SIDE FIXTURES.
# ===========================================================================

scat_mol(name, es, X) = create(name, es, [(X[1, i], X[2, i], X[3, i]) for i in axes(X, 2)])

scat_dimer() = scat_mol("dimer", SCAT_DIMER_E, SCAT_DIMER_X)
scat_cube()  = scat_mol("cube",  SCAT_CUBE_E,  SCAT_CUBE_X)
scat_blob()  = scat_mol("blob",  SCAT_BLOB_E,  SCAT_BLOB_X)

"Water: 3 atoms, entirely convex accessible surface."
scat_water() = create(
    "water", 
    ["o", "h", "h"],
    [(0.0, 0.0, 0.0), (0.9572, 0.0, 0.0), (-0.2400, 0.9266, 0.0)]
)

"Two parallel rows of carbons 6 Å apart: the canyon floor between them is concave."
scat_canyon() = create(
    "canyon", 
    fill("c", 6),
    [(3.0 * x, 3.0 * s, 0.0) for x in -1:1 for s in (-1.0, 1.0)]
)

"Two atoms stacked at the very same point: neither has any accessible surface."
scat_buried() = create("buried", ["c", "c"], [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)])

scat_lmax = 3                         # small: used only by the shape/contract tests
scat_K = (scat_lmax + 1) * (scat_lmax + 2) ÷ 2
scat_chunk = UInt64(2)                # `_CHUNK` is UInt64-typed: `2` would not dispatch
scat_hchunk = UInt64(16)

# lMax at which every ground-truth comparison below is expected to have
# converged. q_max*r_max is at most ~5 for these fixtures (the hydration cloud,
# not the atoms, sets r_max), and the ladders show the error reaching the
# Float64 floor well before 16.
scat_Lconv = 16
scat_wconv = partial_wave_weights(scat_Lconv)
scat_tol = 1e-11                      # measured errors at scat_Lconv are ~1e-13

"""
Section 6's expansion restricted to the two dummy species (no `vacuo`, hence no
Python): with `A = -dns*A_ex + dro*A_sh`,

    I(q) = dns² S_ex,ex + dro² S_sh,sh - 2 dns dro S_ex,sh

This is the assembly the downstream Gram-matrix code will eventually do; it
lives here so the `B_lm` these functions now return can be driven end to end.
"""
scat_toy_I(B_ex, B_sh, w, dns, dro) =
    dns^2 .* self_scatter(B_ex, w) .+
    dro^2 .* self_scatter(B_sh, w) .-
    (2 * dns * dro) .* cross_scatter(B_ex, B_sh, w)

"""
The shell dummies `hydration` would build for `mol`: cartesian positions and
`(M, Q)` amplitudes, obtained by applying the same `classes` filter and
`area * thickness` volume rule the implementation applies.
"""
function scat_shell(mol, qvals; n_target = 80, probe = 1.4,
                    thickness = SHELL_THICKNESS,
                    classes = (SASA.CONVEX, SASA.CONCAVE))
    pts, area, class = SASA.shell_points(mol; probe = probe, n_target = n_target)
    keep = findall(c -> c in classes, class)
    amp = [scat_gauss(area[keep][i] * thickness, qvals[k]) for i in eachindex(keep), k in eachindex(qvals)]
    return pts[:, keep], amp
end

# ===========================================================================
# STUB FORM-FACTOR BACKEND. `vacuo` takes a `form_factor_source`, so the vacuum
# term is drivable without a live xraydb environment. These are test-owned
# types, so defining `Interfaces` methods on them is dispatch, not piracy.
# ===========================================================================

struct ScatStubFF <: Interfaces.FormFactorSource end
struct ScatStubTable
    amp :: Matrix{ComplexF64}   # (n_ions, Q), rows aligned to the ion vector
    ions :: Vector{String}
end

const SCAT_STUB_E = ["a", "b", "c", "b"]
# One row per axis: a newline inside a matrix literal starts a new row, so
# these three lines must stay three lines to keep this (3, 4) and not (1, 12).
const SCAT_STUB_X = [ 0.0  3.2 -1.1  0.4
                      0.0  0.0  2.4 -1.8
                      0.0  0.0  0.9  2.6 ]
"Per-ion `(Re, Im)` weights; `b` is repeated in `SCAT_STUB_E` to pin row lookup."
const SCAT_STUB_W = Dict("a" => 20.0 + 4.0im, "b" => 8.0 + 1.5im, "c" => 1.0 + 0.25im)

scat_stub_f(e, q) = SCAT_STUB_W[e] * exp(-q)

Interfaces.form_factor_table(::ScatStubFF, energy::Real, ions, qvals) =
    ScatStubTable(
        ComplexF64[scat_stub_f(ions[i], qvals[k]) for i in eachindex(ions), k in eachindex(qvals)],
        collect(String, ions)
    )

Interfaces.form_factors(t::ScatStubTable, ions, qvals) = t.amp

"""
Pipeline-side molecule for the stub. The element strings are not real elements,
which is fine here and only here: `create` resolves radii lazily, and `vacuo`
touches only `coords_spherical`, so nothing ever queries the radii backend for
them. Calling `radii`/`vols` on this molecule would throw.
"""
scat_stubmol() = scat_mol("stub", SCAT_STUB_E, SCAT_STUB_X)

"Oracle-side copy of the stub amplitude, built from `SCAT_STUB_W` alone."
const scat_stubamp = ComplexF64[scat_stub_f(SCAT_STUB_E[i], SCAT_Q[k])
                                for i in eachindex(SCAT_STUB_E), k in eachindex(SCAT_Q)]

@testset "Scatterers" begin

    # -----------------------------------------------------------------------
    # Localising guards. These exist so that a change in the radii backend
    # reports itself here rather than as a diffuse Debye mismatch further down.
    # -----------------------------------------------------------------------

    @testset "molecule geometry matches the hardcoded oracle constants" begin
        # If this fails the radii backend has changed: SCAT_RADII (and the
        # volumes derived from it) need regenerating, and every Debye test below
        # will fail for that same single reason, not as separate bugs.
        for (es, X) in (
            (SCAT_DIMER_E, SCAT_DIMER_X),
            (SCAT_CUBE_E, SCAT_CUBE_X),
            (SCAT_BLOB_E, SCAT_BLOB_X)
        )
            mol = scat_mol("guard", es, X)
            @test radii(mol) == [SCAT_RADII[e] for e in es]
            @test all(
                i -> isapprox(vols(mol)[i], 
                scat_vol(es[i]); rtol = 1e-15),
                eachindex(es)
            )
            @test elms(mol) == es
        end
    end

    @testset "molecule frame conventions: centring and the spherical transform" begin
        # The Debye oracle only ever sees pairwise distances, so it is blind to
        # any rigid motion; it cannot catch a bug in `_center` or in
        # `to_spherical`. Both are therefore pinned here directly, against
        # arithmetic done in the test rather than against Molecules' own helpers.
        X = SCAT_BLOB_X
        n = size(X, 2)
        cx = sum(X[1, :]) / n; cy = sum(X[2, :]) / n; cz = sum(X[3, :]) / n
        mol = scat_mol("frame", SCAT_BLOB_E, X)
        C = coords_cartesian(mol)
        @test size(C) == (3, n)
        for i in 1:n
            @test check_float(C[1, i], X[1, i] - cx)
            @test check_float(C[2, i], X[2, i] - cy)
            @test check_float(C[3, i], X[3, i] - cz)
        end
        @test   check_float(sum(C[1, :]), 0.0) && check_float(sum(C[2, :]), 0.0) &&
                check_float(sum(C[3, :]), 0.0)

        # (r, θ, φ) with θ = acos(z/r) from +z and φ = atan(y, x), on the
        # already-centred literals.
        S = coords_spherical(mol)
        @test size(S) == (3, n)
        for i in 1:n
            x = X[1, i] - cx; y = X[2, i] - cy; z = X[3, i] - cz
            r = sqrt(x^2 + y^2 + z^2)
            @test check_float(S[1, i], r)
            @test check_float(S[2, i], acos(z / r))
            @test check_float(S[3, i], atan(y, x))
        end
    end

    @testset "SHELL_THICKNESS is CRYSOL's 3 Å and is hydration's default" begin
        @test SHELL_THICKNESS === 3.0
        mol = scat_water()
        a = hydration(mol, SCAT_Q, scat_lmax, scat_hchunk; n_target = 60)
        b = hydration(mol, SCAT_Q, scat_lmax, scat_hchunk; n_target = 60, thickness = SHELL_THICKNESS)
        @test a == b   # the documented default really is the constant
    end

    # -----------------------------------------------------------------------
    # _gaussian_dummy
    # -----------------------------------------------------------------------

    @testset "_gaussian_dummy shape and element type" begin
        f = _gaussian_dummy([1.0, 8.0, 27.0], SCAT_Q)
        @test f isa Matrix{Float64}
        @test size(f) == (3, length(SCAT_Q))   # (N, Q), compute_B_lm's f_atoms layout
    end

    @testset "_gaussian_dummy matches its closed form entry by entry" begin
        v = [0.5, 4.0, 14.137166941154067, 33.51032163829113]
        f = _gaussian_dummy(v, SCAT_Q)
        for i in eachindex(v), k in eachindex(SCAT_Q)
            @test check_float(f[i, k], scat_gauss(v[i], SCAT_Q[k]))
        end
        # and on the real fixtures, so the oracle's amplitude matrix and the
        # pipeline's are pinned to each other before any Debye comparison
        for es in (SCAT_DIMER_E, SCAT_CUBE_E, SCAT_BLOB_E)
            mol = scat_mol("amp", es, SCAT_BLOB_X[:, 1:length(es)])
            got = _gaussian_dummy(vols(mol), SCAT_Q)
            ref = scat_dummy_amp(es, SCAT_Q)
            @test all(i -> check_float(got[i], ref[i]), eachindex(ref))
        end
    end

    @testset "_gaussian_dummy at q = 0 is the volume itself" begin
        # exp(0) == 1 exactly, so this is an equality, not an approximation:
        # a uniform sphere scatters as its own volume in the forward direction.
        v = [0.5, 4.0, 14.137166941154067]
        @test _gaussian_dummy(v, [0.0])[:, 1] == v
        # the same column appears wherever q = 0 sits in a longer grid
        @test _gaussian_dummy(v, SCAT_Q)[:, 1] == v
    end

    @testset "_gaussian_dummy: a zero-volume dummy is an identically zero row" begin
        f = _gaussian_dummy([0.0, 5.0], SCAT_Q)
        @test all(iszero, f[1, :])          # 0 * exp(...) == 0 at every q
        @test all(>(0.0), f[2, :])          # its neighbour is unaffected
    end

    @testset "_gaussian_dummy is positive and strictly decreasing in q for v > 0" begin
        f = _gaussian_dummy([1.0, 20.0], [0.0, 0.1, 0.4, 1.0, 2.0])
        @test all(>(0.0), f)
        for i in axes(f, 1)
            @test all(k -> f[i, k] > f[i, k + 1], 1:(size(f, 2) - 1))
        end
    end

    @testset "_gaussian_dummy: bigger volumes decay faster in q" begin
        # The Gaussian width is set by v^(1/3), so after dividing out the q = 0
        # value (which is just v) the larger dummy must have fallen further.
        vsmall, vbig = 1.0, 64.0
        q = 0.5
        f = _gaussian_dummy([vsmall, vbig], [0.0, q])
        @test f[2, 2] / f[2, 1] < f[1, 2] / f[1, 1]
        # concretely: v = 64 has v^(2/3) = 16, so its exponent is 16x the v = 1 one
        @test check_float(f[1, 2], exp(-q^2 / (4π)))
        @test check_float(f[2, 2], 64.0 * exp(-16 * q^2 / (4π)))
    end

    @testset "_gaussian_dummy rejects negative volumes" begin
        @test_throws ArgumentError _gaussian_dummy([-1.0], SCAT_Q)
        @test_throws ArgumentError _gaussian_dummy([1.0, 2.0, -1e-12], SCAT_Q)
        @test _gaussian_dummy([0.0], SCAT_Q) isa Matrix{Float64}   # 0 itself is allowed
    end

    @testset "_gaussian_dummy edge cases: empty vols, empty q, single element" begin
        # Broadcasting keeps the (N, Q) contract even when one axis is empty;
        # neither case is an error, so both stay valid `f_atoms` for compute_B_lm.
        @test size(_gaussian_dummy(Float64[], SCAT_Q)) == (0, length(SCAT_Q))
        @test size(_gaussian_dummy([1.0, 2.0], Float64[])) == (2, 0)
        @test size(_gaussian_dummy(Float64[], Float64[])) == (0, 0)
        @test size(_gaussian_dummy([3.0], [0.2])) == (1, 1)
        @test check_float(_gaussian_dummy([3.0], [0.2])[1, 1], scat_gauss(3.0, 0.2))
        # integer volumes are accepted and still come back Float64
        @test _gaussian_dummy([1, 8], SCAT_Q) isa Matrix{Float64}
    end

    # -----------------------------------------------------------------------
    # excluded
    # -----------------------------------------------------------------------

    @testset "excluded returns a bare B_lm array, not a (B_lm, S) tuple" begin
        B = excluded(scat_water(), SCAT_Q, scat_lmax, scat_chunk)
        # This is exactly the contract that changed: the S_ab reduction moved
        # downstream, so nothing is returned here but the multipoles.
        @test !(B isa Tuple)
        @test B isa AbstractArray{<:Complex,3}
        @test eltype(B) <: Complex
        # C == 1 is documented: a dummy sphere's amplitude is real, so there is
        # no anomalous f'' channel to carry.
        @test size(B, 1) == 1
        @test size(B) == (1, scat_K, length(SCAT_Q))
        @test all(isfinite, B)
    end

    @testset "excluded is the PartialWave primitives composed by hand" begin
        mol = scat_water()
        B = excluded(mol, SCAT_Q, scat_lmax, scat_chunk)
        Bhand = compute_B_lm(
            coords_spherical(mol), 
            SCAT_Q,
            _gaussian_dummy(vols(mol), SCAT_Q), 
            scat_lmax, 
            scat_chunk
        )
        @test B == Bhand   # same call, same order: bit-identical
    end

    @testset "excluded's q = 0 closed form is (Σ v_i)/√(4π) in B_00 and 0 above" begin
        # At q = 0, j_l(0) = δ_l0 kills every l > 0, and Y_00 = 1/√(4π), so: B_00(0) = (Σ_i v_i) / √(4π)
        for (es, X) in ((SCAT_DIMER_E, SCAT_DIMER_X), (SCAT_CUBE_E, SCAT_CUBE_X),
                        (SCAT_BLOB_E, SCAT_BLOB_X))
            B = excluded(scat_mol("q0", es, X), [0.0], scat_lmax, scat_chunk)
            tot = sum(scat_vol, es)
            @test isapprox(real(B[1, 1, 1]), tot / sqrt(4π); rtol = 1e-12)
            @test check_float(imag(B[1, 1, 1]), 0.0)
            @test all(k -> check_complex(B[1, k, 1], 0.0), 2:scat_K)
        end
    end

    @testset "excluded is invariant to _CHUNK" begin
        # _CHUNK only sets how many atoms are batched per pass. It changes the
        # summation order, so agreement is tight but not bit-exact.
        mol = scat_blob()
        ref = excluded(mol, SCAT_Q, scat_lmax, UInt64(1))
        for c in (UInt64(2), UInt64(5), UInt64(20), UInt64(1000))
            got = excluded(mol, SCAT_Q, scat_lmax, c)
            @test size(got) == size(ref)
            @test all(i -> check_complex(got[i], ref[i]), eachindex(ref))
        end
    end

    @testset "excluded on a single atom at the origin" begin
        # One atom sits at the centred origin, so r = 0 and j_l(0) = δ_l0 makes
        # every degree above l = 0 vanish at *every* q, not just at q = 0, and
        # B_00(q) is just that atom's amplitude over √(4π).
        mol = create("one", ["c"], [(0.0, 0.0, 0.0)])
        v = scat_vol("c")
        B = excluded(mol, SCAT_Q, scat_lmax, scat_chunk)
        for k in eachindex(SCAT_Q)
            @test isapprox(real(B[1, 1, k]), scat_gauss(v, SCAT_Q[k]) / sqrt(4π); rtol = 1e-12)
            @test all(j -> check_complex(B[1, j, k], 0.0), 2:scat_K)
        end
    end

    # -----------------------------------------------------------------------
    # excluded
    # -----------------------------------------------------------------------

    @testset "excluded reproduces the exact Debye intensity for the dimer" begin
        # Two identical dummies of volume v, 4 Å apart, so the double sum has
        # four terms and collapses by hand to
        #   I(q) = 2 f(q)² (1 + j_0(4q)).
        # Written out here without calling `scat_debye` at all, so the oracle
        # itself is checked against pen and paper before it is trusted below.
        v = scat_vol("c")
        S = self_scatter(excluded(scat_dimer(), SCAT_Q, scat_Lconv, scat_chunk), scat_wconv)
        for k in eachindex(SCAT_Q)
            q = SCAT_Q[k]
            x = 4.0 * q
            j0 = iszero(x) ? 1.0 : sin(x) / x
            @test isapprox(S[k], 2 * scat_gauss(v, q)^2 * (1 + j0); rtol = scat_tol)
        end
        # and the generic oracle agrees with that hand computation
        @test isapprox(scat_debye(
            SCAT_Q, 
            SCAT_DIMER_X,
            scat_dummy_amp(SCAT_DIMER_E, SCAT_Q)), 
            S; 
            rtol = scat_tol
        )
    end

    @testset "excluded converges to Debye as lMax grows" begin
        # The multipole series truncates at lMax with an error set by q*r_max,
        # so this is asserted as a ladder rather than at one lMax: the relative
        # deviation from the exact Debye value has to fall at every rung. Run on
        # the asymmetric blob, the only fixture wide enough that lMax = 2 is
        # genuinely far from converged (the cube's symmetry puts it at the
        # Float64 floor by lMax = 8, which would make a strict ladder meaningless).
        ref = scat_debye(SCAT_Q, SCAT_BLOB_X, scat_dummy_amp(SCAT_BLOB_E, SCAT_Q))
        errs = [scat_relerr(
                self_scatter(
                    excluded(scat_blob(), SCAT_Q, L, UInt64(4)),
                    partial_wave_weights(L)), 
                    ref
                ) for L in (2, 4, 8)]
        @test all(k -> errs[k] > errs[k + 1], 1:(length(errs) - 1))
        @test errs[1] > 1e-3   # the ladder really does start unconverged
    end

    @testset "excluded matches Debye to machine precision at lMax = 16" begin
        for (es, X) in ((SCAT_DIMER_E, SCAT_DIMER_X), (SCAT_CUBE_E, SCAT_CUBE_X),
                        (SCAT_BLOB_E, SCAT_BLOB_X))
            ref = scat_debye(SCAT_Q, X, scat_dummy_amp(es, SCAT_Q))
            got = self_scatter(excluded(scat_mol("dbg", es, X), SCAT_Q,
                                        scat_Lconv, UInt64(4)), scat_wconv)
            @test scat_relerr(got, ref) < scat_tol
        end
    end

    @testset "a deliberately truncated excluded disagrees with Debye loudly" begin
        # Teeth for the tolerance above: at lMax = 1 the series is nowhere near
        # converged for a molecule this wide, so the same comparison must fail
        # by percent-level amounts. If it did not, `scat_tol` would be passing
        # everything rather than pinning anything.
        ref = scat_debye(SCAT_Q, SCAT_BLOB_X, scat_dummy_amp(SCAT_BLOB_E, SCAT_Q))
        bad = self_scatter(
            excluded(scat_blob(), SCAT_Q, 1, UInt64(4)),
            partial_wave_weights(1)
            )
        @test scat_relerr(bad, ref) > 1e-2
        # ... and the disagreement is truncation, not a constant offset: it
        # vanishes at q = 0, where only l = 0 contributes at any lMax.
        @test isapprox(bad[1], ref[1]; rtol = scat_tol)
    end

    # -----------------------------------------------------------------------
    # hydration
    # -----------------------------------------------------------------------

    @testset "hydration returns a bare B_lm array with the documented shape" begin
        B = hydration(scat_water(), SCAT_Q, scat_lmax, scat_hchunk; n_target = 60)
        @test !(B isa Tuple)
        @test B isa AbstractArray{<:Complex,3}
        @test eltype(B) <: Complex
        @test size(B, 1) == 1                       # real dummy amplitude, one channel
        @test size(B) == (1, scat_K, length(SCAT_Q))
        @test all(isfinite, B)
    end

    @testset "hydration argument guards" begin
        mol = scat_water()
        @test_throws ArgumentError hydration(mol, SCAT_Q, scat_lmax, scat_hchunk; thickness = 0.0, n_target = 60)
        @test_throws ArgumentError hydration(mol, SCAT_Q, scat_lmax, scat_hchunk; thickness = -1.0, n_target = 60)
        @test_throws ArgumentError hydration(mol, SCAT_Q, scat_lmax, scat_hchunk; classes = (), n_target = 60)
    end

    @testset "SASA.shell_points is deterministic" begin
        # The plastic sequence carries no RNG, so repeated calls must agree bit-for-bit.
        mol = scat_water()
        a = SASA.shell_points(mol; n_target = 60)
        b = SASA.shell_points(mol; n_target = 60)
        @test a[1] == b[1] && a[2] == b[2] && a[3] == b[3]
        @test   hydration(mol, SCAT_Q, scat_lmax, scat_hchunk; n_target = 60) == 
                hydration(mol, SCAT_Q, scat_lmax, scat_hchunk; n_target = 60)
    end

    @testset "hydration's q = 0 closed form is (total area * thickness)/√(4π)" begin
        # Same l = 0 collapse as `excluded`, with each dummy's volume now
        # `area_i * thickness`; `shell_points` gives every point an equal area,
        # so B_00(0) = (Σ_i area_i * thickness)/√(4π) = total_area*Δ/√(4π).
        mol = scat_water()
        for n in (40, 60, 90)
            area = sum(SASA.shell_points(mol; n_target = n)[2])
            B = hydration(mol, [0.0], scat_lmax, scat_hchunk; n_target = n)
            @test isapprox(real(B[1, 1, 1]), area * SHELL_THICKNESS / sqrt(4π); rtol = 1e-12)
            @test check_float(imag(B[1, 1, 1]), 0.0)
            @test all(k -> check_complex(B[1, k, 1], 0.0), 2:scat_K)
        end
    end

    @testset "hydration's B_00(0) is invariant to n_target" begin
        # Thinning refills `areas = total/length(idx)`, so the total accessible
        # area  is conserved exactly.
        mol = scat_water()
        ref = hydration(mol, [0.0], scat_lmax, scat_hchunk; n_target = 40)[1, 1, 1]
        for n in (60, 90, 120)
            @test isapprox(hydration(mol, [0.0], scat_lmax, scat_hchunk; n_target = n)[1, 1, 1], ref; rtol = 1e-10)
        end
        # n_target really does change the dummy count, so this is not vacuous
        @test length(SASA.shell_points(mol; n_target = 40)[2]) == 40
        @test length(SASA.shell_points(mol; n_target = 90)[2]) == 90
    end

    @testset "hydration scales linearly in thickness at q = 0 only" begin
        # the amplitude scales with `thickness`
        # only where the exponential is exactly 1, i.e. at q = 0. Assert the
        # clean factor there, and only a strictly smaller factor away from it.
        mol = scat_water()
        B1 = hydration(mol, SCAT_Q, scat_lmax, scat_hchunk; n_target = 60, thickness = 1.5)
        B2 = hydration(mol, SCAT_Q, scat_lmax, scat_hchunk; n_target = 60, thickness = 3.0)
        @test isapprox(real(B2[1, 1, 1]), 2 * real(B1[1, 1, 1]); rtol = 1e-10)
        for k in 2:length(SCAT_Q)
            r = real(B2[1, 1, k]) / real(B1[1, 1, k])
            @test 1.0 < r < 2.0    # thicker dummies are also wider, so they decay sooner
        end
    end

    @testset "hydration's classes filter selects real subsets" begin
        # Water's accessible surface is entirely CONVEX, so restricting to
        # CONVEX must be a no-op there ...
        wat = scat_water()
        @test all(==(SASA.CONVEX), SASA.shell_points(wat; n_target = 60)[3])
        @test   hydration(wat, SCAT_Q, scat_lmax, scat_hchunk;
                        n_target = 60, classes = (SASA.CONVEX, SASA.CONCAVE)) ==
                hydration(wat, SCAT_Q, scat_lmax, scat_hchunk;
                        n_target = 60, classes = (SASA.CONVEX,))

        # ... while the canyon really does have recessed CONCAVE beads, so
        # dropping them has to change the answer.
        can = scat_canyon()
        cls = SASA.shell_points(can; n_target = 120)[3]
        @test count(==(SASA.CONCAVE), cls) > 0     # the fixture is doing its job
        both = hydration(   can, [0.0], scat_lmax, scat_hchunk;
                            n_target = 120, classes = (SASA.CONVEX, SASA.CONCAVE))
        conv = hydration(   can, [0.0], scat_lmax, scat_hchunk;
                            n_target = 120, classes = (SASA.CONVEX,))
        conc = hydration(   can, [0.0], scat_lmax, scat_hchunk;
                            n_target = 120, classes = (SASA.CONCAVE,))
        @test both != conv
        # dropping beads drops area, which B_00(0) sees directly ...
        @test real(conv[1, 1, 1]) < real(both[1, 1, 1])
        @test real(conc[1, 1, 1]) > 0.0
        # ... and the two disjoint subsets partition the whole surface, so their
        # forward multipoles add back up to it exactly.
        @test isapprox( real(conv[1, 1, 1]) + real(conc[1, 1, 1]),
                        real(both[1, 1, 1]); rtol = 1e-12)
    end

    @testset "hydration on a molecule with no accessible surface is all zeros" begin
        # Two atoms stacked at one point occlude each other completely; the
        # docstring promises an all-zero B_lm rather than an error.
        mol = scat_buried()
        @test size(SASA.shell_points(mol)[1], 2) == 0
        B = hydration(mol, SCAT_Q, scat_lmax, scat_hchunk)
        @test size(B) == (1, scat_K, length(SCAT_Q))
        @test all(iszero, B)
    end

    @testset "hydration forwards probe to shell_points" begin
        mol = scat_water()
        # The default probe is water's 1.4 Å ...
        @test   hydration(mol, [0.0], scat_lmax, scat_hchunk; n_target = 60) ==
                hydration(mol, [0.0], scat_lmax, scat_hchunk; n_target = 60, probe = 1.4)
        # ... and a bigger probe inflates the accessible surface, which the
        # forward multipole sees exactly via total area * thickness.
        big = hydration(mol, [0.0], scat_lmax, scat_hchunk; n_target = 60, probe = 2.5)
        small = hydration(mol, [0.0], scat_lmax, scat_hchunk; n_target = 60, probe = 1.4)
        @test real(big[1, 1, 1]) > real(small[1, 1, 1])
        area = sum(SASA.shell_points(mol; probe = 2.5, n_target = 60)[2])
        @test isapprox(real(big[1, 1, 1]), area * SHELL_THICKNESS / sqrt(4π); rtol = 1e-12)
    end

    @testset "hydration is invariant to _CHUNK" begin
        # As for `excluded`: chunking only reorders the summation.
        mol = scat_water()
        ref = hydration(mol, SCAT_Q, scat_lmax, UInt64(1); n_target = 40)
        for c in (UInt64(7), UInt64(40), UInt64(10_000))
            got = hydration(mol, SCAT_Q, scat_lmax, c; n_target = 40)
            @test all(i -> check_complex(got[i], ref[i]), eachindex(ref))
        end
    end

    # -----------------------------------------------------------------------
    # hydration
    # -----------------------------------------------------------------------

    @testset "hydration converges to the Debye intensity of its shell cloud" begin
        # The shell cloud sits a probe radius plus an atomic radius outside the
        # atoms, so its r_max is roughly 3-4 Å larger than the molecule's and it
        # needs a correspondingly higher lMax for the same q, which is why the
        # ladder here starts at 4 rather than 2, and why `scat_Lconv = 16` is
        # sized off the shell rather than off the atoms.
        mol = scat_blob()
        P, fh = scat_shell(mol, SCAT_Q; n_target = 80)
        ref = scat_debye(SCAT_Q, P, fh)
        errs = [scat_relerr(
                    self_scatter(
                        hydration(mol, SCAT_Q, L, UInt64(64); n_target = 80),
                        partial_wave_weights(L)
                    ), 
                    ref
                ) for L in (4, 8)]
        @test errs[1] > errs[2]
        @test errs[1] > 1e-3
        @test scat_relerr(self_scatter(
            hydration(mol, SCAT_Q, scat_Lconv, UInt64(64); n_target = 80), 
            scat_wconv), ref) < scat_tol
    end

    @testset "hydration matches Debye on the cube's shell cloud" begin
        mol = scat_cube()
        P, fh = scat_shell(mol, SCAT_Q; n_target = 80)
        got = self_scatter(hydration(mol, SCAT_Q, scat_Lconv, UInt64(64); n_target = 80), scat_wconv)
        @test scat_relerr(got, scat_debye(SCAT_Q, P, fh)) < scat_tol
    end

    @testset "hydration matches Debye under a non-default probe and thickness" begin
        # The keywords have to reach `shell_points`/`_gaussian_dummy` intact for
        # the reference (built from the same keywords by hand) to match.
        mol = scat_blob()
        P, fh = scat_shell(mol, SCAT_Q; n_target = 70, probe = 2.0, thickness = 5.0)
        got = self_scatter(hydration(mol, SCAT_Q, scat_Lconv, UInt64(64);
                                    n_target = 70, probe = 2.0, thickness = 5.0), scat_wconv)
        @test scat_relerr(got, scat_debye(SCAT_Q, P, fh)) < scat_tol
    end

    @testset "hydration matches Debye over a CONVEX-only sub-cloud" begin
        # Filtering by class must drop exactly the beads the reference drops --
        # a mismatched filter would change r_max and the amplitudes together and
        # would not cancel out of the comparison.
        mol = scat_canyon()
        cls = (SASA.CONVEX,)
        P, fh = scat_shell(mol, SCAT_Q; n_target = 120, classes = cls)
        got = self_scatter(hydration(mol, SCAT_Q, scat_Lconv, UInt64(64);
                                    n_target = 120, classes = cls), scat_wconv)
        @test scat_relerr(got, scat_debye(SCAT_Q, P, fh)) < scat_tol
    end

    # -----------------------------------------------------------------------
    # The cross term
    # -----------------------------------------------------------------------

    @testset "cross_scatter converges to the exact Debye cross sum" begin
        # Σ_{i in ex} Σ_{j in sh} f_i conj(f_j) j_0(q r_ij) between two DIFFERENT
        # point sets. This is the comparison that catches a sign or conjugation
        # slip in the cross path.
        #
        # NOTE the coordinate frames: the atom side uses the raw literals while
        # the shell side is in the molecule's centred frame. Both point sets have
        # to share one origin for a cross term, so the atoms are re-centred here
        # by hand (again from the literals, not from `coords_cartesian`).
        mol = scat_blob()
        n = size(SCAT_BLOB_X, 2)
        Xc = SCAT_BLOB_X .- (sum(SCAT_BLOB_X; dims = 2) ./ n)
        fe = scat_dummy_amp(SCAT_BLOB_E, SCAT_Q)
        P, fh = scat_shell(mol, SCAT_Q; n_target = 80)
        ref = scat_debye(SCAT_Q, Xc, fe, P, fh)
        errs = [scat_relerr(cross_scatter(  excluded(mol, SCAT_Q, L, UInt64(4)),
                                            hydration(mol, SCAT_Q, L, UInt64(64);
                                                    n_target = 80),
                                            partial_wave_weights(L)), ref)
                for L in (4, 8)]
        @test errs[1] > errs[2]
        @test errs[1] > 1e-3
        B_ex = excluded(mol, SCAT_Q, scat_Lconv, UInt64(4))
        B_sh = hydration(mol, SCAT_Q, scat_Lconv, UInt64(64); n_target = 80)
        @test scat_relerr(cross_scatter(B_ex, B_sh, scat_wconv), ref) < scat_tol
        # the cross sum is symmetric under swapping the two species
        @test cross_scatter(B_ex, B_sh, scat_wconv) == cross_scatter(B_sh, B_ex, scat_wconv)
    end

    # -----------------------------------------------------------------------
    # The assembled I(q)
    # -----------------------------------------------------------------------

    @testset "assembled I(q) matches the Debye-assembled intensity" begin
        # Section 6 restricted to the two dummy species. Both sides use the same
        # (dns, dro); the left side goes through B_lm and the weighted inner
        # product, the right side through three independent O(N²) real-space
        # sums. Agreement closes the loop from molecule geometry to intensity.
        for (es, X, mkmol) in ((SCAT_CUBE_E, SCAT_CUBE_X, scat_cube),
                                (SCAT_BLOB_E, SCAT_BLOB_X, scat_blob))
            mol = mkmol()
            Xc = X .- (sum(X; dims = 2) ./ size(X, 2))   # shared origin, see above
            fe = scat_dummy_amp(es, SCAT_Q)
            P, fh = scat_shell(mol, SCAT_Q; n_target = 80)
            D_ex = scat_debye(SCAT_Q, Xc, fe)
            D_sh = scat_debye(SCAT_Q, P, fh)
            D_x  = scat_debye(SCAT_Q, Xc, fe, P, fh)
            B_ex = excluded(mol, SCAT_Q, scat_Lconv, UInt64(4))
            B_sh = hydration(mol, SCAT_Q, scat_Lconv, UInt64(64); n_target = 80)
            for (dns, dro) in ((0.334, 0.03), (1.0, 1.0), (0.2, 0.9))
                ref = dns^2 .* D_ex .+ dro^2 .* D_sh .- (2 * dns * dro) .* D_x
                @test scat_relerr(scat_toy_I(B_ex, B_sh, scat_wconv, dns, dro), ref) < 1e-9
            end
        end
    end

    @testset "assembled I(0) is (dns*ΣV_ex - dro*ΣV_sh)², straight from geometry" begin
        # At q = 0 every j_0 is 1 and the Debye sums collapse to products of
        # totals, so the whole model becomes a perfect square:
        #   I(0) = (dns Σ_i v_i - dro * total_area * thickness)².
        # A pen-and-paper endpoint that needs neither the oracle nor B_lm.
        for (es, X, mkmol) in ((SCAT_CUBE_E, SCAT_CUBE_X, scat_cube),
                                (SCAT_BLOB_E, SCAT_BLOB_X, scat_blob))
            mol = mkmol()
            n = 80
            B_ex = excluded(mol, [0.0], scat_Lconv, UInt64(4))
            B_sh = hydration(mol, [0.0], scat_Lconv, UInt64(64); n_target = n)
            V_ex = sum(scat_vol, es)
            V_sh = sum(SASA.shell_points(mol; n_target = n)[2]) * SHELL_THICKNESS
            @test isapprox(self_scatter(B_ex, scat_wconv)[1], V_ex^2; rtol = 1e-12)
            @test isapprox(self_scatter(B_sh, scat_wconv)[1], V_sh^2; rtol = 1e-12)
            @test isapprox(cross_scatter(B_ex, B_sh, scat_wconv)[1], V_ex * V_sh; rtol = 1e-12)
            for (dns, dro) in ((0.334, 0.03), (1.0, 0.0), (0.2, 0.9))
                @test isapprox(scat_toy_I(B_ex, B_sh, scat_wconv, dns, dro)[1],
                               (dns * V_ex - dro * V_sh)^2; rtol = 1e-10)
            end
        end
    end

    @testset "assembled I(q) sanity guards: non-negative, Cauchy-Schwarz, exact limits" begin
        # Cheap structural guards, kept only to catch gross damage; the
        # ground-truth comparisons above are what actually pin the values.
        # I is, mode by mode, |dns*B_ex - dro*B_sh|², so non-negativity follows
        # from the weighted inner product being positive semi-definite.
        mol = scat_blob()
        B_ex = excluded(mol, SCAT_Q, scat_Lconv, UInt64(4))
        B_sh = hydration(mol, SCAT_Q, scat_Lconv, UInt64(64); n_target = 80)
        S_ex = self_scatter(B_ex, scat_wconv)
        S_sh = self_scatter(B_sh, scat_wconv)
        X = cross_scatter(B_ex, B_sh, scat_wconv)
        for k in eachindex(SCAT_Q)
            @test abs(X[k]) <= sqrt(S_ex[k] * S_sh[k]) * (1 + 1e-10)
        end
        for (dns, dro) in ((0.334, 0.03), (0.5, -0.2), (-1.7, 2.3))
            I = scat_toy_I(B_ex, B_sh, scat_wconv, dns, dro)
            @test I isa AbstractVector{<:Real}
            @test length(I) == length(SCAT_Q)
            @test all(isfinite, I) && all(>=(0.0), I)
        end
        # switching a contrast off drops two terms exactly (0.0 * x == 0.0)
        @test scat_toy_I(B_ex, B_sh, scat_wconv, 0.7, 0.0) == 0.7^2 .* S_ex
        @test scat_toy_I(B_ex, B_sh, scat_wconv, 0.0, 0.4) == 0.4^2 .* S_sh
        @test all(iszero, scat_toy_I(B_ex, B_sh, scat_wconv, 0.0, 0.0))
        # and with the two species identical the three terms cancel identically
        @test cross_scatter(B_ex, B_ex, scat_wconv) == S_ex
        @test all(iszero, scat_toy_I(B_ex, B_ex, scat_wconv, 0.75, 0.75))
    end

    # -----------------------------------------------------------------------
    # vacuo
    # -----------------------------------------------------------------------

    let
        scat_ffmol = scat_mol("ffmol", SCAT_FF_E, SCAT_FF_X)
        # oracle-side per-atom amplitude, assembled from the hardcoded table only
        scat_ffamp = reduce(vcat, [permutedims(SCAT_FF_TABLE[e]) for e in SCAT_FF_E])

        @testset "the live form factors match the hardcoded oracle table" begin
            # Same role as the radii guard: if the xraydb tables or the backend
            # change, that reports itself here rather than as a Debye mismatch,
            # and SCAT_FF_TABLE is what needs regenerating.
            tbl = Interfaces.form_factor_table(SCAT_FF_ENERGY, SCAT_FF_E, SCAT_Q)
            amp = Interfaces.form_factors(tbl, SCAT_FF_E, SCAT_Q)
            @test size(amp) == size(scat_ffamp)
            @test all(i -> check_complex(amp[i], scat_ffamp[i]), eachindex(scat_ffamp))
            @test any(x -> imag(x) != 0, scat_ffamp)   # the fixture really is anomalous
        end

        @testset "vacuo returns a bare B_lm array with two channels near an edge" begin
            B = vacuo(scat_ffmol, SCAT_Q, scat_lmax, SCAT_FF_E, SCAT_FF_ENERGY, scat_chunk)
            @test !(B isa Tuple)
            @test B isa AbstractArray{<:Complex,3}
            @test size(B) == (2, scat_K, length(SCAT_Q))   # f'' != 0 at 8 keV
            @test all(isfinite, B)
            # ... and a real amplitude collapses back to one channel: at 1 eV
            # the tabulated f'' is zero (test_formfactor.jl pins that too).
            @test size(vacuo(scat_ffmol, SCAT_Q, scat_lmax, SCAT_FF_E, 1.0, scat_chunk), 1) == 1
        end

        @testset "vacuo is compute_B_lm composed by hand from the facade" begin
            tbl = Interfaces.form_factor_table(SCAT_FF_ENERGY, SCAT_FF_E, SCAT_Q)
            amp = Interfaces.form_factors(tbl, SCAT_FF_E, SCAT_Q)
            @test vacuo(scat_ffmol, SCAT_Q, scat_lmax, SCAT_FF_E,
                        SCAT_FF_ENERGY, scat_chunk) ==
                        compute_B_lm(coords_spherical(scat_ffmol), SCAT_Q, amp,
                            scat_lmax, scat_chunk)
        end

        @testset "vacuo's q = 0 limit is Σ_i f_i(0)/√(4π), channel by channel" begin
            # Same l = 0 collapse as the dummy species, but Re(f) and Im(f) land
            # in separate channels rather than being folded together. Both sums
            # come from the hardcoded table.
            B = vacuo(scat_ffmol, [0.0], scat_lmax, SCAT_FF_E, SCAT_FF_ENERGY, scat_chunk)
            @test isapprox(real(B[1, 1, 1]),
                            sum(real, scat_ffamp[:, 1]) / sqrt(4π); rtol = 1e-10)
            @test isapprox(real(B[2, 1, 1]),
                            sum(imag, scat_ffamp[:, 1]) / sqrt(4π); rtol = 1e-10)
            for c in axes(B, 1)
                @test all(k -> check_complex(B[c, k, 1], 0.0), 2:scat_K)
            end
        end

        @testset "vacuo reproduces the exact complex-f Debye intensity" begin
            # The general complex Debye form, Re Σ_ij f_i conj(f_j) j_0(q r_ij),
            # against the two-channel partial-wave sum. `partial_wave_weights`' 
            # docstring records a MEASURED 3% error from folding a complex f directly 
            # instead of splitting it, so a mismatch here that does not shrink with lMax
            # would mean the channel split is wrong, not that lMax is too small.
            ref = scat_debye(SCAT_Q, SCAT_FF_X, scat_ffamp)
            errs = [scat_relerr(self_scatter(vacuo(scat_ffmol, SCAT_Q, L, SCAT_FF_E,
                                                    SCAT_FF_ENERGY, UInt64(4)),
                                                partial_wave_weights(L)), ref)
                    for L in (1, 2, 4)]
            @test all(k -> errs[k] > errs[k + 1], 1:(length(errs) - 1))
            @test errs[1] > 1e-3
            @test scat_relerr(self_scatter(vacuo(scat_ffmol, SCAT_Q, scat_Lconv,
                                                SCAT_FF_E, SCAT_FF_ENERGY, UInt64(4)),
                                            scat_wconv), ref) < scat_tol
        end
    end

    # -----------------------------------------------------------------------
    # vacuo
    # -----------------------------------------------------------------------

    # `vacuo`'s `form_factor_source` keyword exists so the vacuum term can be
    # driven without a live xraydb environment. Nothing else in the suite uses
    # it, so without these tests the seam could stop dispatching and every
    # xraydb-backed assertion above would still pass. The stub is defined on
    # test-owned types, so adding `Interfaces` methods for it is not piracy.
    @testset "vacuo drives an injected form-factor backend" begin
        got = vacuo(scat_stubmol(), SCAT_Q, scat_Lconv, SCAT_STUB_E, 1234.0,
                    scat_chunk; form_factor_source = ScatStubFF())

        @testset "the stub really is what answered" begin
            # If the keyword were ignored and the default xraydb backend ran
            # instead, the amplitudes would be Fe/O/H form factors rather than
            # the stub's, and the Debye check below would fail.
            @test Interfaces.form_factor_table(ScatStubFF(), 1234.0, SCAT_STUB_E,
                                                SCAT_Q) isa ScatStubTable
            @test size(got, 1) == 2          # Im(f) != 0, so the C = 2 branch
            @test size(got) == (2, length(scat_wconv), length(SCAT_Q))
        end

        @testset "an injected backend needs no xraydb call" begin
            # The stub path must not reach the extension at all. `energy` is
            # 1234.0 eV, which is below Fe's tabulated range: had the real
            # backend answered, it would have logged an F0-ONLY downgrade and
            # returned a purely real f, collapsing this to one channel.
            @test size(got, 1) == 2
        end

        @testset "vacuo reproduces the exact Debye intensity through the stub" begin
            # Same ground-truth comparison as the xraydb fixture, but with an
            # amplitude this file defines outright, so the oracle shares nothing
            # with the backend under test.
            ref = scat_debye(SCAT_Q, SCAT_STUB_X, scat_stubamp)
            @test scat_relerr(self_scatter(got, scat_wconv), ref) < scat_tol
        end

        @testset "the channel split survives an f'' as large as f'" begin
            # `partial_wave_weights`' docstring records a MEASURED 3% error from
            # folding a complex f directly instead of splitting it. Real 8 keV
            # f''/f' is ~0.15; here it is 1.0, which magnifies that failure mode
            # by ~an order of magnitude. Converging to the Float64 floor anyway
            # is strong evidence the two channels are combined correctly.
            heavy = ComplexF64[(1.0 + 1.0im) * scat_stubamp[i, k]
                                for i in axes(scat_stubamp, 1), k in axes(scat_stubamp, 2)]
            B = compute_B_lm(coords_spherical(scat_stubmol()), SCAT_Q, heavy,
                                scat_Lconv, scat_chunk)
            @test size(B, 1) == 2
            @test scat_relerr(self_scatter(B, scat_wconv),
                                scat_debye(SCAT_Q, SCAT_STUB_X, heavy)) < scat_tol
        end

        @testset "a real-amplitude stub collapses to one channel" begin
            # The same seam with Im(f) identically zero must take the C = 1
            # branch, and still match Debye.
            re = ComplexF64.(real.(scat_stubamp))
            B = compute_B_lm(coords_spherical(scat_stubmol()), SCAT_Q, re,
                            scat_Lconv, scat_chunk)
            @test size(B, 1) == 1
            @test scat_relerr(self_scatter(B, scat_wconv),
                                scat_debye(SCAT_Q, SCAT_STUB_X, re)) < scat_tol
        end
    end
end
