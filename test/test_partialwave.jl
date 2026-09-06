# Exercises src/Scattering/PartialWave.jl at the unit level: the packed (l, m)
# ordering shared by `partial_wave_weights` and `compute_B_lm`, the `_deg_contrib`
# contraction, the real/imaginary two-channel split, chunking invariance,
# linearity/additivity of B_lm, and the algebraic identities relating
# `self_scatter` and `cross_scatter`. 
using .Scattering: compute_B_lm, partial_wave_weights, self_scatter, cross_scatter, _deg_contrib
using .SphFuncs: sphHarm, sphBess
using .Molecules: create, coords_spherical, to_spherical

# packed row offsets, duplicated here by hand so the tests pin the convention
# rather than inheriting it from the code under test
pw_idx(l, m) = l * (l + 1) ÷ 2 + m + 1
pw_nrows(lMax) = (lMax + 1) * (lMax + 2) ÷ 2

# Straight transcription of the definition in the module docstring:
#   B_lm(q) = Σ_i f_i(q) * j_l(q r_i) * conj(Y_lm(θ_i, φ_i))
function pw_naive_B(coords_sph, qvals, f_atoms, lMax)
    N = size(coords_sph, 2)
    Q = length(qvals)
    B = zeros(ComplexF64, pw_nrows(lMax), Q)
    for i in 1:N
        Y = sphHarm(lMax, [coords_sph[2, i]], [coords_sph[3, i]])   # (K, 1)
        j = sphBess([coords_sph[1, i]], qvals, lMax)                # (lMax+1, Q, 1)
        for l in 0:lMax, m in 0:l
            k = pw_idx(l, m)
            for qi in 1:Q
                B[k, qi] += f_atoms[i, qi] * j[l + 1, qi, 1] * conj(Y[k, 1])
            end
        end
    end
    return B
end

# A fixed 5-atom geometry plus a fixed q-grid and amplitude matrix.
pw_cart = [ 0.0  1.3 -2.1  0.4  1.7;
            0.0 -0.7  1.1  2.2 -1.9;
            1.5  0.9  0.3 -1.4  0.6]
pw_sph = to_spherical(pw_cart)
pw_q = [0.0, 0.15, 0.4, 0.9]
# amplitudes decaying with q, distinct per atom, all real
pw_f = [(1.0 + 0.5i) * exp(-0.3 * q^2) for i in 1:5, q in pw_q]

@testset "PartialWave" begin

    @testset "partial_wave_weights: exact values and packed layout" begin
        # (0,0), (1,0), (1,1), (2,0), (2,1), (2,2), m = 0..l
        @test partial_wave_weights(0) == [1.0]
        @test partial_wave_weights(1) == [1.0, 1.0, 2.0]
        @test partial_wave_weights(2) == [1.0, 1.0, 2.0, 1.0, 2.0, 2.0]
        @test partial_wave_weights(3) == [1.0, 1.0, 2.0, 1.0, 2.0, 2.0, 1.0, 2.0, 2.0, 2.0]
        @test partial_wave_weights(2) isa Vector{Float64}
    end

    @testset "partial_wave_weights: the 1.0 entries land exactly on m = 0" begin
        # this is the load-bearing indexing claim; the weight vector is
        # consumed positionally against B_lm's packed second axis
        for lMax in 0:6
            w = partial_wave_weights(lMax)
            @test length(w) == pw_nrows(lMax)
            for l in 0:lMax, m in 0:l
                @test w[pw_idx(l, m)] == (m == 0 ? 1.0 : 2.0)
            end
            # equivalently: exactly lMax+1 ones, the rest twos
            @test count(==(1.0), w) == lMax + 1
            @test count(==(2.0), w) == length(w) - (lMax + 1)
        end
        # a larger lMax only appends; the prefix is untouched
        @test partial_wave_weights(6)[1:pw_nrows(3)] == partial_wave_weights(3)
    end

    @testset "partial_wave_weights exception contract" begin
        @test_throws ArgumentError partial_wave_weights(-1)
        @test_throws ArgumentError partial_wave_weights(-4)
        # lMax = 0 is the boundary of the valid range, not an error
        @test partial_wave_weights(0) == [1.0]
    end

    @testset "_deg_contrib equals the explicit triple loop" begin
        # the docstring promises Y_l * (f_t' .* j_deg)', i.e.
        #   out[m, q] = Σ_i Y_l[m, i] * f_t[i, q] * j_deg[q, i]
        chunk, Q, l = 4, 3, 2
        f_t = ComplexF64[(0.7i - 0.2q) + 0.0im for i in 1:chunk, q in 1:Q]
        j_deg = Float64[0.1 * q + 0.05 * i for q in 1:Q, i in 1:chunk]
        Y_l = ComplexF64[cis(0.3m + 0.11i) * (1.0 + 0.1i) for m in 0:l, i in 1:chunk]

        out = _deg_contrib(f_t, j_deg, Y_l)
        @test size(out) == (l + 1, Q)
        for m in 1:(l + 1), q in 1:Q
            ref = sum(Y_l[m, i] * f_t[i, q] * j_deg[q, i] for i in 1:chunk)
            @test check_complex(out[m, q], ref)
        end
    end

    @testset "_deg_contrib is linear in f_t and handles complex amplitude" begin
        chunk, Q, l = 3, 2, 1
        j_deg = Float64[0.2q + 0.3i for q in 1:Q, i in 1:chunk]
        Y_l = ComplexF64[cis(0.4m - 0.2i) for m in 0:l, i in 1:chunk]
        f1 = ComplexF64[0.5i + 0.1q for i in 1:chunk, q in 1:Q]
        # a genuinely complex amplitude, to show the contraction is agnostic to it
        f2 = ComplexF64[0.3i - 0.7q + 0.9im * i for i in 1:chunk, q in 1:Q]
        a, b = -1.5, 2.25

        lhs = _deg_contrib(a .* f1 .+ b .* f2, j_deg, Y_l)
        rhs = a .* _deg_contrib(f1, j_deg, Y_l) .+ b .* _deg_contrib(f2, j_deg, Y_l)
        @test all(check_complex(lhs[k], rhs[k]) for k in eachindex(lhs))
        @test eltype(_deg_contrib(f2, j_deg, Y_l)) <: Complex
    end

    @testset "compute_B_lm matches the naive Σ_i f j_l conj(Y_lm) definition" begin
        # the single most important test in this file: an independent
        # transcription of the defining formula, atom by atom
        for lMax in 0:3
            B = compute_B_lm(pw_sph, pw_q, pw_f, lMax, UInt64(2))
            ref = pw_naive_B(pw_sph, pw_q, pw_f, lMax)
            @test size(B) == (1, pw_nrows(lMax), length(pw_q))
            for k in 1:pw_nrows(lMax), qi in eachindex(pw_q)
                @test check_complex(B[1, k, qi], ref[k, qi])
            end
        end
    end

    @testset "compute_B_lm output shape and element type" begin
        for lMax in 0:4
            B = compute_B_lm(pw_sph, pw_q, pw_f, lMax, UInt64(3))
            @test size(B) == (1, pw_nrows(lMax), length(pw_q))
            @test eltype(B) <: Complex
            @test B isa AbstractArray{<:Complex,3}
        end
    end

    @testset "compute_B_lm: B_00 = (1/sqrt(4pi)) * Σ_i f_i j_0(q r_i)" begin
        # Y_00 = 1/sqrt(4π) is real and constant, so the l = 0 multipole has a
        # closed form independent of the harmonics machinery
        B = compute_B_lm(pw_sph, pw_q, pw_f, 0, UInt64(5))
        j0(x) = x == 0.0 ? 1.0 : sin(x) / x
        for (qi, q) in enumerate(pw_q)
            ref = sum(pw_f[i, qi] * j0(q * pw_sph[1, i]) for i in 1:size(pw_sph, 2)) / sqrt(4π)
            @test check_complex(B[1, 1, qi], complex(ref))
        end
    end

    @testset "compute_B_lm channel count follows the imaginary part" begin
        f_real = real.(pw_f)
        f_cplx = pw_f .+ 0.4im .* pw_f          # genuinely complex
        f_zeroimag = ComplexF64.(f_real)        # complex type, all-zero imaginary part

        @test size(compute_B_lm(pw_sph, pw_q, f_real, 2, UInt64(2)), 1) == 1
        @test size(compute_B_lm(pw_sph, pw_q, f_cplx, 2, UInt64(2)), 1) == 2
        # the `any(x -> imag(x) != 0, ...)` branch: complex *type* is not enough,
        # a numerically zero imaginary part still collapses to one channel
        @test size(compute_B_lm(pw_sph, pw_q, f_zeroimag, 2, UInt64(2)), 1) == 1
        # a single nonzero imaginary entry anywhere flips it to two channels
        f_one = ComplexF64.(f_real); f_one[3, 2] += 1e-8im
        @test size(compute_B_lm(pw_sph, pw_q, f_one, 2, UInt64(2)), 1) == 2
    end

    @testset "compute_B_lm channels are B(Re f) and B(Im f) separately" begin
        # load-bearing claim of the whole two-channel design: the channels are
        # not a repacking of one complex B_lm, they are two independent
        # real-amplitude expansions
        f_cplx = pw_f .+ [(0.2i - 0.05q) * im for i in 1:5, q in eachindex(pw_q)]
        lMax = 3
        B = compute_B_lm(pw_sph, pw_q, f_cplx, lMax, UInt64(2))
        Bre = compute_B_lm(pw_sph, pw_q, real.(f_cplx), lMax, UInt64(2))
        Bim = compute_B_lm(pw_sph, pw_q, imag.(f_cplx), lMax, UInt64(2))
        @test size(B, 1) == 2
        @test all(check_complex(B[1, k, q], Bre[1, k, q]) for k in 1:pw_nrows(lMax), q in eachindex(pw_q))
        @test all(check_complex(B[2, k, q], Bim[1, k, q]) for k in 1:pw_nrows(lMax), q in eachindex(pw_q))
    end

    @testset "compute_B_lm is invariant to _CHUNK" begin
        # chunking changes only the order the per-atom terms are accumulated in,
        # so results are equal up to floating-point summation order
        N = size(pw_sph, 2)
        lMax = 3
        ref = compute_B_lm(pw_sph, pw_q, pw_f, lMax, UInt64(N))
        for c in (1, 2, 3, N - 1, N, N + 10, 1 << 40)
            B = compute_B_lm(pw_sph, pw_q, pw_f, lMax, UInt64(c))
            @test size(B) == size(ref)
            @test all(check_complex(B[i], ref[i]) for i in eachindex(B))
        end
    end

    @testset "compute_B_lm: _CHUNK = 1 does not underflow the UInt64 range" begin
        # `stop = start + _CHUNK - 1` with an unsigned _CHUNK, and the nested
        # `view(coords_sph, 2:3, idx)` reindexing, both used to wrap around
        # typemax(UInt64); the code converts the index to Int to avoid it
        N = size(pw_sph, 2)
        @test N > 1
        B1 = compute_B_lm(pw_sph, pw_q, pw_f, 2, UInt64(1))
        Ball = compute_B_lm(pw_sph, pw_q, pw_f, 2, UInt64(N))
        @test all(check_complex(B1[i], Ball[i]) for i in eachindex(B1))
    end

    @testset "compute_B_lm is linear in f_atoms" begin
        # both operands are real so both sides have one channel; mixing a real
        # and a complex f would compare a (1, K, Q) against a (2, K, Q)
        f1 = real.(pw_f)
        f2 = [cos(0.7i + 1.3q) for i in 1:5, q in eachindex(pw_q)]
        a, b = 2.5, -0.75
        lMax = 2
        lhs = compute_B_lm(pw_sph, pw_q, a .* f1 .+ b .* f2, lMax, UInt64(2))
        rhs = a .* compute_B_lm(pw_sph, pw_q, f1, lMax, UInt64(2)) .+ b .* compute_B_lm(pw_sph, pw_q, f2, lMax, UInt64(2))
        @test size(lhs) == size(rhs)
        @test all(check_complex(lhs[i], rhs[i]) for i in eachindex(lhs))
    end

    @testset "compute_B_lm is additive over disjoint atom sets" begin
        # B_lm is a plain sum over atoms, so splitting the molecule in two and
        # adding the parts must reproduce the whole
        lMax = 3
        left = 1:2
        right = 3:5
        Bwhole = compute_B_lm(pw_sph, pw_q, pw_f, lMax, UInt64(2))
        Bleft = compute_B_lm(pw_sph[:, left], pw_q, pw_f[left, :], lMax, UInt64(2))
        Bright = compute_B_lm(pw_sph[:, right], pw_q, pw_f[right, :], lMax, UInt64(2))
        @test all(check_complex(Bwhole[i], Bleft[i] + Bright[i]) for i in eachindex(Bwhole))
    end

    @testset "compute_B_lm: a larger lMax only appends degree blocks" begin
        # mirrors the sphHarm packing test: degree l's rows depend on l alone,
        # so truncating at a lower lMax must give the shared prefix back
        Bbig = compute_B_lm(pw_sph, pw_q, pw_f, 5, UInt64(2))
        Bsmall = compute_B_lm(pw_sph, pw_q, pw_f, 2, UInt64(2))
        @test all(check_complex(Bbig[1, k, q], Bsmall[1, k, q]) for k in 1:pw_nrows(2), q in eachindex(pw_q))
    end

    @testset "compute_B_lm: atoms on the z axis populate only the m = 0 rows" begin
        # θ = 0 kills every m > 0 harmonic, so the nonzero rows of B are exactly
        # the packed offsets partial_wave_weights assigns the weight 1.0 to --
        # an independent, behavioural check that the two agree on the layout
        lMax = 4
        zsph = to_spherical([0.0 0.0 0.0; 0.0 0.0 0.0; 1.0 2.0 3.5])
        f = ones(3, length(pw_q))
        B = compute_B_lm(zsph, pw_q, f, lMax, UInt64(2))
        w = partial_wave_weights(lMax)
        for l in 0:lMax, m in 0:l
            k = pw_idx(l, m)
            # q = 0 kills every l > 0 too, so use the largest q for the
            # "is nonzero" half of the claim
            v = B[1, k, length(pw_q)]
            if m == 0
                @test w[k] == 1.0
                @test abs(v) > 1e-6
            else
                @test w[k] == 2.0
                @test check_float(abs(v), 0.0)
            end
        end
    end

    @testset "compute_B_lm zero cases" begin
        # a zero amplitude contributes nothing at any (l, m, q)
        Bz = compute_B_lm(pw_sph, pw_q, zeros(5, length(pw_q)), 3, UInt64(2))
        @test size(Bz) == (1, pw_nrows(3), length(pw_q))
        @test all(iszero, Bz)
        # no atoms at all: the chunk loop never runs, so the (still correctly
        # shaped) accumulator comes back untouched rather than erroring
        B0 = compute_B_lm(zeros(3, 0), pw_q, zeros(0, length(pw_q)), 2, UInt64(4))
        @test size(B0) == (1, pw_nrows(2), length(pw_q))
        @test all(iszero, B0)
    end

    @testset "compute_B_lm accepts a Molecule's own spherical coordinates" begin
        # the documented input layout is exactly what Molecules.coords_spherical
        # returns so it must go straight through
        mol = create("tri", ["C", "O", "N"], [(0.0, 0.0, 1.0), (1.2, -0.3, 0.5), (-0.8, 0.9, -1.1)])
        sph = coords_spherical(mol)
        @test size(sph) == (3, 3)
        f = [1.0 + 0.1i for i in 1:3, _ in eachindex(pw_q)]
        B = compute_B_lm(sph, pw_q, f, 2, UInt64(2))
        ref = pw_naive_B(sph, pw_q, f, 2)
        @test all(check_complex(B[1, k, q], ref[k, q]) for k in 1:pw_nrows(2), q in eachindex(pw_q))
    end

    @testset "compute_B_lm exception contract" begin
        # _CHUNK = 0 would make `1:_CHUNK:N` step by zero and loop forever
        @test_throws DomainError compute_B_lm(pw_sph, pw_q, pw_f, 2, UInt64(0))
        @test_throws ArgumentError compute_B_lm(pw_sph, pw_q, pw_f, -1, UInt64(2))
        @test_throws ArgumentError compute_B_lm(pw_sph, pw_q, pw_f, -3, UInt64(2))
        # coords must be the (3, N) spherical layout, not (N, 3) or cartesian-ish
        @test_throws ArgumentError compute_B_lm(zeros(2, 5), pw_q, pw_f, 2, UInt64(2))
        @test_throws ArgumentError compute_B_lm(zeros(4, 5), pw_q, pw_f, 2, UInt64(2))
        @test_throws ArgumentError compute_B_lm(permutedims(pw_sph), pw_q, pw_f, 2, UInt64(2))
        # f_atoms is (N, Q): one row per atom, one column per q
        @test_throws ArgumentError compute_B_lm(pw_sph, pw_q, pw_f[1:4, :], 2, UInt64(2))
        @test_throws ArgumentError compute_B_lm(pw_sph, pw_q, vcat(pw_f, pw_f), 2, UInt64(2))
        @test_throws ArgumentError compute_B_lm(pw_sph, pw_q, pw_f[:, 1:2], 2, UInt64(2))
        @test_throws ArgumentError compute_B_lm(pw_sph, pw_q, hcat(pw_f, pw_f), 2, UInt64(2))
        # lMax = 0 and a single q are on the allowed side of every boundary
        @test size(compute_B_lm(pw_sph, [0.3], pw_f[:, 1:1], 0, UInt64(1))) == (1, 1, 1)
    end

    @testset "self_scatter equals its defining sum by hand" begin
        lMax = 3
        w = partial_wave_weights(lMax)
        f_cplx = pw_f .+ [(0.15i + 0.02q) * im for i in 1:5, q in eachindex(pw_q)]
        for B in (
            compute_B_lm(pw_sph, pw_q, real.(pw_f), lMax, UInt64(2)),
            compute_B_lm(pw_sph, pw_q, f_cplx, lMax, UInt64(2))
        )
            S = self_scatter(B, w)
            @test length(S) == length(pw_q)
            @test eltype(S) <: Real
            for qi in eachindex(pw_q)
                ref = 4π * sum(
                    w[k] * abs2(B[c, k, qi])
                    for c in 1:size(B, 1), k in 1:length(w)
                    )
                @test check_float(S[qi], ref)
            end
            # a weighted sum of squared magnitudes is non-negative by construction
            @test all(>=(0.0), S)
        end
    end

    @testset "self_scatter scales quadratically in B" begin
        lMax = 2
        w = partial_wave_weights(lMax)
        B = compute_B_lm(pw_sph, pw_q, real.(pw_f), lMax, UInt64(2))
        S = self_scatter(B, w)
        for α in (0.0, 0.5, -2.0, 3.0)
            Sα = self_scatter(α .* B, w)
            @test all(check_float(Sα[q], α^2 * S[q]) for q in eachindex(S))
        end
    end

    @testset "self_scatter adds the two channels incoherently" begin
        # the channels are summed, never mixed: a two-channel B's self_scatter
        # must equal the sum of the two single-channel ones
        lMax = 3
        w = partial_wave_weights(lMax)
        f_cplx = pw_f .+ [(0.3i - 0.06q) * im for i in 1:5, q in eachindex(pw_q)]
        B = compute_B_lm(pw_sph, pw_q, f_cplx, lMax, UInt64(2))
        @test size(B, 1) == 2
        S1 = self_scatter(B[1:1, :, :], w)
        S2 = self_scatter(B[2:2, :, :], w)
        S = self_scatter(B, w)
        @test all(check_float(S[q], S1[q] + S2[q]) for q in eachindex(S))
    end

    @testset "cross_scatter of a B with itself is self_scatter" begin
        lMax = 3
        w = partial_wave_weights(lMax)
        f_cplx = pw_f .+ [(0.2i + 0.01q) * im for i in 1:5, q in eachindex(pw_q)]
        for B in (
            compute_B_lm(pw_sph, pw_q, real.(pw_f), lMax, UInt64(2)),
            compute_B_lm(pw_sph, pw_q, f_cplx, lMax, UInt64(2))
        )
            X = cross_scatter(B, B, w)
            S = self_scatter(B, w)
            @test length(X) == length(pw_q)
            @test eltype(X) <: Real   # Re(B conj(B)) is real even for complex B
            @test all(check_float(X[q], S[q]) for q in eachindex(S))
        end
    end

    @testset "cross_scatter is symmetric and matches its defining sum" begin
        lMax = 2
        w = partial_wave_weights(lMax)
        fa = real.(pw_f)
        fb = [sin(0.9i - 0.4q) for i in 1:5, q in eachindex(pw_q)]
        A = compute_B_lm(pw_sph, pw_q, fa, lMax, UInt64(2))
        Bb = compute_B_lm(pw_sph, pw_q, fb, lMax, UInt64(2))
        X = cross_scatter(A, Bb, w)
        # Re(a conj(b)) = Re(b conj(a)), so the order of the operands is free
        @test all(check_float(X[q], cross_scatter(Bb, A, w)[q]) for q in eachindex(X))
        for qi in eachindex(pw_q)
            ref = 4π * sum(w[k] * real(A[1, k, qi] * conj(Bb[1, k, qi])) for k in 1:length(w))
            @test check_float(X[qi], ref)
        end
        # unlike self_scatter it is signed: negating one operand negates it
        @test all(check_float(cross_scatter(A, -1.0 .* Bb, w)[q], -X[q]) for q in eachindex(X))
    end

    @testset "cross_scatter is bilinear" begin
        lMax = 2
        w = partial_wave_weights(lMax)
        A = compute_B_lm(pw_sph, pw_q, real.(pw_f), lMax, UInt64(2))
        B1 = compute_B_lm(
            pw_sph, pw_q, [cos(0.5i + q) for i in 1:5, q in eachindex(pw_q)],
            lMax, UInt64(2)
        )
        B2 = compute_B_lm(
            pw_sph, pw_q, [0.3i - 0.2q for i in 1:5, q in eachindex(pw_q)],
            lMax, UInt64(2)
        )
        a, b = 1.75, -0.5
        lhs = cross_scatter(A, a .* B1 .+ b .* B2, w)
        rhs = a .* cross_scatter(A, B1, w) .+ b .* cross_scatter(A, B2, w)
        @test all(check_float(lhs[q], rhs[q]) for q in eachindex(lhs))
        # linearity in the first slot follows from the same identity
        lhs2 = cross_scatter(a .* B1 .+ b .* B2, A, w)
        @test all(check_float(lhs2[q], rhs[q]) for q in eachindex(lhs2))
    end

    @testset "cross_scatter with mismatched channel counts uses min(C_a, C_b)" begin
        # documented behaviour: a real amplitude has no imaginary channel, and
        # its absent channel must contribute zero rather than raise
        lMax = 2
        w = partial_wave_weights(lMax)
        f_cplx = pw_f .+ [(0.25i - 0.03q) * im for i in 1:5, q in eachindex(pw_q)]
        A2 = compute_B_lm(pw_sph, pw_q, f_cplx, lMax, UInt64(2))       # 2 channels
        B1 = compute_B_lm(pw_sph, pw_q, real.(pw_f), lMax, UInt64(2))  # 1 channel
        @test size(A2, 1) == 2 && size(B1, 1) == 1
        X = cross_scatter(A2, B1, w)
        @test length(X) == length(pw_q)
        # equals the cross of the first channels alone
        Xfirst = cross_scatter(A2[1:1, :, :], B1, w)
        @test all(check_float(X[q], Xfirst[q]) for q in eachindex(X))
        # and the operand order does not change which channels survive
        @test all(check_float(cross_scatter(B1, A2, w)[q], X[q]) for q in eachindex(X))
    end
end
