using .SphFuncs: sphHarm, sphBess

"""
    _deg_contrib(f_t, j_deg, Y_l) -> AbstractMatrix{<:Number}

Per-degree `B_lm` contribution: `Y_l * (f_t' .* j_deg)'`.

Computes the contraction of the spherical-harmonic matrix for degree `l`
against the elementwise product of the (transposed) per-atom amplitude and
the degree-`l` spherical Bessel factors.

# Arguments
    -   `f_t::AbstractMatrix`: `(chunk, Q)` — per-atom amplitude for this chunk
        (real or imaginary channel), evaluated on the `q`-grid. May be complex
        (both channels are cast to a common complex element type upstream so
        the contraction is uniform regardless of channel).
    -   `j_deg::AbstractMatrix`: `(Q, chunk)` — `j_l(q·r)` for this degree, laid
        out the same way as `f_t'` (i.e. `(Q, chunk)`). Always real-valued.
    -   `Y_l::AbstractMatrix`: `(l+1, chunk)` — the `l+1` orders `m = 0, …, l` of
        `Y_l^m` for every atom in the chunk. Complex-valued.

# Returns
    -   `AbstractMatrix{<:Number}` of size `(l+1, Q)`: `Y_l * (f_t' .* j_deg)'`.
        Complex-valued whenever `f_t` or `Y_l` is complex, matching their
        promoted element type.

# Details
    - `W = f_t' .* j_deg` — elementwise product, `(Q, chunk)`.
    - `W'` — transpose back to `(chunk, Q)`.
    - `Y_l * W'` — matrix product, `(l+1, chunk) * (chunk, Q) → (l+1, Q)`.

Kept as a free function (not a closure/method on a struct) purely for
symmetry with the code this was ported from. Unlike that version, there is
no separate "compiled" variant here: Julia specializes each method on
argument *types*, not array *shapes*, so this already compiles once per
element-type combination and is reused across every chunk size and every
degree without a "one graph per degree" blowup.
"""
function _deg_contrib(
    f_t::AbstractMatrix, j_deg::AbstractMatrix, Y_l::AbstractMatrix
)::AbstractMatrix{<:Number}
    W = transpose(f_t) .* j_deg   # (Q, chunk) .* (Q, chunk) -> (Q, chunk)
    return Y_l * transpose(W)     # (l+1, chunk) * (chunk, Q) -> (l+1, Q)
end

"""
    partial_wave_weights(lMax) -> Vector{Float64}

±m symmetry weights: `m = 0 -> 1`, `m > 0 -> 2`.

Valid only against a `B_lm` built from a REAL per-atom amplitude, where
`B_{l,-m} = (-1)^m * conj(B_lm)`. `compute_B_lm` guarantees that by
splitting a complex `f` into real/imaginary channels, each of which is
separately real-amplitude; `self_scatter`/`cross_scatter` then sum the
channels. Folding a complex `f` directly is NOT rotationally invariant
(measured 3% error with Fe at 8 keV).

# Arguments
- `lMax::Integer`: maximum spherical harmonic degree. Must be non-negative.

# Returns
    -   `Vector{Float64}` of length `(lMax+1)(lMax+2)/2`, one weight per stored
        `(l, m)` pair with `m = 0, …, l`, flattened in the same
        `(0,0), (1,0), (1,1), (2,0), (2,1), (2,2), …` degree-major order used to
        index `B_lm` elsewhere (i.e. `k0 = l*(l+1)/2` gives the offset of degree `l`'s block).

# Mathematical derivation

We want:
```
S(q) = 4π * Σ_l Σ_{m=-l}^{l} |B_lm(q)|²
```

with
```
B_lm(q) = Σ_i f_i(q) * j_l(q r_i) * conj(Y_lm(θ_i, φ_i))
```

Spherical harmonics satisfy `Y_{l,-m} = (-1)^m * conj(Y_lm)`. Plugging into
`B_{l,-m}`:
```
B_{l,-m} = Σ_i f_i * j_l * conj(Y_{l,-m})
         = Σ_i f_i * j_l * conj[(-1)^m * conj(Y_lm)]
         = (-1)^m * Σ_i f_i * j_l * Y_lm
```

If `f_i` is real, that last sum is just `conj(B_lm)`, since conjugating a
real-coefficient sum of `Y_lm`'s conjugates flips it back to `Y_lm`:
```
B_{l,-m} = (-1)^m * conj(B_lm)   =>   |B_{l,-m}|² = |B_lm|²
```

Every negative-`m` term is a free duplicate of its positive-`m` partner, so:
```
Σ_{m=-l}^{l} |B_lm|²    = |B_l0|² + Σ_{m=1}^{l} (|B_lm|² + |B_{l,-m}|²)
                        = 1·|B_l0|² + Σ_{m=1}^{l} 2·|B_lm|²
```
which is exactly the `1`-for-`m=0`, `2`-for-`m>0` weighting returned here.
"""
function partial_wave_weights(lMax::Integer)::Vector{Float64}
    lMax < 0 && throw(ArgumentError("partial_wave_weights: lMax must be non-negative"))
    return [2.0 - (m == 0) for deg in 0:lMax for m in 0:deg]
end

"""
    compute_B_lm(coords_sph, qvals, f_atoms, lMax; _CHUNK, backend=Array)
    -> AbstractArray{<:Complex,3}

Compute `B_lm(q) = Σ_i f_atoms[i](q) * j_l(q*r_i) * conj(Y_lm(θ_i, φ_i))`.

`f_atoms` carries whatever per-atom complex scattering amplitude the
caller wants (element form factors for the vacuum term, dummy-atom
excluded-volume/shell amplitudes for the other terms) — this is the
shared low-level primitive every `S_(a,b)` term downstream is built from.

# Arguments
    -   `coords_sph::AbstractMatrix{<:Real}`, size `(3, N)`: per-atom spherical
        coordinates in the column-per-atom layout `Molecules.coords_spherical`
        returns — row 1 is `r`, row 2 is `θ`, row 3 is `φ`. Passed straight
        through; the `(θ, φ)` rows go to `sphHarm` as a `(2, chunk)` block and
        the `r` row to `sphBess`, with no unpacking into loose vectors.
    -   `qvals::AbstractVector{<:Real}`, length `Q`: momentum transfer grid.
    -   `f_atoms::AbstractMatrix{<:Number}`, size `(N, Q)`: per-atom scattering
        amplitude, already evaluated on `qvals`. Real or complex — a purely
        real matrix is handled directly, with no need to pre-cast it to complex.
    -   `lMax::Integer`: maximum spherical harmonic degree. Must be non-negative.
    -   `_CHUNK::UInt64`: number of atoms to process in one pass.
    -   `backend::Type{<:AbstractArray} = Array`: array type used for the
        per-chunk compute buffers. Pass e.g. `CUDA.CuArray` to run on GPU.

# Returns
    -   `AbstractArray{<:Complex,3}` of size `(C, (lMax+1)(lMax+2)÷2, Q)`, backed
        by `backend`. `C = 1` when `f_atoms` is real, `C = 2` when it has a
        nonzero imaginary part (anomalous `f''`): channel 1 from `Re(f_atoms)`,
        channel 2 from `Im(f_atoms)`. The two channels add incoherently in the
        `m`-summed invariant (see [`self_scatter`](@ref)/[`cross_scatter`](@ref))
    —   they must not be recombined into one complex `B_lm`.
"""
function compute_B_lm(
    coords_sph::AbstractMatrix{<:Real},
    qvals::AbstractVector{<:Real},
    f_atoms::AbstractMatrix{<:Number},
    lMax::Integer,
    _CHUNK::UInt64,
    backend::Type{<:AbstractArray}=Array
)::AbstractArray{<:Complex,3}

    # _CHUNK == 0 would make the `1:_CHUNK:N` range below step by zero,
    # looping forever instead of raising.
    if _CHUNK == 0
        throw(DomainError("_CHUNK must be > 0"))
    end

    lMax < 0 && throw(ArgumentError("compute_B_lm: lMax must be non-negative"))

    # `coords_sph` is the column-per-atom spherical form straight from
    # `Molecules.coords_spherical`: 3 rows, `(r, θ, φ)` in that order.
    size(coords_sph, 1) == 3 || throw(ArgumentError(
        "compute_B_lm: coords_sph must be a (3, N) matrix with rows (r, θ, φ), " *
        "as returned by Molecules.coords_spherical; got $(size(coords_sph, 1)) rows"))
    N = size(coords_sph, 2)
    r = view(coords_sph, 1, :)   # (θ, φ) are sliced per-chunk straight from `coords_sph`

    Q = length(qvals)
    size(f_atoms, 1) == N ||
        throw(ArgumentError("compute_B_lm: f_atoms must have N rows matching coords_sph's columns"))
    size(f_atoms, 2) == Q ||
        throw(ArgumentError("compute_B_lm: f_atoms must have Q columns matching qvals"))

    # Packed (l,m) row count for m = 0..l only; see `partial_wave_weights` for
    # why the m < 0 half never needs to be stored.
    N_reduced = (lMax + 1) * (lMax + 2) ÷ 2
    
    # A purely real f_atoms needs only the Re(f) channel; an imaginary part
    # (anomalous f'') needs a second channel for the ±m symmetry.
    n_chan = any(x -> imag(x) != 0, f_atoms) ? 2 : 1
    B_lm = backend(zeros(ComplexF64, n_chan, N_reduced, Q))

    # Atoms are processed `_CHUNK` at a time rather than all at once: `Y`/`j`
    # below are (N_reduced, N) / (lMax+1, Q, N) in the worst case, so doing
    # every atom in one shot would allocate an O(N*Q*lMax) buffer per array
    # even though the actual output `B_lm` is only O(lMax^2 * Q). Chunking
    # bounds that intermediate memory to O(_CHUNK*Q*lMax) regardless of N.
    for start in 1:_CHUNK:N
        stop = min(start + _CHUNK - 1, N)
        # `_CHUNK` is a `UInt64`, so `start:stop` would be a UInt64 range;
        # reindexing the nested `view(coords_sph, 2:3, idx)` against one
        # underflows to `typemax(UInt64)` and throws. Keep the index Int.
        idx = Int(start):Int(stop)

        # rows 2:3 of `coords_sph` are (θ, φ): handed to `sphHarm` as a
        # (2, chunk) angle block, not split into loose vectors.
        Y = sphHarm(Int(lMax), view(coords_sph, 2:3, idx))  # (N_reduced, chunk), complex
        j = sphBess(view(r, idx), qvals, Int(lMax))         # (lMax+1, Q, chunk), real

        # `_deg_contrib` is called lMax+1 times per channel per chunk,
        # so hoisting the conjugation out of that loop saves (lMax+1)x work.
        Y_dev = backend(ComplexF64.(conj.(Y)))
        j_dev = backend(Float64.(j))

        # Re(f) and Im(f) are each a real amplitude, so each channel's
        # B_lm obeys the ±m conjugate symmetry `partial_wave_weights` assumes.
        for chan in 1:n_chan
            part = chan == 1 ? real.(view(f_atoms, idx, :)) : imag.(view(f_atoms, idx, :))
            f_t = backend(ComplexF64.(part))  # (chunk, Q)

            for deg in 0:lMax
                
                # k0 is the packed-row offset of degree `deg`'s (m=0..deg)
                # block; matches the offset `partial_wave_weights` assumes.
                k0 = deg * (deg + 1) ÷ 2
                rows = (k0 + 1):(k0 + deg + 1)
                Y_l = view(Y_dev, rows, :)          # (l+1, chunk)
                j_deg = view(j_dev, deg + 1, :, :)  # (Q, chunk)
                B_lm[chan, rows, :] .+= _deg_contrib(f_t, j_deg, Y_l)
            end
        end
    end

    return B_lm
end

"""
    self_scatter(B_lm, weights) -> AbstractVector{<:Real}

`S(q) = 4π * Σ_c Σ_lm w_lm * |B_lm(q)|²`.

Channels (real/imaginary amplitude) add incoherently: the cross terms
between them cancel identically once summed over the full `-l..l` range of`m`.

# Arguments
- `B_lm::AbstractArray{<:Complex,3}`, size `(C, K, Q)`: as returned by [`compute_B_lm`](@ref).
- `weights::AbstractVector{<:Real}`, length `K`: as returned by [`partial_wave_weights`](@ref).

# Returns
- `AbstractVector{<:Real}` of length `Q`.
"""
function self_scatter(
    B_lm::AbstractArray{<:Complex,3}, weights::AbstractVector{<:Real}
)::AbstractVector{<:Real}
    
    # weights is (K,); reshape to (1, K, 1) so it broadcasts against B_lm's
    # (C, K, Q) over the channel and q axes without being repeated by hand.
    w = reshape(weights, 1, :, 1)
    
    # abs2.(B_lm) is |B_lm|^2 per (channel, l/m, q) entry; summing over dims
    # (1, 2) collapses channel and l/m, leaving one number per q.
    # `vec` drops the resulting (1, 1, Q) down to a plain (Q,) vector.
    return vec(sum(w .* abs2.(B_lm); dims=(1, 2))) .* (4π)
end

"""
    cross_scatter(B_lm_a, B_lm_b, weights) -> AbstractVector{<:Real}

`S(q) = 4π * Σ_c Σ_lm w_lm * Re(B_a(q) * conj(B_b(q)))`.

Only channels present in both operands contribute; a real amplitude has no imaginary 
channel (`C = 1`), so its missing channel contributes zero rather than erroring.

# Arguments
- `B_lm_a::AbstractArray{<:Complex,3}`, size `(C_a, K, Q)`.
- `B_lm_b::AbstractArray{<:Complex,3}`, size `(C_b, K, Q)`.
- `weights::AbstractVector{<:Real}`, length `K`: as returned by [`partial_wave_weights`](@ref).

# Returns
- `AbstractVector{<:Real}` of length `Q`.
"""
function cross_scatter(
    B_lm_a::AbstractArray{<:Complex,3},
    B_lm_b::AbstractArray{<:Complex,3},
    weights::AbstractVector{<:Real},
)::AbstractVector{<:Real}
    
    # Only the channels both operands actually have can be paired up; a real
    # amplitude's missing second channel would otherwise have nothing to mul against.
    n_chan = min(size(B_lm_a, 1), size(B_lm_b, 1))
    
    # Re(B_a * conj(B_b)) per (channel, l/m, q) entry, matching `self_scatter`
    # with |B_lm|^2 (= Re(B_lm * conj(B_lm))) generalised to two operands.
    cross = real.(view(B_lm_a, 1:n_chan, :, :) .* conj.(view(B_lm_b, 1:n_chan, :, :)))
    w = reshape(weights, 1, :, 1)  # see `self_scatter` for this broadcast shape
    return vec(sum(w .* cross; dims=(1, 2))) .* (4π)
end