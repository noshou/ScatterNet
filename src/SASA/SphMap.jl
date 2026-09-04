using Roots
using GLMakie

"Cubic `x^3 - x - 1`; its real root in `(1, 2)` is the plastic ratio."
f(x) = x^3 - x - 1

"Fractional part of `x`, i.e. `x - floor(x)`, in `[0, 1)`."
frac(x) = x - floor(x)

# plastic ratio ≈ 1.3247
const PLASTIC_RATIO 	= find_zero(f, (1, 2))
const PLASTIC_RATIO_SQR	= PLASTIC_RATIO ^ 2

"""
    plastic_seq_sph(start::Int64, n::Int64) -> Matrix{Float64}

Terms `start+1 .. start+n` of the 1-based 2-D plastic (R2) sequence, mapped onto
the unit sphere: the 2-D box is rolled into a cylinder, then squashed inward to
wrap the sphere. Returns an `n x 2` matrix; column 1 is azimuth `phi` in
`[0, 2pi)`, column 2 is height `z` in `[-1, 1)`.

# Arguments
- `start`: number of leading sequence terms to skip; must be `>= 0`.
- `n`: number of terms to return; must be `>= 0`.
"""
function plastic_seq_sph(start::Int64, n::Int64)::Matrix{Float64}
    
    # 0. Assertion checks 
    if start < 0
        throw(DomainError(start, "start must be >= 0"))
    elseif n < 0
        throw(DomainError(n, "n must be >= 0"))
    end
    
    # 1. Pre-allocate an n x 2 matrix
    seq = Matrix{Float64}(undef, n, 2)
    
    # 2. Populate the matrix
    for (row_idx, i) in enumerate(start+1:start+n)
        seq[row_idx, 1] = 2 * pi * (frac(i / PLASTIC_RATIO))      
        seq[row_idx, 2] = 2 * frac(i / PLASTIC_RATIO_SQR) - 1  
    end

    return seq
end

"""
    vis_plastic_seq_2D(n_points) -> Nothing

Scatter-plot `n_points` of the plastic sequence on the unit sphere with GLMakie,
blocking until the window is closed.

Run with:
    ```
        julia --project=. -e 'include("src/SASA/SphMap.jl"); vis_plastic_seq_2D(parse(Int, ARGS[1]))' <num_points>
    ```

# Arguments
- `n_points`: number of sequence points to generate and draw.
"""
function vis_plastic_seq_2D(n_points)
    
    # 1. Generate the 2D sequence data
    seq = plastic_seq_sph(0, n_points)

    # 2. Project from cylinder to unit sphere 3D coordinates
    phi = @view seq[:, 1]
    z = @view seq[:, 2]
    r = @. sqrt(max(0.0, 1.0 - z^2))

    x = @. r * cos(phi)
    y = @. r * sin(phi)

    # 3. Plot in 3D using GLMakie
    fig = Figure()
    ax = Axis3(fig[1, 1], title = "Plastic Sequence on Unit Sphere", aspect = :data)
    scatter!(ax, x, y, z, color = z, colormap = :viridis, markersize = 6)
    wait(display(fig))
end