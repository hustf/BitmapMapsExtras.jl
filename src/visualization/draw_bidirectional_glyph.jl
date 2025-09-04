# Drawing functionality,
# builds on `draw_direct.jl` and `BitmapMaps.jl/mark_utils.jl`.
# Tied in with how we define a matrix format for bidirectional quantity 
# pairs in `differential_geom/curvature.jl`
# Contains `plot_principal_directions_glyph!` and `draw_bidirectional_quantity_glyph!` 
# and callee `draw_two_arrow_glyph!`

"""
    plot_principal_directions_glyph!(cov, pt, K, gs::GSTensor)
    plot_principal_directions_glyph!(vcov::Vector{Matrix{Float32}}, pt, K, gs::GSTensor)

Also see `draw_bidirectional_quantity_glyph!`
"""
function plot_principal_directions_glyph!(cov, pt, K, gs::GSTensor)
     @assert length(gs.directions) == 1
    # Extract primary and / or secondary principal components
    if is_in_limits(gs, K)
        v = gs.multip .* view(K, :, first(gs.directions))
        draw_bidirectional_quantity_glyph!(cov, pt, v, gs.strength)
    else
        # A dash instead of the out-of-limits glyph
        if gs.dashsize > 0
            spray!(cov, pt, gs.dashsize, gs.strength)
        end
    end
    cov
end

function plot_principal_directions_glyph!(vcov::Vector{Matrix{Float32}}, pt, K, gs::GSTensor)
     @assert gs.directions == 1:2
     @assert length(vcov) == 2
    # Extract primary and / or secondary principal components
    if is_in_limits(gs, K)
        v1 = gs.multip .* view(K, :, 1)
        draw_bidirectional_quantity_glyph!(vcov[1], pt, v1, gs.strength)
        v2 = gs.multip .* view(K, :, 2)
        draw_bidirectional_quantity_glyph!(vcov[2], pt, v2, gs.strength)
    else
        # A dash instead of the out-of-limits glyph
        if gs.dashsize > 0
            spray!(vcov[1], pt, gs.dashsize, gs.strength)
            spray!(vcov[2], pt, gs.dashsize, gs.strength)
        end
    end
    vcov
end


"""
    draw_bidirectional_quantity_glyph!(cov, p, v, strength)

This function modifies `cov` in-place. Called by `plot_principal_directions_glyph!`

Spray a two-arrow glyph 'v' on an image at position `p`.

# Arguments
- `cov`: The coverage matrix. Each element (pixel) accumulate coverage, how fast depends on strength.
- `p`: The location of the glyph in 2D space.
- `v`: A 2d vector. First and second quadrant: Positive sign. The opposite direction is always drawn, too. 
    See [`draw_two_arrow_glyph!`](@ref).
- `strength`: the maximum coverage applied to pixels per 'spray dash'. The centre is often hit with more such dashes.

# Notes
- Due to symmetry, only 180° (π) is needed for the direction of each quantity. The additional
information is used to determine if the quantity is negative (arrows toward `p`) or positive
(arrows from `p`).
-If an orthonormal glyph (representing e.g. planar stress or surface curvature) is projected onto a plane that is
  rotated with respect to the (stress, curvature) plane, the projected axes will not be orthonormal. Each principal 
  component is still 180° symmetric.
- The sign convention for the y-axis is what it is... See the example.

# Examples
```julia
julia> begin # The shape is hardly visble on this small scale.
    cov = zeros(Float32, 11, 11); # Example coverage matrix
    p = CartesianIndex((6, 6))    # Center position
    v = [4.75, 3.75]              # Example bi-directional vector
    strength = 0.5f0
    draw_bidirectional_quantity_glyph!(cov, p, v, strength)
end
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.5
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.5  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.5  0.5  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.5  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.5  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.5  0.5  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.5  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.5  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
draw_bidirectional_quantity_glyph!(cov, p, v, strength) =  draw_two_arrow_glyph!(cov, p, v, strength)

"""
    draw_two_arrow_glyph!(cov, p, v::AbstractVector, strength::Float32)

See `draw_bidirectional_quantity_glyph!`.

v[1:2] specifies both direction and sign. 
    0-π:  Positive sign, arrows outward.   
    π-2π: Negative sign: arrows inward.   
"""
function draw_two_arrow_glyph!(cov, p, v::AbstractVector, strength::Float32)
    Δj = Int(round(v[1]))
    Δi = Int(round(v[2]))
    if is_bidirec_vect_positive(v)
        # First or second quadrant. Positive.
        draw_vector!(cov, p, -Δi, Δj, strength)
        draw_vector!(cov, p, Δi, -Δj, strength)
    else
        # Third or fourth quadrant. Negative.
        draw_vector!(cov, p + CartesianIndex(-Δi, Δj), Δi, -Δj, strength)
        draw_vector!(cov, p + CartesianIndex(Δi, -Δj), -Δi, Δj, strength)
    end
end