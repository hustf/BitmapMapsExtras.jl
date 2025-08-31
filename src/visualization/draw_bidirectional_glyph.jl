# Drawing functionality,
# builds on `draw_direct.jl` and `BitmapMaps.jl/mark_utils.jl`.
# Tied in with how we define a matrix format for bidirectional quantity 
# pairs in `differential_geom/curvature.jl`
# Contains `plot_principal_directions_glyph!` and `draw_bidirectional_quantity_glyph!` 
# and callee `draw_two_arrow_glyph!`

"""
    plot_principal_directions_glyph!(cov, pt, K, gs::GSTensor)

Also see `draw_bidirectional_quantity_glyph!`
"""
function plot_principal_directions_glyph!(cov, pt, K, gs::GSTensor)
    # Extract primary and / or secondary principal components
    #dK =  K[:, gs.directions]
    if is_in_limits(gs, K)
        draw_bidirectional_quantity_glyph!(cov, pt, K, gs)
    else
        # A dash instead of the out-of-limits glyph
        if gs.dashsize > 0
            spray!(cov, pt, gs.dashsize, gs.strength)
        end
    end
    cov
end






"""
    draw_bidirectional_quantity_glyph!(cov, p, K, strength)
    draw_bidirectional_quantity_glyph!(cov, p, v::AbstractVector, strength)

Draws two-arrow glyphs on an image at position `p` for principal bi-directional quantities
specified by the matrix `K`.

# Arguments
- `cov`: The coverage matrix to render the glyphs onto.
- `p`: The location of the glyphs in 2D space.
- `K`: A principal bi-directional 4x4 matrix where `K[:,1]` specifies the first principal direction,
  magnitude, and sign, and `K[:,2]` specifies the second. See [`draw_two_arrow_glyph!`](@ref).
- `strength`: the maximum coverage applied to pixels

# Details
Due to symmetry, only 180° (π) is needed for the direction of each quantity. The additional
information is used to determine if the quantity is negative (arrows toward `p`) or positive
(arrows from `p`).

# Notes
- As an example, 2D principal stress directions are always orthonormal, whereas in `K`, the
  directions need not be orthonormal. If an orthonormal glyph is projected onto a plane that is
  rotated with respect to the stress plane, the projected axes would not be orthonormal. This
  function could be used to draw such a projection.
- This function modifies `img` in-place.

# Examples
```julia
img = zeros(100, 100)  # Example image matrix
p = [50.0, 50.0]       # Center position
K = rand(4, 4)         # Example principal bi-directional matrix
draw_bidirectional_quantity_glyph!(img, p, K)  # Draws two glyphs
```
"""
function draw_bidirectional_quantity_glyph!(cov, p, K, gs)
    @assert size(K) == (2, 2)
    if 1 ∈ gs.directions
        draw_two_arrow_glyph!(cov, p, gs.multip .* view(K, :, 1), gs.strength)
    end
    if 2 ∈ gs.directions
        draw_two_arrow_glyph!(cov, p, gs.multip .* view(K, :, 2), gs.strength)
    end
    cov
end
function draw_bidirectional_quantity_glyph!(cov, p, v::AbstractVector, gs)
    throw("Dead")
    @assert size(v, 1) == 2
    draw_two_arrow_glyph!(cov, p, gs.multip .* v, strength)
end

"""
    draw_two_arrow_glyph!(cov, p, v::AbstractVector, strength::Float32)

img is the image matrix
p is the location
F is the length of both arrows. If negative, arrows end at p.
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