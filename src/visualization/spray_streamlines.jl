# Functions that plot streamlines.
#
# Calls 
# - `solve_streamlines.jl`
# - `extract_streamlines.jl`
# - DrawAndSpray.jl

###################
# Spray streamlines
###################

"""
    plot_streamlines(fxy::AbstractXYFunctor, pts; 
        stroke = Stroke(), odekws...)
    --> Matrix{RGBA{N0f8}}

# Arguments

- `fxy` determines which curve to follow, and also the size of the output image. 
  See ``
- `pts` is a vector of starting points. See `indices_on_grid`.
- `stroke` specifies the appearance of each streamline.
- `odekws` is a list of keywords for the ordinary differential equation. `dtmax = 1.0` is an example: A smaller timestep can solve many problems.

# Example, 2d vector field

```
julia> using BitmapMapsExtras: TestMatrices

julia> vaxy = Vec2AtXY(𝐧ₚ!, z_cylinder(π / 6));

julia> img = background(vaxy)

julia> img = plot_streamlines!(img, vaxy, pts; stroke = Stroke(color = PALETTE_GRGB[2]), dtmax = 1.0)
```

# Example, 2d vector selected from tensor map
"""
function plot_streamlines(fxy::AbstractXYFunctor, pts; 
        stroke = Stroke(), odekws...)
    # Allocate an empty color image (since user didn't supply one)
    img = zeros(RGBA{N0f8}, size(fxy)...)
    # Modify the image
    plot_streamlines!(img, fxy, pts; stroke, odekws...)
end

"""
    plot_streamlines!(img, fxy::AbstractXYFunctor, pts; 
        stroke = Stroke(), odekws...)

See `plot_streamlines`. Use this method when you already have an image.
"""
function plot_streamlines!(img, fxy::AbstractXYFunctor, pts; 
    stroke = Stroke(), odekws...)
    # Coverage buffer
    cov = zeros(Float32, size(img)...)
    # Modify coverage
    plot_streamlines!(cov, fxy, pts; stroke, odekws...)
    # Apply color by coverage
    apply_color_by_coverage!(img, cov, stroke.color)
    img
end

"""
    plot_streamlines!(cov::Array{Float32}, fxy::AbstractXYFunctor, pts; 
        stroke = Stroke(), odekws...)
"""
function plot_streamlines!(cov::Array{Float32}, fxy::AbstractXYFunctor, pts; 
        stroke = Stroke(), odekws...)
    # Indexed points on the streamlines
    spts = get_streamlines_points(fxy, pts, stroke.tdensity; odekws...)
    # Plot points
    spray_along_nested_indices!(cov, spts, stroke.r, stroke.strength)
    cov
end

