# Functions that both calculate and draw glyphs,
# one at a time instead of allocating result matrices first.
# Contains 
# - plot_glyphs 
# - plot_glyphs!
# - plot_glyphs_given_values!
# - apply_color_by_coverage!
# - coverage_fitting_image
# Also, extensions for 
# - coverage_fitting_image
# - apply_color_by_coverage!

#######################
# Plot glyphs
#######################

"""
    plot_glyphs(z::Matrix{<:AbstractFloat}, pts, gs::AbstractGlyphSpec)
    plot_glyphs(b::AbstractIJFunctor, pts, gs::AbstractGlyphSpec)

Plot the glyph specificed by `gs` at every `pts`. 

Values and image size are based on `z` or `b` (namely, `b.z`).
"""
function plot_glyphs(z::Matrix{<:AbstractFloat}, pts, gs::AbstractGlyphSpec)
    # Allocate an empty color image (since user didn't supply one)
    img = zeros(RGBA{N0f8}, size(z)...)
    # Modify the image
    plot_glyphs!(img, z, pts, gs)
end
function plot_glyphs(b::AbstractIJFunctor, pts, gs::AbstractGlyphSpec)
    # Allocate an empty color image (since user didn't supply one)
    img = zeros(RGBA{N0f8}, size(b.z)...)
    # Modify the image
    plot_glyphs!(img, b, pts, gs)
end
"""
    plot_glyphs!(img::Matrix{<:Colorant}, z::Matrix{<:AbstractFloat}, pts, gs::AbstractGlyphSpec)
    plot_glyphs!(img::Matrix{<:Colorant}, b::AbstractIJFunctor, pts, gs::AbstractGlyphSpec)
    plot_glyphs!(cov::T, z::Matrix{<:AbstractFloat}, pts, gs::GSTangentBasis) where 
        {T<:Vector{Matrix{Float32}}}
    plot_glyphs!(cov::T, z::Matrix{<:AbstractFloat}, pts, gs) where 
        T<:Union{Matrix{Float32}, Vector{Matrix{Float32}}}
    plot_glyphs!(cov::T, b::AbstractIJFunctor, pts, gs::AbstractGlyphSpec) where 
        {T<:Union{Matrix{Float32}, Vector{Matrix{Float32}}}}

Modify the image `img` by plotting the glyph specificed at every `pts`. 

Values that the glyps represent are based on `z`, or by `b`.

Methods with the `cov` argument can be considered internal. `cov` stands for coverage in a specific color.
"""
function plot_glyphs!(img::Matrix{<:Colorant}, z::Matrix{<:AbstractFloat}, pts, gs::AbstractGlyphSpec)
    # This layer in the call hierarchy: img-> cov -> img
    cov = coverage_fitting_image(img, gs)
    # Modify cover buffer(s)
    plot_glyphs!(cov, z, pts, gs)
    # Cover ('unlimited coverage') maped to colors
    apply_color_by_coverage!(img, cov, gs)
    img
end
function plot_glyphs!(img::Matrix{<:Colorant}, b::AbstractIJFunctor, pts, gs::AbstractGlyphSpec)
    # This layer in the call hierarchy: img-> cov -> img
    cov = coverage_fitting_image(img, gs)
    # Modify cover buffer(s)
    plot_glyphs!(cov, b, pts, gs)
    # Cover ('unlimited coverage') maped to colors
    apply_color_by_coverage!(img, cov, gs)
    img
end
function plot_glyphs!(cov::T, z::Matrix{<:AbstractFloat}, pts, gs::GSTangentBasis) where 
    {T<:Vector{Matrix{Float32}}}
    # We don't define a functor for tangent basis, unlike what we do for other glyphs.
    Ri, Ω, v, P, _, _, _, _, _ = allocations_curvature(CartesianIndices(z))
    # Plot tangent basis for internal points one at a time
    for pt in pts
        # Find P in-place
        tangent_basis!(P, v, view(z, Ω .+ pt))
        # Plot the single glyph
        plot_glyph_given_value!(cov, pt, P, gs)
    end
    cov
end
function plot_glyphs!(cov::T, z::Matrix{<:AbstractFloat}, pts, gs) where 
    T<:Union{Matrix{Float32}, Vector{Matrix{Float32}}}
    # This makes a default functor for z (containing the z values) and passes it on.
    plot_glyphs!(cov, default_ij_functor(z, gs), pts, gs)
end
function plot_glyphs!(cov::T, b::AbstractIJFunctor, pts, gs::AbstractGlyphSpec) where 
    {T<:Union{Matrix{Float32}, Vector{Matrix{Float32}}}}
    # This plots a glyph at each pt without storing the values.
    # No domain checking.
    for pt in pts
        value = b(pt.I...)
        plot_glyph_given_value!(cov, pt, value, gs)
    end
    cov
end



##########################
# Plot glyphs given values
##########################

"""
    plot_glyphs_given_values!(img, pts, values, gs::AbstractGlyphSpec)
    plot_glyphs_given_values!(cov::T, pts, values, gs::AbstractGlyphSpec) where 
        {T<:Union{Matrix{Float32}, Vector{Matrix{Float32}}}}
"""
function plot_glyphs_given_values!(img, pts, values, gs::AbstractGlyphSpec)
    # This layer in the call hierarchy: img-> cov -> img
    cov = coverage_fitting_image(img, gs)
    # Modify cover buffer(s)
    plot_glyphs_given_values!(cov, pts, values, gs)
    # Cover ('unlimited coverage') maped to colors
    apply_color_by_coverage!(img, cov, gs)
    img
end
function plot_glyphs_given_values!(cov::T, pts, values, gs::AbstractGlyphSpec) where 
    {T<:Union{Matrix{Float32}, Vector{Matrix{Float32}}}}
    @assert length(values) == length(pts)
    # This plots a glyph at each pt without storing the value
    for (pt, value) in zip(pts, values)
        plot_glyph_given_value!(cov, pt, value, gs)
    end
    cov
end


# These methods are placed here because they use an 'internal' type,
# whereas `draw_direct.jl` is more general and might be moved to 
# a separate package.
"""
    apply_color_by_coverage!(img, cov, gs::AbstractGlyphSpec)
    apply_color_by_coverage!(img, cov, gs::GSTensor{<:Any, 1})
    apply_color_by_coverage!(img, cov, gs::GSTensor{<:Any, 2})
    apply_color_by_coverage!(img, cov, gs::GSTangentBasis)
"""
function apply_color_by_coverage!(img, cov, gs::AbstractGlyphSpec)
    apply_color_by_coverage!(img, cov, gs.color)
end
function apply_color_by_coverage!(img, cov, gs::GSTensor{<:Any, 1})
    apply_color_by_coverage!(img, cov, gs.colors[1])
end
function apply_color_by_coverage!(img, cov, gs::GSTensor{<:Any, 2})
    apply_color_by_coverage!(img, cov[1], gs.colors[1])
    apply_color_by_coverage!(img, cov[2], gs.colors[2])
end
function apply_color_by_coverage!(img, cov, gs::GSTangentBasis)
    # Note that we don't 'spray' this glyph type. It's full coverage
    # or none.
    apply_color_by_any_coverage!(img, cov[1], gs.colors[1])
    apply_color_by_any_coverage!(img, cov[2], gs.colors[2])
    apply_color_by_any_coverage!(img, cov[3], gs.colors[3])
    img
end



"""
    coverage_fitting_image(img, gs::AbstractGlyphSpec)
    coverage_fitting_image(img, gs::GSTensor{<:Any, 2})
    coverage_fitting_image(img, gs::GSTangentBasis)

Output one, two or three color-separate matrices, depending on the
glyph specification. 
"""
coverage_fitting_image(img, gs::AbstractGlyphSpec) = zeros(Float32, size(img)...)
function coverage_fitting_image(img, gs::GSTensor{<:Any, 2})
    [zeros(Float32, size(img)...), zeros(Float32, size(img)...)]
end
function coverage_fitting_image(img, gs::GSTangentBasis)
    [zeros(Float32, size(img)...), zeros(Float32, size(img)...), zeros(Float32, size(img)...)]
end
