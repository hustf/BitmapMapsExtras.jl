# Functions that both calculate and draw glyphs,
# one at a time instead of allocating result matrices 
# first.
# Contains `plot_tangent_basis_glyphs`,
# `plot_tangent_basis_glyphs!`
# `plot_curvature_glyphs`,
# `plot_curvature_glyphs!`
# Refers  RED_GREEN_BLUE and COLOR_CURVGLYPH.


####################
# Plot tangent basis
####################
#=
"""
    plot_tangent_basis_glyphs(z, pts;  gs = GSTangentBasis())
"""
function plot_tangent_basis_glyphs(z, pts; gs = GSTangentBasis())
    # Allocate an empty colorful image
    img = zeros(RGBA{N0f8}, size(z)...)
    # Modify the image
    plot_tangent_basis_glyphs!(img, z, pts; gs)
end

"""
    plot_tangent_basis_glyphs!(img, z, pts; gs = GSTangentBasis())
"""
function plot_tangent_basis_glyphs!(img, z, pts; gs = GSTangentBasis())
    throw("dead")
    Ri, Ω, v, P, _, _, _, _, _ = allocations_curvature(CartesianIndices(z))
    # Black-white buffer
    bbuf = Array{Gray{Bool}}(falses( size(img)...))
    # Plot tangent basis for internal points one at a time
    for pt in filter(pt -> pt ∈ Ri, sort(vec(pts)))
        # Find P in-place
        tangent_basis!(P, v, view(z, Ω .+ pt))
        # Plot the single glyph
        plot_orthonormal_basis_glyph!(img, bbuf, pt, P, gs)
    end
    img
end

"""
    plot_orthonormal_basis_glyph!(img, bbuf::Array{Gray{Bool}}, pt, P, gs)
"""
function plot_orthonormal_basis_glyph!(img, bbuf::Array{Gray{Bool}}, pt, P, gs)
    throw("dead")
    @assert size(img) == size(bbuf) 
    for i in 1:3
        # Blank bw-buffer
        fill!(bbuf, Gray{Bool}(false))
        perps = setdiff([1,2,3], i)
        𝐞_perp1 = P[:, perps[1]]
        𝐞_perp2 = P[:, perps[2]]
        plot_plane!(bbuf, pt, 𝐞_perp1, 𝐞_perp2; gs.halfsize)
        # Function converting Gray{Bool} to color index i
        f_color = x -> RGBA{N0f8}(RED_GREEN_BLUE[i, :]..., N0f8(x == 1))
        # Overlay bbuf on img in the proper color
        map!(BlendLighten, img, img, f_color.(bbuf)) # Consider TODO Not here, move this out if speed is important
    end
    img
end
=#

#######################
# Plot glyphs
#######################

"""
    plot_glyphs(z, pts, gs::GlyphSpec)
"""
function plot_glyphs(z::Matrix{<:AbstractFloat}, pts, gs::GlyphSpec)
    # Allocate an empty color image (since user didn't supply one)
    img = zeros(RGBA{N0f8}, size(z)...)
    # Modify the image
    plot_glyphs!(img, z, pts, gs)
end
function plot_glyphs!(img::Matrix{<:Colorant{N0f8, 4}}, z::Matrix{<:AbstractFloat}, pts, gs::GlyphSpec) 
    # This layer in the call hierarchy: img-> cov -> img
    cov = coverage_fitting_image(img, gs)
    # Modify cover buffer(s)
    plot_glyphs!(cov, z, pts, gs)
    # Cover ('unlimited coverage') maped to colors
    apply_color_by_coverage!(img, cov, gs)
    img
end
function plot_glyphs!(img::Matrix{<:Colorant{N0f8, 4}}, b::U, pts, gs::GlyphSpec) where U <:Union{BidirectionOnGrid, DirectionOnGrid}
    # This layer in the call hierarchy: img-> cov -> img
    cov = coverage_fitting_image(img, gs)
    # Modify cover buffer(s)
    plot_glyphs!(cov, b, pts, gs)
    # Cover ('unlimited coverage') maped to colors
    apply_color_by_coverage!(img, cov, gs)
    img
end
function plot_glyphs!(cov::T, z::Matrix{<:AbstractFloat}, pts, gs::GlyphSpec) where T<:Union{Matrix{Float32}, Vector{Matrix{Float32}}}
    # This makes a functor for z (containing the z values) and passes it on.
    throw("No default functor implemented for $(typeof(gs))")
end
function plot_glyphs!(cov::T, z::Matrix{<:AbstractFloat}, pts, gs::GSTensor) where T<:Union{Matrix{Float32}, Vector{Matrix{Float32}}}
    # This makes a default functor for z (containing the z values) and passes it on.
    b = BidirectionOnGrid(𝐊!, z)
    plot_glyphs!(cov, b, pts, gs)
end
function plot_glyphs!(cov::Matrix{Float32}, z::Matrix{<:AbstractFloat}, pts, gs::GSVector)
    # This makes a default functor for z (containing the z values) and passes it on.
    d = DirectionOnGrid(𝐧ₚ!, z)
    plot_glyphs!(cov, d, pts, gs)
end
function plot_glyphs!(cov::T, b::U, pts, gs::GlyphSpec) where 
    {T<:Union{Matrix{Float32}, Vector{Matrix{Float32}}}, U <:Union{BidirectionOnGrid, DirectionOnGrid}}
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
    plot_glyphs_given_values(pts, values, gs::GlyphSpec)
"""
function plot_glyphs_given_values(pts, values, gs::GlyphSpec,)
    # Allocate an empty color image (since user didn't supply one)
    img = zeros(RGBA{N0f8}, size(z)...)
    # Modify the image
    plot_glyphs_given_values!(img, pts, values, gs)
end

"""
    plot_glyphs_given_values!(img, pts, values, gs::U) where U<:GlyphSpec
    plot_glyphs_given_values!(cov::T, pts, values, gs::U) where 
    {T<:Union{Matrix{Float32}, Vector{Matrix{Float32}}}, U<:GlyphSpec}
"""
function plot_glyphs_given_values!(img, pts, values, gs::U) where U<:GlyphSpec
    # This layer in the call hierarchy: img-> cov -> img
    cov = coverage_fitting_image(img, gs)
    # Modify cover buffer(s)
    plot_glyphs_given_values!(cov, pts, values, gs)
    # Cover ('unlimited coverage') maped to colors
    apply_color_by_coverage!(img, cov, gs.color)
    img
end
function plot_glyphs_given_values!(cov::T, pts, values, gs::U) where 
    {T<:Union{Matrix{Float32}, Vector{Matrix{Float32}}}, U<:GlyphSpec}
    @assert length(values) == length(pts)
    # This plots a glyph at each pt without storing the value
    for (pt, value) in zip(pts, values)
        plot_glyph_given_value!(cov, pt, value, gs)
    end
    cov
end














#=

# This method is called by 'pack_curvature_glyphs`, since 
# values are already calculated.
function plot_curvature_glyphs!(img, gs::GSTensor, pts, vK)
    throw("dead")
    if length(gs.directions) == 1
        # Coverage buffer
        cov = zeros(Float32, size(img)...)
    else
        # Coverage buffers
        cov = [zeros(Float32, size(img)...), zeros(Float32, size(img)...)]
    end
    # Modify cover(s)
    plot_curvature_glyphs!(cov, gs, pts, vK)
    # Apply to img in specified color(s)
    apply_color_by_coverage!(img, cov, gs)
    img
end
function plot_curvature_glyphs!(img, gs::GSVector, pts, vK)
    throw("dead")
    cov = zeros(Float32, size(img)...)
    # Modify cover(s)
    plot_curvature_glyphs!(cov, gs, pts, vK)
    # Apply to img in specified color(s)
    apply_color_by_coverage!(img, cov, gs)
    img
end



# This is the two-colour method
#=
 [1] plot_curvature_glyphs!(cov::Vector{…}, gs::GSTensor, pts::Vector{…}, checked_values::Vector{…})
   @ BitmapMapsExtras C:\Users\f\.julia\dev\BitmapMapsExtras\src\calculate_and_draw_glyphs.jl:128
 [2] plot_curvature_glyphs!(img::Matrix{…}, gs::GSTensor, pts::Vector{…}, vK::Vector{…})
   @ BitmapMapsExtras C:\Users\f\.julia\dev\BitmapMapsExtras\src\calculate_and_draw_glyphs.jl:122
 [3] pack_curvature_glyphs!(img::Matrix{…}, z::Matrix{…}, gs::GSTensor; scatterdist::Float32, seed::MersenneTwister)
   @ Main c:\Users\f\.julia\dev\BitmapMapsExtras\test\t_calculate_and_draw_glyphs.jl:392
 [4] pack_curvature_glyphs(z::Matrix{Float64}, gs::GSTensor; scatterdist::Float32, seed::MersenneTwister)
   @ Main c:\Users\f\.julia\dev\BitmapMapsExtras\test\t_calculate_and_draw_glyphs.jl:397
 [5] pack_curvature_glyphs(z::Matrix{Float64}, gs::GSTensor)
=#
function plot_curvature_glyphs!(cov::T, gs::GSTensor, pts, checked_values) where T<:Vector{Matrix{Float32}}
    throw("dead")
    @assert length(pts) == length(checked_values)
    # Plot curvature glyphs for internal points one at a time
    for (pt, K) in zip(pts, checked_values)
        # Scale and plot the single glyph
        v1 = gs.multip .* view(K, :, 1)
        draw_bidirectional_quantity_glyph!(cov[1], pt, v1, gs.strength)
        v2 = gs.multip .* view(K, :, 2)
        draw_bidirectional_quantity_glyph!(cov[2], pt, v2, gs.strength)
    end
    cov
end

# This is the one-colour method
#=
 [1] plot_curvature_glyphs!(cov::Vector{…}, gs::GSTensor, pts::Vector{…}, checked_values::Vector{…})
   @ BitmapMapsExtras C:\Users\f\.julia\dev\BitmapMapsExtras\src\calculate_and_draw_glyphs.jl:128
 [2] plot_curvature_glyphs!(img::Matrix{…}, gs::GSTensor, pts::Vector{…}, vK::Vector{…})
   @ BitmapMapsExtras C:\Users\f\.julia\dev\BitmapMapsExtras\src\calculate_and_draw_glyphs.jl:122
 [3] pack_curvature_glyphs!(img::Matrix{…}, z::Matrix{…}, gs::GSTensor; scatterdist::Float32, seed::MersenneTwister)
   @ Main c:\Users\f\.julia\dev\BitmapMapsExtras\test\t_calculate_and_draw_glyphs.jl:392
 [4] pack_curvature_glyphs(z::Matrix{Float64}, gs::GSTensor; scatterdist::Float32, seed::MersenneTwister)
   @ Main c:\Users\f\.julia\dev\BitmapMapsExtras\test\t_calculate_and_draw_glyphs.jl:397
 [5] pack_curvature_glyphs(z::Matrix{Float64}, gs::GSTensor)
=#
function plot_curvature_glyphs!(cov::Matrix{<:AbstractFloat}, gs::GSTensor, pts, checked_values)
    throw("dead")
    # Plot curvature glyphs for internal points one at a time
    for (pt, K) in zip(pts, checked_values)
        # Scale and plot the single glyph
        v = gs.multip .* view(K, :, first(gsdirections))
        draw_bidirectional_quantity_glyph!(cov, pt, v, gs.strength)
    end
    cov
end


####################################
# Plot normal unit vector projection
####################################


"""
    oldlot_glyphs(z, pts; gs = GSVector())
"""
function oldoldoldlot_glyphs(z, pts; gs = GSVector())
    throw("dead")
    # Allocate an empty color image (since user didn't supply one)
    img = zeros(RGBA{N0f8}, size(z)...)
    # Modify the image
    oldlot_glyphs!(img, z, pts; gs)
end


"""
    oldlot_glyphs!(img, z, pts; gs = GSVector())
    oldlot_glyphs!(cov::Matrix{Float32}, z, pts; gs = GSVector())
    oldlot_glyphs!(cov::Matrix{Float32}, z, pts, Ri, Ω, v, gs)
"""
function oldlot_glyphs!(img, z, pts; gs = GSVector())
    throw("dead")
    # Coverage buffer
    cov = zeros(Float32, size(img)...)
    # Modify coverage 
    oldlot_glyphs!(cov, z, pts, gs)
    # To color
    apply_color_by_coverage!(img, cov, gs.color)
    img
end
function oldlot_glyphs!(cov::Matrix{Float32}, z, pts, gs::GSVector)
    throw("dead. We should have the function passed.0
    + ")
    # Prepare
    Ri, Ω, v, _, _, _, _, _, _ = allocations_curvature(CartesianIndices(z))
    # Plot projected vector glyphs for internal points one at a time
    for pt in filter(pt -> pt ∈ Ri, sort(vec(pts)))
        # Find 𝐧ₚ in-place
        𝐧ₚ!(v, view(z, Ω .+ pt))
        # Scale and plot the single glyph
        plot_glyph_given_value!(cov, pt, v, gs)
    end
    cov
end
=#

# These methods are placed here because it uses an 'internal' type,
# whereas draw_direct is more general and might be moved to 
# a separate package.
"""
    apply_color_by_coverage!(img, cov, gs::GlyphSpec)
    apply_color_by_coverage!(img, cov, gs::GSTensor)
"""
function apply_color_by_coverage!(img, cov, gs::GlyphSpec)
    apply_color_by_coverage!(img, cov, gs.color)
end
function apply_color_by_coverage!(img, cov, gs::GSTensor)
    if length(gs.directions) == 1
        color = first(gs.directions) == 1 ? gs.color1 : gs.color2
        apply_color_by_coverage!(img, cov, color)
    else
        apply_color_by_coverage!(img, cov[1], gs.color1)
        apply_color_by_coverage!(img, cov[2], gs.color2)
    end
    img
end


"""
    coverage_fitting_image(img, gs::GlyphSpec)
    coverage_fitting_image(img, gs::GSTensor)
"""
coverage_fitting_image(img, gs::GlyphSpec) = zeros(Float32, size(img)...)
function coverage_fitting_image(img, gs::GSTensor)
    if length(gs.directions) == 1
        # Coverage buffer
        zeros(Float32, size(img)...)
    else
        # Coverage buffers (potentially for two different colors)
        [zeros(Float32, size(img)...), zeros(Float32, size(img)...)]
    end
end
