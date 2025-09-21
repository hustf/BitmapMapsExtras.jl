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


#######################
# Plot curvature glyphs
#######################

"""
    plot_curvature_glyphs(z, pts; gs = GSTensor())
"""
function plot_curvature_glyphs(z, pts; gs = GSTensor())
    # Allocate an empty color image (since user didn't supply one)
    img = zeros(RGBA{N0f8}, size(z)...)
    # Modify the image
    plot_curvature_glyphs!(img, z, pts; gs)
end

"""
    plot_curvature_glyphs!(img, z, pts; gs = GSTensor())
    plot_curvature_glyphs!(cov::Matrix{Float32}, z, pts, directions, gs::GSTensor)
    plot_curvature_glyphs!(img, gs::GSTensor, pts, values)
"""
function plot_curvature_glyphs!(img, z, pts; gs = GSTensor())
    if length(gs.directions) == 1
        # Coverage buffer
        cov = zeros(Float32, size(img)...)
    else
        # Coverage buffers
        cov = [zeros(Float32, size(img)...), zeros(Float32, size(img)...)]
    end
    # Modify cover buffer(s)
    plot_curvature_glyphs!(cov, z, pts, gs)
    apply_color_by_coverage!(img, cov, gs)
    img
end
function  plot_curvature_glyphs!(cov::T, z, pts, gs::GSTensor) where T<:Union{Matrix{Float32}, Vector{Matrix{Float32}}}
    # Prepare
    Ri, Ω, v, P, K, vα, vκ, vβ, lpc = allocations_curvature(CartesianIndices(z))
    # Plot curvature glyphs for internal points one at a time
    for pt in filter(pt -> pt ∈ Ri, sort(vec(pts)))
        # Find P in-place
        tangent_basis!(P, v, view(z, Ω .+ pt))
        # Update K etc. 
        principal_curvature_components!(K, vα, vβ, vκ, P, view(z, Ω .+ pt), VΦ, lpc)
        # Scale and plot the single glyph
        plot_principal_directions_glyph!(cov, pt, K, gs)
    end
    cov
end
# This method is called by 'pack_curvature_glyphs`, since 
# values are already calculated.
function plot_curvature_glyphs!(img, gs::GSTensor, pts, vK)
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
function plot_curvature_glyphs!(cov::Matrix{Float32}, gs::GSTensor, pts, checked_values)
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
    plot_𝐧ₚ__glyphs(z, pts; gs = GSVector())
"""
function plot_𝐧ₚ_glyphs(z, pts; gs = GSVector())
    # Allocate an empty color image (since user didn't supply one)
    img = zeros(RGBA{N0f8}, size(z)...)
    # Modify the image
    plot_𝐧ₚ_glyphs!(img, z, pts; gs)
end


"""
    plot_𝐧ₚ_glyphs!(img, z, pts; gs = GSVector())
    plot_𝐧ₚ_glyphs!(cov::Matrix{Float32}, z, pts; gs = GSVector())
    plot_𝐧ₚ_glyphs!(cov::Matrix{Float32}, z, pts, Ri, Ω, v, gs)
"""
function plot_𝐧ₚ_glyphs!(img, z, pts; gs = GSVector())
    # Coverage buffer
    cov = zeros(Float32, size(img)...)
    # Modify coverage 
    plot_𝐧ₚ_glyphs!(cov, z, pts, gs)
    # To color
    apply_color_by_coverage!(img, cov, gs.color)
    img
end
function plot_𝐧ₚ_glyphs!(cov::Matrix{Float32}, z, pts, gs::GSVector)
    # Prepare
    Ri, Ω, v, _, _, _, _, _, _ = allocations_curvature(CartesianIndices(z))
    # Plot projected vector glyphs for internal points one at a time
    for pt in filter(pt -> pt ∈ Ri, sort(vec(pts)))
        # Find 𝐧ₚ in-place
        𝐧ₚ!(v, view(z, Ω .+ pt))
        # Scale and plot the single glyph
        plot_vector_glyph!(cov, pt, v, gs)
    end
    cov
end


# This method is placed here because it uses an 'internal' type,
# whereas draw_direct is more general and might be moved to 
# a separate package.
"""
    apply_color_by_coverage!(img, cov, gs::GSTensor)
"""
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