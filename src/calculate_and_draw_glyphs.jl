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
    Ri, Ω, v, P, _, _, _, _ = allocations_curvature(CartesianIndices(z))
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
        map!(BlendLighten, img, img, f_color.(bbuf)) # TODO Not here, move this out
    end
    img
end


#######################
# Plot curvature glyphs
#######################

"""
    plot_curvature_glyphs(z, pts; gs = GSTensor())
    plot_curvature_glyphs!(img, z, pts; gs = GSTensor())
    plot_curvature_glyphs!(cov::Matrix{Float32}, z, pts, directions, gs::GSTensor)
"""
function plot_curvature_glyphs(z, pts; gs = GSTensor())
    # Allocate an empty color image (since user didn't supply one)
    img = zeros(RGBA{N0f8}, size(z)...)
    # Modify the image
    plot_curvature_glyphs!(img, z, pts; gs)
end
function plot_curvature_glyphs!(img, z, pts; gs = GSTensor())
    # Coverage buffer
    cov = zeros(Float32, size(img)...)
    # Modify cover
    plot_curvature_glyphs!(cov, z, pts, gs)
    # Apply color by coverage
    apply_color_by_coverage!(img, cov, gs.color)
end
function  plot_curvature_glyphs!(cov::Matrix{Float32}, z, pts, gs::GSTensor)
    # Prepare
    Ri, Ω, v, P, K, vα, vκ, vβ = allocations_curvature(CartesianIndices(z))
    # Plot curvature glyphs for internal points one at a time
    for pt in filter(pt -> pt ∈ Ri, sort(vec(pts)))
        # Find P in-place
        tangent_basis!(P, v, view(z, Ω .+ pt))
        # Update K etc. 
        principal_curvature_components!(K, vα, vβ, vκ, P, view(z, Ω .+ pt), VΦ)
        # Scale and plot the single glyph
        plot_principal_directions_glyph!(cov, pt, K, gs)
    end
    cov
end

####################################
# Plot normal unit vector projection
####################################

"""
    plot_𝐧ₚ_glyphs!(img, z, pts; gs = GSVector())
    plot_𝐧ₚ_glyphs!(cov::Matrix{Float32}, z, pts; gs = GSVector())
    plot_𝐧ₚ_glyphs!(cov::Matrix{Float32}, z, pts, Ri, Ω, v, gs)
"""
function plot_𝐧ₚ_glyphs!(img, z, pts; gs = GSVector())
    # Coverage buffer
    cov = zeros(Float32, size(img)...)
    # Modify coverage 
    plot_𝐧ₚ_glyphs!(cov, z, pts; gs)
    # To color
    apply_color_by_coverage!(img, cov, gs.color)
    img
end
function plot_𝐧ₚ_glyphs!(cov::Matrix{Float32}, z, pts; gs = GSVector())
    # Allocate
    Ri, Ω, v, _, _, _, _, _ = allocations_curvature(CartesianIndices(z))
    # Plot projected vector glyphs for internal points one at a time
    plot_𝐧ₚ_glyphs!(cov, z, pts, Ri, Ω, v, gs)
    cov
end
function plot_𝐧ₚ_glyphs!(cov::Matrix{Float32}, z, pts, Ri, Ω, v, gs)
    f = v -> is_in_limits(gs, v)
    # Plot projected vector glyphs for internal points one at a time
    for pt in filter(pt -> pt ∈ Ri, sort(vec(pts)))
        # Find 𝐧ₚ in-place, mutates v.
        # v is in the format (dz/dx, dz/dy)
        𝐧ₚ!(v, view(z, Ω .+ pt))
        # Scale and plot the single glyph
        plot_vector!(cov, pt, f, gs.dashsize, gs.multip * v, gs.strength)
    end
    cov
end
