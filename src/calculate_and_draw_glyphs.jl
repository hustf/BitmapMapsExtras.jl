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
    plot_tangent_basis_glyphs(z, pts; halfsize = 30)
"""
function plot_tangent_basis_glyphs(z, pts; halfsize = 30)
    # Allocate an empty colorful image
    img = zeros(RGBA{N0f8}, size(z)...)
    # Modify the image
    plot_tangent_basis_glyphs!(img, z, pts; halfsize)
end

"""
    plot_tangent_basis_glyphs!(img, z, pts; halfsize = 30)
"""
function plot_tangent_basis_glyphs!(img, z, pts; halfsize = 30)
    Ri, Ω, v, P, _, _, _, _, _ = allocations_curvature(CartesianIndices(z), [])
    # Black-white buffer
    bbuf = Array{Gray{Bool}}(falses( size(img)...))
    # Plot tangent basis for internal points one at a time
    for pt in filter(pt -> pt ∈ Ri, sort(vec(pts)))
        # Find P in-place
        tangent_basis!(P, v, view(z, Ω .+ pt))
        # Plot the single glyph
        plot_orthonormal_basis_glyph!(img, bbuf, pt, P, halfsize)
    end
    img
end

"""
    plot_orthonormal_basis_glyph!(img, bbuf::Array{Gray{Bool}}, pt, P, halfsize)
"""
function plot_orthonormal_basis_glyph!(img, bbuf::Array{Gray{Bool}}, pt, P, halfsize)
    @assert size(img) == size(bbuf) 
    for i in 1:3
        # Blank bw-buffer
        fill!(bbuf, Gray{Bool}(false))
        perps = setdiff([1,2,3], i)
        𝐞_perp1 = P[:, perps[1]]
        𝐞_perp2 = P[:, perps[2]]
        plot_plane!(bbuf, pt, 𝐞_perp1, 𝐞_perp2; halfsize)
        # Function converting Gray{Bool} to color index i
        f_color = x -> RGBA{N0f8}(RED_GREEN_BLUE[i, :]..., N0f8(x == 1))
        # Overlay bbuf on img in the proper color
        map!(BlendLighten, img, img, f_color.(bbuf))
    end
    img
end


#######################
# Plot curvature glyphs
#######################

"""
    plot_curvature_glyphs(z, pts; directions = 1:2, multglyph = 50, minglyph = -50, maxglyph = 50)
"""
function plot_curvature_glyphs(z, pts; directions = 1:2, multglyph = 50, minglyph = -50, maxglyph = 50)
    # Allocate an empty color image (since user didn't supply one)
    img = zeros(RGBA{N0f8}, size(z)...)
    # Modify the image
    plot_curvature_glyphs!(img, z, pts; directions, multglyph, minglyph, maxglyph)
end

"""
    plot_curvature_glyphs!(img, z, pts; directions = 1:2, 
        multglyph = 50, minglyph = -50, maxglyph = 50, dashsize = maxglyph ÷ 10)
"""
function plot_curvature_glyphs!(img, z, pts; directions = 1:2, 
    multglyph = 50, minglyph = -50, maxglyph = 50, dashsize = maxglyph ÷ 10)
    # Black-white buffer
    bbuf = Array{GrayA{N0f8}}(falses( size(img)...))
    # Modify the image
    plot_curvature_glyphs!(bbuf, img, z, pts, directions; multglyph, minglyph, maxglyph, dashsize)
end

"""
    plot_curvature_glyphs!(bbuf, img, z, pts, directions; 
        multglyph = 50, maxglyph = 50, minglyph = -50, dashsize = maxglyph ÷ 10)
"""
function  plot_curvature_glyphs!(bbuf, img, z, pts, directions; 
    multglyph = 50, maxglyph = 50, minglyph = -50, dashsize = maxglyph ÷ 10)
    # Prepare
    Ri, Ω, v, P, K, vα, vκ, vβ, f_is_within_limits = allocations_curvature(CartesianIndices(z), directions; maxglyph, minglyph)
    # Plot curvature glyphs for internal points one at a time
    for pt in filter(pt -> pt ∈ Ri, sort(vec(pts)))
        # Find P in-place
        tangent_basis!(P, v, view(z, Ω .+ pt))
        # Update K 
        principal_curvature_components!(K, vα, vβ, vκ, P, view(z, Ω .+ pt), VΦ)
        # Scale and plot the single glyph
        plot_principal_directions_glyph!(bbuf, pt, directions, f_is_within_limits, dashsize, multglyph * K)
    end
    # Function converting GrayA{N0f8} to proper color
    f = x -> RGBA{N0f8}(COLOR_CURVGLYPH.r, COLOR_CURVGLYPH.g, COLOR_CURVGLYPH.b, 
                            x.val)
    # Composite bbuf over img
    # Overlay bbuf on img in the proper color
    map!(BlendLighten, img, img, f.(bbuf))
    img
end

####################################
# Plot normal unit vector projection
####################################

"""
    plot_𝐧ₚ_glyphs(z, pts; 
        multglyph = 50, minglyph = -50, maxglyph = 50, dashsize = maxglyph ÷ 10)
"""
function plot_𝐧ₚ_glyphs(z, pts; 
    multglyph = 50, minglyph = -50, maxglyph = 50, dashsize = maxglyph ÷ 10)
    # Allocate an empty color image (since user didn't supply one)
    img = zeros(RGBA{N0f8}, size(z)...)
    # Modify the image
    plot_𝐧ₚ_glyphs!(img, z, pts; multglyph, minglyph, maxglyph, dashsize)
end

"""
    plot_𝐧ₚ_glyphs!(img, z, pts; 
        multglyph = 50, minglyph = -50, maxglyph = 50, dashsize = maxglyph ÷ 10)
"""
function plot_𝐧ₚ_glyphs!(img, z, pts; 
    multglyph = 50, minglyph = -50, maxglyph = 50, dashsize = maxglyph ÷ 10)
    # Black-white buffer
    bbuf = Array{GrayA{N0f8}}(falses( size(img)...))
    # Modify the image
    plot_𝐧ₚ_glyphs!(bbuf, img, z, pts; multglyph, minglyph, maxglyph, dashsize)
end

"""
    plot_𝐧ₚ_glyphs!(bbuf, img, z, pts; 
        multglyph = 50, maxglyph = 50, minglyph = -50, dashsize = maxglyph ÷ 10)
"""
function plot_𝐧ₚ_glyphs!(bbuf, img, z, pts; 
        multglyph = 50, maxglyph = 50, minglyph = -50, dashsize = maxglyph ÷ 10)
    # Allocate
    Ri, Ω, v, _, _, _, _, _, _ = allocations_curvature(CartesianIndices(z), [])
    # Captures maxglyph and minglyph, limitations on vector length (negative limits
    # are irrelevant here)
    f_is_within_limits = let minglyph = float(minglyph), maxglyph = float(maxglyph)
        v -> minglyph ≤ norm(v) ≤ maxglyph
    end
    # Plot projected vector glyphs for internal points one at a time
    for pt in filter(pt -> pt ∈ Ri, sort(vec(pts)))
        # Find 𝐧ₚ in-place, mutates v.
        # v is in the format (dz/dx, dz/dy)
        𝐧ₚ!(v, view(z, Ω .+ pt))
        # Scale and plot the single glyph
        plot_vector!(bbuf, pt, f_is_within_limits, dashsize, multglyph * v)
    end
    # Function converting GrayA{N0f8} to proper color.
    # CONSIDER TODO: Separate color constant for vectors.
    #     Also, use mutable containers for color definitions.
    f = x -> RGBA{N0f8}(COLOR_CURVGLYPH.r, COLOR_CURVGLYPH.g, COLOR_CURVGLYPH.b, 
                            x.val)
    # Composite bbuf over img
    # Overlay bbuf on img in the proper color
    map!(BlendLighten, img, img, f.(bbuf))
    img
end