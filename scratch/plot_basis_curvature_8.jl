# Functions with callees in
# `utilites_graphic_6.jl`
# `utilties_tangent_5.jl`
# `utilties_curvature_7.jl`

"""
Low-effort definition of distinguishable plane colors in roughly similar luminance.
"""
const RED_GREEN_BLUE = SMatrix{3, 3, N0f8}([ 0.95  0.0  0.0
                                        0.0  0.5  0.0
                                        0.0  0.0  0.9])

"""
Low-effort definition of curvature glyph color, used where overlain a colorful picture.
"""
const COLOR_CURVGLYPH = RGB{N0f8}(0.85, 0.5, 0.9)
##############################
# Plot tangent basis functions
##############################

function plot_tangent_basis_glyphs(z, pts; halfsize = 30)
    # Allocate an empty colorful image
    img = zeros(RGBA{N0f8}, size(z)...)
    # Modify the image
    plot_tangent_basis_glyphs!(img, z, pts; halfsize)
end

function plot_tangent_basis_glyphs!(img, z, pts; halfsize = 30)
    Ri, Ω, v, P, _, _, _, _, _ = _allocations_curvature(CartesianIndices(z), [])
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




function plot_plane!(img, pt, e_perp1, e_perp2; halfsize)
    # We're plotting a 2d projection of the plane that is
    # - centered at the origin
    # - has basis directions along e_perp1 and e_perp2, in standard basis
    # Define four corners in 3D space
    p1 = halfsize * (-e_perp1 - e_perp2)
    p2 = halfsize * ( e_perp1 - e_perp2)
    p3 = halfsize * ( e_perp1 + e_perp2)
    p4 = halfsize * (-e_perp1 + e_perp2)
    # Project to 2D by dropping z-component
    # Translate to pt
    # Compensate for image-y down
    p1_2d = pt + CartesianIndex(-Int(round(p1[2])), Int(round(p1[1])))
    p2_2d = pt + CartesianIndex(-Int(round(p2[2])), Int(round(p2[1])))
    p3_2d = pt + CartesianIndex(-Int(round(p3[2])), Int(round(p3[1])))
    p4_2d = pt + CartesianIndex(-Int(round(p4[2])), Int(round(p4[1])))
    # Draw lines to form the quadrilateral
    line!(img, p1_2d, p2_2d)
    line!(img, p2_2d, p3_2d)
    line!(img, p3_2d, p4_2d)
    line!(img, p4_2d, p1_2d)
    img
end

##########################
# Plot curvature functions
##########################
function plot_curvature_glyphs(z, pts; directions = 1:2, multglyph = 50, minglyph = -50, maxglyph = 50)
    # Allocate an empty color image (since user didn't supply one)
    img = zeros(RGBA{N0f8}, size(z)...)
    # Modify the image
    plot_curvature_glyphs!(img, z, pts; directions, multglyph, minglyph, maxglyph)
end
function plot_curvature_glyphs!(img, z, pts; directions = 1:2, multglyph = 50, minglyph = -50, maxglyph = 50)
    # Black-white buffer
    bbuf = Array{GrayA{N0f8}}(falses( size(img)...))
    # Modify the image
    plot_curvature_glyphs!(bbuf, img, z, pts, directions; multglyph, minglyph, maxglyph)
end

function  plot_curvature_glyphs!(bbuf, img, z, pts, directions; multglyph = 50, maxglyph = 50, minglyph = -50)
    # Prepare
    Ri, Ω, v, P, K, vα, vκ, vβ, f_is_within_limits = _allocations_curvature(CartesianIndices(z), directions; maxglyph, minglyph)
    # Size of dash placed where the glyph size is outside limits
    dashsize = maxglyph ÷ 10
    # Plot curvature glyps for internal points one at a time
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


function plot_principal_directions_glyph!(bbuf, pt, directions, f_is_within_limits, dashsize, K)
    # Extract primary and / or secondary principal components
    dK =  K[:, directions]
    if f_is_within_limits(dK)
        draw_bidirectional_quantity_glyph!(bbuf, pt, dK)
    else
        # A dash instead of the out-of-limits curvature glyph
        color_neighbors!(bbuf, CartesianIndices(bbuf), pt, dashsize)
    end
    bbuf
end









