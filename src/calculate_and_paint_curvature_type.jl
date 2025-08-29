# We can simplify the quantified curvature K a bit by
# categorizing each of the principal direction to +, 0 and -.
# Like we do here, see `curvature_type`. The categorization
# would be interesting for similar bidirectional quantities 
# such as elongation or stress as well.
# Here is also an effective graphical function, `paint_curvature_types`.
#
# Relies on constant PALETTE_GRGB

"""
    curvature_type(K, maxcurv_flat)

Classifies surface curvature based on principal curvatures.

Args:
    K (2x2 Matrix): Principal curvature vectors in columns.
    maxcurv_flat (Float): Threshold for flat curvature.

Returns:
    Int: 1 (flat), 2 (convex), 3 (concave), 4 (saddle).
"""
function curvature_type(K, maxcurv_flat)
    # Compute norms and signs of principal curvatures
    s1 = norm(K[:, 1])
    s2 = norm(K[:, 2])
    pos1 = is_bidirec_vect_positive(K[:, 1])
    pos2 = is_bidirec_vect_positive(K[:, 2])
    # Flat (enough) in both principal directions
    if s1 < maxcurv_flat && s2 < maxcurv_flat     
        return 1
    end
    # Single flat direction: curvature depends on the non-flat direction's sign
    if s1 < maxcurv_flat
        return pos2 ? 2 : 3 # Convex or concave
    elseif s2 < maxcurv_flat
        return pos1 ? 2 : 3 # Convex or concave
    end
    # No flat principal direction: check signs for convex, concave, or saddle
    if pos1 && pos2
        return 2  # Convex
    elseif !pos1 && !pos2
        return 3  # Concave
    else
        return 4  # Saddle
    end
end

# Graphical output 

"""
paint_curvature_types!(img, z, pts;  
    maxcurv_flat = 0.00001, mode::BlendMode = BlendHue, palette = PALETTE_GRGB)

Overlays a color-map onto the image img.

See callees `color_point_by_curvature_type!` and  `curvature_type`!
"""
function paint_curvature_types!(img, z, pts;  
    maxcurv_flat = 0.00001, mode::BlendMode = BlendHue, palette = PALETTE_GRGB)
    #
    buf = zeros(RGB{N0f8}, size(img)...)
    _paint_curvature_types!(buf, z, pts, maxcurv_flat, palette)
    # Composite buf over img in-place.
    # mode = BlendHue: The result is a color with the hue of the source color and the saturation and luminosity of the backdrop color.
    #      --> overwrite img with buf
    @. img = blend(img, buf, mode = mode)
    img
end 

function _paint_curvature_types!(buf, z, pts, maxcurv_flat, palette)
    # Prepare
    Ri, Ω, v, P, K, vα, vκ, vβ, f_is_within_limits = allocations_curvature(CartesianIndices(z), [1, 2])
    # Color points one at a time, hopefully fast
    for pt in filter(pt -> pt ∈ Ri, sort(vec(pts)))
        color_point_by_curvature_type!(buf, z, pt, Ω, v, P, K, vα, vκ, vβ, maxcurv_flat, palette)
    end
    buf
end


"""
    color_point_by_curvature_type!(buf, z, pt, Ω, v, P, K, vα, vκ, vβ, maxcurv_flat, palette)
    color_point_by_curvature_type!(buf, pt,  K, maxcurv_flat, palette)

K is a 2x2 matrix with two bidirectional vectors, returned from `components_matrix!`.

If `palette` is gray, red, green, blue, such as:

```
palette = RGB(0.467, 0.467, 0.467), RGB(0.957, 0.0, 0.078), RGB(0.0, 0.549, 0.0), RGB(0.0, 0.443, 1.0)]

```
...then the CartesianIndex point pt is colored as:

 1: Gray, flat (enough)
 2: Red, convex
 3: Green, concave
 4: Blue, convex-concave, hyperbolic
"""
function color_point_by_curvature_type!(buf, z, pt, Ω, v, P, K, vα, vκ, vβ, maxcurv_flat, palette)
    # Update K etc. 
    principal_curvature_components!(K, vα, vβ, vκ, P, view(z, Ω .+ pt), VΦ)
    # Color by curvature 
    color_point_by_curvature_type!(buf, pt,  K, maxcurv_flat, palette)
    buf
end
function color_point_by_curvature_type!(buf, pt,  K, maxcurv_flat, palette)
    i = curvature_type(K, maxcurv_flat)
    buf[pt] = palette[i]
end


