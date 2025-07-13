# Drawing functionality,
# builds on `BitmapMaps.jl/mark_utils.jl`.
# Contains `plot_plane!`, called from `calculate_and_draw.jl`

"""
    plot_plane!(img, pt, e_perp1, e_perp2; halfsize)

Black-and-white image modified in-place.
"""
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