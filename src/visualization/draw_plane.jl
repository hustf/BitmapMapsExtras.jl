# Drawing functionality,
# builds on `BitmapMaps.jl/mark_utils.jl`.
# Contains `draw_plane!`, called from `calculate_and_draw.jl`

"""
    draw_plane!(img, pt, e_perp1, e_perp2; halfsize)

Black-and-white image (or 'coverage') modified in-place.
"""
function draw_plane!(img, pt, e_perp1, e_perp2; halfsize)
    # We're plotting a 2d projection of the plane that is
    # - centered at the origin
    # - has basis directions along e_perp1 and e_perp2, in standard basis
    # Define four corners in 3D space
    p1 = halfsize * (-e_perp1 - e_perp2)
    p2 = halfsize * ( e_perp1 - e_perp2)
    p3 = halfsize * ( e_perp1 + e_perp2)
    p4 = halfsize * (-e_perp1 + e_perp2)
    # Define mid-way points in 3D space
    m12 = halfsize * (- e_perp2)
    m23 = halfsize * e_perp1
    m34 = halfsize * e_perp2
    m41 = halfsize * (- e_perp1)
    # Project to 2D by dropping z-component
    # Translate to pt
    # Compensate for image-y down
    pt1   = pt + CartesianIndex(-Int(round(p1[2])), Int(round(p1[1])))
    pt2   = pt + CartesianIndex(-Int(round(p2[2])), Int(round(p2[1])))
    pt3   = pt + CartesianIndex(-Int(round(p3[2])), Int(round(p3[1])))
    pt4   = pt + CartesianIndex(-Int(round(p4[2])), Int(round(p4[1])))
    mpt12 = pt + CartesianIndex(-Int(round(m12[2])), Int(round(m12[1])))
    mpt23 = pt + CartesianIndex(-Int(round(m23[2])), Int(round(m23[1])))
    mpt34 = pt + CartesianIndex(-Int(round(m34[2])), Int(round(m34[1])))
    mpt41 = pt + CartesianIndex(-Int(round(m41[2])), Int(round(m41[1])))
    # Draw lines to form the quadrilateral
    line!(img, pt1, mpt12)
    line!(img, pt2, mpt23)
    line!(img, pt3, mpt34, 2.0)
    line!(img, pt4, mpt41)
    line!(img, mpt12, pt2)
    line!(img, mpt23, pt3, 2.0)
    line!(img, mpt34, pt4)
    line!(img, mpt41, pt1)
    img
end

