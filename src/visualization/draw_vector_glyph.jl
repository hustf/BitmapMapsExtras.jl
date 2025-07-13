# Drawing functionality,
# builds on `draw_direct.jl` and `BitmapMaps.jl/mark_utils.jl`.
# Contains `plot_vector!`, which should have similar scaling and limitation 
# functionality as 'plot_principal_directions_glyph!` (for which the complexity
# is more justified than here).
# One callee: `draw_vector!`

"""
    plot_vector!(bbuf, pt, f_is_within_limits, dashsize, v)

This checks if the glyph is within limits. v is a (scaled for graphics) 2d vector,
so `f_is_within_limits` would typically be this simple function:

    v -> minglyph ≤ norm(v) ≤ maxglyph

In keeping with similar plot functions, pt is v in 'image-space',
where [1,2] is upper left, one pixel to the right.

'v' is, awkwardly, in (x, y) coordinates. 
"""
function plot_vector!(bbuf, pt, f_is_within_limits, dashsize, v)
    if f_is_within_limits(v)
        Δj = Int(round(v[1]))
        Δi = -Int(round(v[2]))
        draw_vector!(bbuf, pt, Δi, Δj)
    else
        # A dash instead of the out-of-limits vector glyph
        spray_neighbors!(bbuf, CartesianIndices(bbuf), pt, dashsize)
    end
    bbuf
end
