# Drawing functionality,
# builds on `draw_direct.jl` and `BitmapMaps.jl/mark_utils.jl`.
# Contains `plot_vector!`, which should have similar scaling and limitation 
# functionality as 'plot_principal_directions_glyph!` (for which the complexity
# is more justified than here).
# One callee: `draw_vector!`

"""
    plot_vector!(cov::Matrix{Float32}, pt, f_is_within_limits, dashsize, v, strength::Float32)

'v' is, awkwardly, expected to be in (x, y) coordinates. 

`f_is_within_limits` would typically be this simple function:

    v -> ming ≤ norm(v) ≤ maxg

In keeping with similar plot functions, pt is v in 'image-space',
where [1,2] is upper left, one pixel to the right.

"""
function plot_vector!(cov::Matrix{Float32}, pt, f_is_within_limits, dashsize, v, strength::Float32)
    if f_is_within_limits(v)
        Δj = Int(round(v[1]))
        Δi = -Int(round(v[2]))
        draw_vector!(cov, pt, Δi, Δj, strength)
    else
        if dashsize > 0
            # A dash instead of the out-of-limits vector glyph
            spray!(cov, pt, Float32(dashsize), strength)
        end
    end
    cov
end
