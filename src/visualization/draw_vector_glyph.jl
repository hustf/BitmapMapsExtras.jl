# Drawing functionality,
# builds on `draw_direct.jl` and `BitmapMaps.jl/mark_utils.jl`.
# Contains `plot_vector_glyph!`

"""
    plot_vector_glyph!(cov, pt, v, gs::GSVector)
"""
function plot_vector_glyph!(cov, pt, v, gs::GSVector)
    #if norm(v) > 0
    #        @show v gs
    #        throw("ok")
    #end
    if is_in_limits(gs, v)
        Δj = Int(round(gs.multip * v[1]))
        Δi = -Int(round(gs.multip * v[2]))
        draw_vector!(cov, pt, Δi, Δj, gs.strength)
    else
        if gs.dashsize > 0
            # A dash instead of the out-of-limits glyph
            spray!(cov, pt, gs.dashsize, gs.strength)
        end
    end
    cov
end
