"""
    plot_glyph_given_value!(cov, pt, v, gs::GSVector)

pt is in CartesianIndex (i is image row)
v is a 2d vector in screen coordinates (x, y)
"""
function plot_glyph_given_value!(cov, pt, v, gs::GSVector)
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


function plot_glyph_given_value!(cov, pt, K, gs::GSTensor)
    throw("this here")
    if is_in_limits(gs, K)
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