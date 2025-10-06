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


# This is the method for one principal direction
function plot_glyph_given_value!(cov::Matrix{Float32}, pt, K, gs::GSTensor)
    if is_in_limits(gs, K)
        # We know that just one direction is wanted since 'cov' contains one matrix
        dirno = first(gs.directions)
        # Scale the single glyph value
        v = gs.multip .* view(K, :, dirno)
        draw_bidirectional_vector_glyph_given_value!(cov, pt, v, gs.strength)
    else
        if gs.dashsize > 0
            # A dash instead of the out-of-limits glyph
            spray!(cov, pt, gs.dashsize, gs.strength)
        end
    end
    cov
end
# This is the method for two principal directions
function plot_glyph_given_value!(cov::Vector{Matrix{Float32}}, pt, K, gs::GSTensor)
    if is_in_limits(gs, K)
        # We know that two directions are wanted since 'cov' contains two matrices
        # Note that this is not possible: gs.directions == 2:-1:1
        @assert gs.directions == 1:2
        # Scale the two glyph values
        v1 = gs.multip .* view(K, :, 1)
        v2 = gs.multip .* view(K, :, 2)
        # Plot the two double arrows
        draw_bidirectional_vector_glyph_given_value!(cov[1], pt, v1, gs.strength)
        draw_bidirectional_vector_glyph_given_value!(cov[2], pt, v2, gs.strength)
    else
        if gs.dashsize > 0
            # A dash instead of the out-of-limits glyph
            spray!(cov[1], pt, gs.dashsize, gs.strength)
            spray!(cov[2], pt, gs.dashsize, gs.strength)
        end
    end
    cov
end