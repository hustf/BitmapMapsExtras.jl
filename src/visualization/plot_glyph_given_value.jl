"""
    plot_glyph_given_value!(cov, pt, v, gs::GSVector)
    plot_glyph_given_value!(cov::Matrix{Float32}, pt, K, gs::GSTensor{D, 1}) where D
    plot_glyph_given_value!(cov::Vector{Matrix{Float32}}, pt, K, gs::GSTensor{D12, 2})
    plot_glyph_given_value!(cov, pt, P, gs::GSTangentBasis)

- `cov`: Coverage matrix or matrices. Each element (pixel) accumulate coverage, how fast depends on strength. See `coverage_fitting_image`.
- `pt` is in CartesianIndex (i is image row)
- `v`` is a 2d vector in screen coordinates (x-y up)
- `K` is a TENSORMAP.
- `gs` is a type <: `AbstractGlyphSpec`.

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
function plot_glyph_given_value!(cov::Matrix{Float32}, pt, K, gs::GSTensor{D, 1}) where D
    # Single direction method
    if is_in_limits(gs, K)
        # Scale the single glyph value
        v = gs.multip .* view(K, :, D)
        draw_bidirectional_vector!(cov, pt, v, gs.strength)
    else
        if gs.dashsize > 0
            # A dash instead of the out-of-limits glyph
            spray!(cov, pt, gs.dashsize, gs.strength)
        end
    end
    cov
end
function plot_glyph_given_value!(cov::Vector{Matrix{Float32}}, pt, K, gs::GSTensor{D12, 2})
    # Two directions method
    if is_in_limits(gs, K)
        # Scale the two values to 
        vs1 = gs.multip .* view(K, :, 1)
        vs2 = gs.multip .* view(K, :, 2)
        # Plot the two values in the two directions
        draw_bidirectional_vector!(cov[1], pt, vs1, gs.strength)
        draw_bidirectional_vector!(cov[2], pt, vs2, gs.strength)
    else
        if gs.dashsize > 0
            # A dash instead of the out-of-limits glyph
            spray!(cov[1], pt, gs.dashsize, gs.strength)
            spray!(cov[2], pt, gs.dashsize, gs.strength)
        end
    end
    cov
end
function plot_glyph_given_value!(cov, pt, P, gs::GSTangentBasis)
    for i in 1:3
        perps = setdiff([1,2,3], i)
        𝐞_perp1 = P[:, perps[1]]
        𝐞_perp2 = P[:, perps[2]]
        draw_plane!(cov[i], pt, 𝐞_perp1, 𝐞_perp2; gs.halfsize)
    end
    cov
end