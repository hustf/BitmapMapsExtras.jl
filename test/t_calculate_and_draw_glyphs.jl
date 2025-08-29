using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: plot_tangent_basis_glyphs, plot_tangent_basis_glyphs!
using BitmapMapsExtras: plot_curvature_glyphs, plot_curvature_glyphs!
using BitmapMapsExtras: plot_𝐧ₚ_glyphs, plot_𝐧ₚ_glyphs!, PALETTE_GRGB
using BitmapMapsExtras: BidirectionOnGrid, 𝐊!, Domain, TENSORMAP
import StatsBase
using StatsBase: Weights, sample, norm
import Random
using Random: MersenneTwister
import ImageContrastAdjustment
using ImageContrastAdjustment: LinearStretching, ContrastStretching, adjust_histogram!
import ImageMorphology
using ImageMorphology: erode!

include("common.jl")

# Tangent basis
z = z_paraboloid()
img = background(z)
plot_tangent_basis_glyphs!(img, z, [CartesianIndex((315, 215))])

########################################################
# Format for 'args', argument to 
# 'fcall', 'grid_fcall' and 'grid_fcall_with_background'
########################################################
# Case 1-3: Already in canonical form (Vector of Tuples)
args = Vector{Tuple{Tuple{Function, Vararg}, NamedTuple}}()
grid_fcall_with_background(plot_tangent_basis_glyphs!, args; z = z_sphere())
# Case 8: Nothing
args = nothing
grid_fcall_with_background(plot_tangent_basis_glyphs!, args; z = z_paraboloid(a = TestMatrices.r, b = 0.25* TestMatrices.r ))
# Case 9: Keywords
args = (; halfsize = 40)
grid_fcall_with_background(plot_tangent_basis_glyphs!, args; z = z_cos())
# Case 10: Vector of Nothing
args = [nothing, nothing, nothing]
grid_fcall_with_background(plot_tangent_basis_glyphs!, args; z = z_ellipsoid(; tilt = π / 6))
# Case 11: Vector of keywords
args = [(; halfsize = 35), (; halfsize = 40), (; halfsize = 45)]
grid_fcall_with_background(plot_tangent_basis_glyphs!, args; z = z_ridge_peak_valleys())
args = nothing
grid_fcall_with_background(plot_tangent_basis_glyphs!, args)



# Projected normal vector
args = (; multip = 100, maxg = 100)
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = z_cylinder_offset(π / 6), Δ = 50)

args = (; multip = 500, maxg = 100, rgb = 0.5 * PALETTE_GRGB[3])
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = z_ellipsoid(; tilt = π / 4), Δ = 50)

args = (; multip = 100, maxg = 100)
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = z_paraboloid())

args = (; multip = 100, maxg = 100)
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = z_paraboloid(; a = 500, b = 500))

args = (; multip = 100, maxg = 100)
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = z_sphere(), Δ = 50)

args = (; multip = 50)
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = z_cos(), Δ = 75)

args = (; multip = 50)
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = z_ridge_peak_valleys(), Δ = 25)

args = (; multip = 50, rgb = 0.5 * PALETTE_GRGB[3])
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = -z_paraboloid(; a = 500, b = 40))

args = (; multip = 50, rgb = 0.5 * PALETTE_GRGB[3])
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = z_exp3())



# Curvature glyphs
args = (; multip = 30000, rgb = 0.7 * PALETTE_GRGB[3])
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_ellipsoid(; tilt = π / 4), Δ = 50)

# Curvature glyphs (separate colors)
args = [(; directions = 1, multip = 12000, ming = -100, rgb = 0.7 * PALETTE_GRGB[3]),
        (; directions = 2, multip = 12000, ming = -100, rgb = PALETTE_GRGB[4])]
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_cylinder_offset(π / 6), Δ = 50)

args = [(; directions = 1, multip = 50000, ming = -50, rgb = 0.7 * PALETTE_GRGB[3]),
        (; directions = 2, multip = 50000, ming = -50, rgb = PALETTE_GRGB[4])]
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_ellipsoid(; tilt = π / 4), Δ = 50)

args = [(; directions = 1,  multip = 15000, maxg = 250, rgb = 0.7 * PALETTE_GRGB[3]),
        (; directions = 2,  multip = 15000, maxg = 250, rgb = PALETTE_GRGB[4])]
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_paraboloid())

args = [(; directions = 1,  multip = 10000, maxg = 250, rgb = 0.7 * PALETTE_GRGB[3]),
        (; directions = 2,  multip = 10000, maxg = 250, rgb = PALETTE_GRGB[4])]
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_paraboloid(; a = 400, b = 600), Δ = 50)

args = [(; multip = 12000, directions = 1, rgb = 0.7 * PALETTE_GRGB[3]),
        (; multip = 12000, directions = 2, rgb = PALETTE_GRGB[4])]
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_paraboloid())

args = [(; multip = 12000, directions = 1, rgb = 0.7 * PALETTE_GRGB[3]),
        (; multip = 12000, directions = 2, rgb = PALETTE_GRGB[4])]
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_sphere(), Δ = 50)

args = [(; multip = 5000, directions = 1, rgb = 0.7 * PALETTE_GRGB[3], maxg = 300),
        (; multip = 5000, directions = 2, rgb = PALETTE_GRGB[4], ming = -300)]
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_exp3(), Δ = 40)

args = [(; multip = 25000, directions = 1, rgb = 0.7 * PALETTE_GRGB[3], maxg = 25, ming = -25),
        (; multip = 25000, directions = 2, rgb = PALETTE_GRGB[4], maxg = 25, ming = -25)]
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_cos(;mult = 50), Δ = 15)

args = (; multip = 4000)
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_ridge_peak_valleys(), Δ = 25)

args = [(; multip = 7500, directions = 1, rgb = 0.99 * PALETTE_GRGB[3], maxg = 25, ming = -25),
        (; multip = 7500, directions = 2, rgb = PALETTE_GRGB[4], maxg = 25, ming = -25)]
grid_fcall_with_background(plot_curvature_glyphs!, args, z = -z_paraboloid(; a = 500, b = -300), Δ = 50)


# Curvature glyphs, secondary (smallest value) direction only
args = (; multip = 12000, ming = -30, directions = 2, rgb = 0.5 * PALETTE_GRGB[3])
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_cylinder(π / 6), Δ = 50)

args = (; multip = 10000, directions = 2, rgb = 0.2 * PALETTE_GRGB[3], ming = -30)
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_ellipsoid(; tilt = π / 4), Δ = 30)

args = (; multip = 20000, directions = 2, rgb = 0.5 * PALETTE_GRGB[3])
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_paraboloid())

args = (; multip = 10000, directions = 2, rgb = 0.5 * PALETTE_GRGB[3])
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_paraboloid(; a = 500, b = 500))

args = (; multip = 1000, directions = 2, rgb = 0.8 * PALETTE_GRGB[3], ming = -50)
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_cos(), Δ = 15)

args = (; multip = 4000, directions = 2, rgb = 0.5 * PALETTE_GRGB[3])
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_ridge_peak_valleys(), Δ = 25)




function sidelength(K::TENSORMAP)
    2 * maximum(abs.(K))
end
# Fallback
sidelength(glyphspec) = throw("unsupported for now: $(typeof(glyphspec))")

function side_for_glyph(z)
    obj = BidirectionOnGrid(𝐊!, z)
    d = Domain(TestMatrices.R, CartesianIndices((-2:2, -2:2)))
    glyph_side = map(CartesianIndices(z)) do I
        if d(I.I...)
            sidelength(obj(I))
        else
            0.0
        end
    end
    # Normalize
    #m = maximum(glyph_side)
    # TEMP debug
    #glyph_side[100:400, 100:500] .= 0.0
    #glyph_side[100:400, 500:900] .= 1.0
    #glyph_side ./ m
end


"""
    adjust_distribution(p; dst_minval = 0.01, dst_maxval = 1.0, src_minval = 0.2, src_maxval = 1.0, t = 0.5, slope = 1.2, r = 40
"""
function adjust_distribution!(p; dst_minval = 0.01, dst_maxval = 1.0, src_minval = 0.2, src_maxval = 1.0, t = 0.5, slope = 1.2, r = 40)
    # Remember that the sum of two opposite interpolation is a const plus a linear interpolation.
    @assert 0.0 <= maximum(p)
    @assert minimum(p) >= 0.0
    @assert 1.0 >= dst_maxval >= 0.0
    adjust_histogram!(p, ContrastStretching(;t, slope))
    adjust_histogram!(p, LinearStretching(; dst_minval, dst_maxval, src_minval, src_maxval))
    #=
    adjust_histogram!(p, GammaCorrection(; gamma))
    src_minval, src_maxval = extrema(p)
    adjust_histogram!(p, LinearStretching(;src_minval, src_maxval, dst_minval, dst_maxval))
    =#
    # If a 'large glyph' pixel is adjacent to a 'small glyph' pixel,
    # we probably shouldn't place a glyph on the 'small glyph' pixel.
    # Hence, we extend the low probability areas.
    #erode!(p; r)
end

"""
    glyph_probability(z; dst_minval = 0.01, dst_maxval = 1.0, src_minval = 0.2, src_maxval = 1.0, t = 0.5, slope = 1.2, r = 40)

Note that r is quite expensive. It decreases propbability 
in radius r from a low probability pixel.


TODO: Pick a better algorithm:
1) Select (too many) pixels
2) For each pixel, calculate probability.
3) Toss the dice and plot

Possibly, check if the vincinity is already covered by an earlier glyph.
"""
function glyph_probability(z; dst_minval = 0.01, dst_maxval = 1.0, src_minval = 0.2, src_maxval = 1.0, t = 0.5, slope = 1.2, r = 40)
    # Area needed for a glyph here. Values on the edges won't matter.
    p = side_for_glyph(z).^2
    # Likelihood (not normalized) that we will place a glyph on each pixel
    m = maximum(p)
    p .= (m .- p)
    # But there's not zero space left, even where the glyphs are largest
    adjust_distribution!(p; dst_minval, dst_maxval, src_minval, src_maxval, t, slope, r)
end

function pts_for_glyphs(z, n; dst_minval = 0.01, dst_maxval = 1.0, src_minval = 0.2, src_maxval = 1.0, t = 0.5, slope = 1.2, r = 40)
    spac = glyph_probability(z; dst_minval, dst_maxval, src_minval, src_maxval, t, slope, r)
    sample(MersenneTwister(123), CartesianIndices(spac), Weights(vec(spac)), n)
end


z = z_exp3()
args = [(; multip = 5000, directions = 2, rgb = PALETTE_GRGB[4], ming = -300),
        (; multip = 5000, directions = 1, rgb = PALETTE_GRGB[3], maxg = 300)]
begin
    img = background(z)
    dst_minval = 0.07
    pts = pts_for_glyphs(z, 8000; dst_minval)
    fcall(plot_curvature_glyphs!, args, img, z, pts)
end

z = z_cos(;mult = 50)
args = [(; multip = 10000, directions = 1, rgb = PALETTE_GRGB[3], maxg = 300),
        (; multip = 10000, directions = 2, rgb = PALETTE_GRGB[4], ming = -300)]
begin
    img = background(z)
    dst_minval = 0.03
    r = 40
    pts = pts_for_glyphs(z, 8000; dst_minval, r)
    fcall(plot_curvature_glyphs!, args, img, z, pts)
end

# The variation here is too large for this type of plot
z = z_ridge_peak_valleys()
args = [(; multip = 1000, directions = 1, rgb = PALETTE_GRGB[3], maxg = 30, ming = -30),
        (; multip = 1000, directions = 2, rgb = PALETTE_GRGB[4], maxg = 30, ming = -30)]
begin
    img = background(z)
    dst_minval = 0.03
    r = 40
    pts = pts_for_glyphs(z, 8000; dst_minval, r)
    fcall(plot_curvature_glyphs!, args, img, z, pts)
end

# TODO somehow limit the size by sending arguments to pts_for_glyphs
z = z_ellipsoid(; tilt = π / 4)
args = [(; multip = 1000, directions = 1, rgb = PALETTE_GRGB[3], maxg = 30, ming = -30),
        (; multip = 1000, directions = 2, rgb = PALETTE_GRGB[4], maxg = 30, ming = -30)]
begin
    img = background(z)
    dst_minval = 0.03
    r = 40
    pts = pts_for_glyphs(z, 8000; dst_minval, r)
    fcall(plot_curvature_glyphs!, args, img, z, pts)
end
