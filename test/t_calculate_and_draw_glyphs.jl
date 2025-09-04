using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: plot_tangent_basis_glyphs, plot_tangent_basis_glyphs!
using BitmapMapsExtras: plot_curvature_glyphs, plot_curvature_glyphs!
using BitmapMapsExtras: plot_𝐧ₚ_glyphs!, PALETTE_GRGB
using BitmapMapsExtras: BidirectionOnGrid, 𝐊!, Domain, TENSORMAP
using BitmapMapsExtras: GSTensor, GSTangentBasis, GSVector
import StatsBase
using StatsBase: Weights, sample, norm
import Random
using Random: MersenneTwister
import ImageContrastAdjustment
using ImageContrastAdjustment: LinearStretching, ContrastStretching, adjust_histogram!

!@isdefined(hashstr) && include("common.jl")

# Tangent basis
z = z_paraboloid()
img = background(z)
plot_tangent_basis_glyphs!(img, z, [CartesianIndex((315, 215))])


# DEV /
#=
using BenchmarkTools
# 36.886 ms (63 allocations: 12.50 MiB)
@btime plot_tangent_basis_glyphs!(img, z, [CartesianIndex((315, 215))])

@profview plot_tangent_basis_glyphs!(img, z, [CartesianIndex((315, 215))])

@profview_allocs plot_tangent_basis_glyphs!(img, z, [CartesianIndex((315, 215))])
=#
# / DEV


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
args = (; gs = GSTangentBasis(;halfsize = 40))
grid_fcall_with_background(plot_tangent_basis_glyphs!, args; z = z_cos())
# Case 10: Vector of Nothing
args = [nothing, nothing, nothing]
grid_fcall_with_background(plot_tangent_basis_glyphs!, args; z = z_ellipsoid(; tilt = π / 6))
# Case 11: Vector of keywords

args = [(; gs = GSTangentBasis(;halfsize = 35)), (; gs = GSTangentBasis(;halfsize = 40)), (; gs = GSTangentBasis(;halfsize = 40))]
grid_fcall_with_background(plot_tangent_basis_glyphs!, args; z = z_ridge_peak_valleys())

args = nothing
grid_fcall_with_background(plot_tangent_basis_glyphs!, args)



# Projected normal vector
args = (; gs = GSVector(;multip = 100, maxg = 100))
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = z_cylinder_offset(π / 6), Δ = 50)

args = (; gs = GSVector(multip = 500, maxg = 100, color = 0.5 * PALETTE_GRGB[3]))
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = z_ellipsoid(; tilt = π / 4), Δ = 50)

args = (; gs = GSVector(multip = 100, maxg = 100))
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = z_paraboloid())

args = (; gs = GSVector(multip = 100, maxg = 100))
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = z_paraboloid(; a = 500, b = 500))

args = (; gs = GSVector(multip = 100, maxg = 100))
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = z_sphere(), Δ = 50)

args = (; gs = GSVector(multip = 50))
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = z_cos(), Δ = 75)

args = (; gs = GSVector(multip = 50))
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = z_ridge_peak_valleys(), Δ = 25)

args = (; gs = GSVector(multip = 50, color = 0.5 * PALETTE_GRGB[3]))
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = -z_paraboloid(; a = 500, b = 40))

args = (; gs = GSVector(multip = 50, color = 0.5 * PALETTE_GRGB[3]))
grid_fcall_with_background(plot_𝐧ₚ_glyphs!, args, z = z_exp3())

# Curvature glyphs, direct call

plot_curvature_glyphs(z_ellipsoid(; tilt = π / 4), grid_indices((999, 999)); gs = GSTensor(;multip = 30000, directions = 1))
# 46.051 ms (570 allocations: 15.26 MiB)
plot_curvature_glyphs(z_ellipsoid(; tilt = π / 4), grid_indices((999, 999)); gs = GSTensor(;multip = 30000, directions = 2))
# 56.357 ms (705 allocations: 19.07 MiB)
plot_curvature_glyphs(z_ellipsoid(; tilt = π / 4), grid_indices((999, 999)); gs = GSTensor(;multip = 30000, directions = 1:2))

# Curvature glyphs, indirect call

args = (; gs = GSTensor(;multip = 30000, color1 = PALETTE_GRGB[1])) 
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_ellipsoid(; tilt = π / 4), Δ = 50)

# Curvature glyphs (separate colors)
args = (; gs = GSTensor(multip = 12000))
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_cylinder_offset(π / 6), Δ = 50)

args = (; gs = GSTensor(multip = 50000, strength = 10))
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_ellipsoid(; tilt = π / 4), Δ = 50)

args = [(; gs = GSTensor(directions = 1,  multip = 15000, maxg = 250, strength = 10)),
        (; gs = GSTensor(directions = 2,  multip = 15000, maxg = 250, strength = 10))]
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_paraboloid())

args = [(; gs = GSTensor(multip = 10000))]
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_paraboloid(; a = 400, b = 600), Δ = 50)

args = (; gs = GSTensor(multip = 12000))
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_paraboloid())

args = (; gs = GSTensor(multip = 12000))
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_sphere(), Δ = 50)

args = (; gs = GSTensor(multip = 5000, maxg = 300, dashsize = 3))
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_exp3(), Δ = 40)

args = (; gs = GSTensor(multip = 25000, maxg = 25, ming = -25))
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_cos(;mult = 50), Δ = 15)

args = (; gs = GSTensor(multip = 4000))
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_ridge_peak_valleys(), Δ = 25)

args = (; gs = GSTensor(multip = 7500))
grid_fcall_with_background(plot_curvature_glyphs!, args, z = -z_paraboloid(; a = 500, b = -300), Δ = 50)


# Curvature glyphs, secondary (smallest value) direction only
args = (; gs = GSTensor(multip = 12000, ming = -30, directions = 2, strength = 10))
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_cylinder(π / 6), Δ = 50)

args = (; gs = GSTensor(; multip = 10000, directions = 2, ming = -30, strength = 10))
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_ellipsoid(; tilt = π / 4), Δ = 30)

args = (; gs = GSTensor(; multip = 20000, directions = 2, strength = 10))
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_paraboloid())

args = (; gs = GSTensor(multip = 22000, directions = 2, strength = 10))
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_paraboloid(; a = 500, b = 500))

args = (; gs = GSTensor(multip = 1000, directions = 2, ming = -50, strength = 10))
grid_fcall_with_background(plot_curvature_glyphs!, args, z = z_cos(), Δ = 15)

args = (; gs = GSTensor(multip = 4000, directions = 2))
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
args = [(; gs = GSTensor(multip = 5000, directions = 2, ming = -300)),
        (; gs = GSTensor(multip = 5000, directions = 1, maxg = 300))]
begin
    img = background(z)
    dst_minval = 0.07
    pts = pts_for_glyphs(z, 8000; dst_minval)
    fcall(plot_curvature_glyphs!, args, img, z, pts)
end

z = z_cos(;mult = 50)
args = [(; gs = GSTensor(multip = 10000, directions = 1, maxg = 300)),
        (; gs = GSTensor(multip = 10000, directions = 2, ming = -300))]
begin
    img = background(z)
    dst_minval = 0.03
    r = 40
    pts = pts_for_glyphs(z, 8000; dst_minval, r)
    fcall(plot_curvature_glyphs!, args, img, z, pts)
end

# The variation here is too large for this type of plot
z = z_ridge_peak_valleys()
args = [(; gs = GSTensor(multip = 1000, directions = 1, maxg = 30, ming = -30)),
        (; gs = GSTensor(multip = 1000, directions = 2, maxg = 30, ming = -30))]
begin
    img = background(z)
    dst_minval = 0.03
    r = 40
    pts = pts_for_glyphs(z, 8000; dst_minval, r)
    fcall(plot_curvature_glyphs!, args, img, z, pts)
end

# TODO somehow limit the size by sending arguments to pts_for_glyphs
z = z_ellipsoid(; tilt = π / 4)
args = [(; gs = GSTensor(multip = 1000, directions = 1, maxg = 30, ming = -30)),
        (; gs = GSTensor(multip = 1000, directions = 2, maxg = 30, ming = -30))]
begin
    img = background(z)
    dst_minval = 0.03
    r = 40
    pts = pts_for_glyphs(z, 8000; dst_minval, r)
    fcall(plot_curvature_glyphs!, args, img, z, pts)
end
