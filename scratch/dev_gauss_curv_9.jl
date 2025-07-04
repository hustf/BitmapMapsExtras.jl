# Improving on dev_gauss_curv_7.jl
#    Make colorful curvature plotting
#    Enable different max and min glyph limits
using Test
import ImageFiltering
import ImageFiltering.Kernel
import BitmapMaps
using BitmapMaps: scaleminmax00, line!, display_if_vscode, open_as_temp_in_imgedit
import ImageTransformations
using ImageTransformations: recenter, center, RotMatrix, WarpedView
import Interpolations
using Interpolations: OnGrid, Flat
import ImageCore
using ImageCore: N0f8, GrayA, RGBA, colormap, Gray, RGB
import ColorBlendModes
using ColorBlendModes: BlendLighten
import ColorTypes
using ColorTypes: alpha
import LinearAlgebra
using LinearAlgebra: ⋅, ×, norm, normalize!
using StaticArrays
using BenchmarkTools
import LinearSolve
using LinearSolve: init, solve!, LinearProblem, NormalCholeskyFactorization
using LinearSolve: OperatorCondition, OperatorAssumptions

######################
# Utilties
######################
includet("utilties_test_matrix.jl")
using .TestMatrices
using .TestMatrices: r, I0
includet("utilties_graphic_6.jl")
includet("utilties_tangent_5.jl")
includet("utilties_curvature_7.jl")
includet("plot_basis_curvature_8.jl")
"""
These sample angles effectively estimate curvature in the tangent plane.
Values are a compromise between tests on a sphere and on a cylinder with r=499.
"""
const VΦ = MVector{4, Float64}([4, 98, 133, 179] .* (π / 180 ))
# Possible test candidate: [11.25, 56.25, 101.25, 146.25]
# where Δ = 90 / 4:
#  Δ[0.5, 2.5, 4.5, 6.5] - we mirrored two of the four angles around 90°
# in order to avoid bias.

#########
# Testing 
#########
function test_plot_basis(; z = z_sphere(),
         pts = [I0 + CartesianIndex(-y , x ) for x in (-r - 1):100:r, y in (-r - 1):100:r],
         mi = -400.0, ma = 900.0)
    img = bluesc(z; mi, ma) .* 0.6
    plot_tangent_basis_glyphs!(img, z, pts)
end



#####################
# Tangent basis plots
#####################
test_plot_basis()
test_plot_basis(; z = -z_sphere())
extrema(z_ellipsoid())
test_plot_basis(; z = z_ellipsoid(), mi = 0.0, ma = 500.0)
extrema(z_paraboloid(; a = 1000))
test_plot_basis(; z = z_paraboloid(; a = 1000), mi = 0.0, ma = 373.)
extrema(z_paraboloid(; a= 1000, b = -500))
test_plot_basis(; z = z_paraboloid(; a= 1000, b = -500), mi = -250.0, ma = 125.0)

#################################
# Single curvature value
#    lenient acceptance criteria)
#################################
# 
@testset "Curvature at centre cylinder sphere" begin
    z = z_cylinder(0)
    K = principal_curvature_components(z, I0)
    @test isapprox(norm(K[:, 1]), 0, atol = 1e-4)
    @test isapprox(norm(K[:, 2]), 1 / r, atol = 1e-4)
    @test !is_quantity_positive(K[:,2])
    z = z_sphere()
    K = principal_curvature_components(z, I0)
    @test isapprox(norm(K[:, 1]), 1 / r, atol = 1e-3)
    @test isapprox(norm(K[:, 2]), 1 / r, atol = 1e-3)
    @test !is_quantity_positive(K[:,1])
    @test !is_quantity_positive(K[:,2])
end
@testset "Curvature at centre sadde and elliptic paraboloid" begin
    a = 1.0
    b = 0.5
    z = z_paraboloid(; a, b = -b)
    K = principal_curvature_components(z, I0)
    @test isapprox(norm(K[:, 1]), 1 / a, atol = 0.05)
    @test isapprox(norm(K[:, 2]), 1 / b, atol = 0.06)
    @test is_quantity_positive(K[:,1])
    @test !is_quantity_positive(K[:,2])
    z = z_paraboloid(; a, b)
    K = principal_curvature_components(z, I0)
    # Remember that the order is determined by sorting quantity values.
    @test isapprox(norm(K[:, 2]), 1 / a, atol = 0.03)
    @test isapprox(norm(K[:, 1]), 1 / b, atol = 0.05)
    @test is_quantity_positive(K[:,1])
    @test is_quantity_positive(K[:,2])
end


#############################
# Plot single curvature glyph
#############################

# Negative curvature
let pt = I0
    z = -z_paraboloid()
    multglyph = 100
    minglyph = -300; maxglyph = 300
    dashsize = maxglyph ÷ 10
    directions = 1:2
    K = principal_curvature_components(z, pt)
    @show norm(K[:, 1]) norm(K[:, 2])
    @test !is_quantity_positive(K[:,1])
    @test !is_quantity_positive(K[:,2])
    # Prepare plot
    f_is_within_limits = func_is_glyph_within_limits(directions, maxglyph, minglyph)
    bbuf = zeros(GrayA{N0f8}, size(z)...)
    # Scale and plot the single glyph
    plot_principal_directions_glyph!(bbuf, pt, directions, f_is_within_limits, dashsize, multglyph * K)
end

# Positive - zero curvature
let ∠ = π / 6, pt = I0
    z = -z_cylinder(∠)
    multglyph = 100000
    minglyph = -300; maxglyph = 300
    dashsize = maxglyph ÷ 10
    directions = 1:2
    K = principal_curvature_components(z, pt)
    @show norm(K[:, 1]) norm(K[:, 2])
    @test is_quantity_positive(K[:,1])
    # Prepare plot
    f_is_within_limits = func_is_glyph_within_limits(directions, maxglyph, minglyph)
    bbuf = zeros(GrayA{N0f8}, size(z)...)
    # Scale and plot the single glyph
    plot_principal_directions_glyph!(bbuf, pt, directions, f_is_within_limits, dashsize, multglyph * K)
end


#######################
# Plot curvature glyphs
#######################
   
function test_curvature_plot(; rngxy = (-r + 5):100:0,
        multglyph = 10000, minglyph = -100, maxglyph = 100,     z = z_cylinder(π/4), directions = 1:2)
    # Points
    vxy = sort(vcat(rngxy, 0, -rngxy))
    pts = [TestMatrices.I0 + CartesianIndex(x, y) for x in vxy, y in vxy]
    #    
    mi, ma = extrema(z)
    img = bluesc(z; mi, ma) .* 0.6
    # Add simple contour lines, too
    Δc = 20
    wc = Δc / 10
    map!(img, z, img) do zz, pix
        mod(zz, Δc) < wc ? RGBA{N0f8}(0.1, 0.1, 0.1, 1.0) : pix 
    end
    # Overlay the glyphs
    plot_curvature_glyphs!(img, z, pts; multglyph, minglyph, maxglyph, directions)
end

test_curvature_plot()
test_curvature_plot(; z = z_cylinder(π / 6))
test_curvature_plot(; z = -z_cylinder(π / 6))
test_curvature_plot(; z = -z_cylinder(0))
test_curvature_plot(; z = z_sphere())
test_curvature_plot(; z = z_sphere(), directions = 1:2) 
test_curvature_plot(; z = z_sphere(), directions = 2) # The smallest-curvature-value direction only
test_curvature_plot(; z = z_paraboloid(; a = 1000, b = -500), multglyph = 20000)
test_curvature_plot(; z = z_paraboloid(; a = 1000, b = 500), multglyph = 20000)


test_curvature_plot(; z = -z_ellipsoid(;rel_half_axis = 0.5), multglyph = 5000, directions = 2,
    rngxy = (-r + 5):10:0)
test_curvature_plot(; z = -z_cylinder(π / 6), rngxy = (-r + 5):20:0, multglyph = 5000)


# 66.158 ms (73012 allocations: 56.60 MiB)
@btime test_curvature_plot()




