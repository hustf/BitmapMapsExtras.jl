using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: plot_glyphs!, plot_glyphs, PALETTE_GRGB
using BitmapMapsExtras: BidirectionOnGrid, 𝐊!, Domain, TENSORMAP
using BitmapMapsExtras: GSTensor, GSTangentBasis, GSVector, GlyphSpec
using BitmapMapsExtras: is_in_limits, MAG_EPS, VECTOR_REL_HALFWIDTH
using BitmapMapsExtras: DirectionOnGrid, 𝐧ₚᵤ!, 𝐧ₚ!, VECTOR_REL_HALFWIDTH
using BitmapMapsExtras: potential_scattered_placements, intersection_distance_ellipse
using BitmapMapsExtras: intersection_distance_glyph, pack_glyphs!
using BitmapMapsExtras.BitmapMaps: mark_at!
import StatsBase
using StatsBase:  norm
import Random
using Random: MersenneTwister


!@isdefined(hashstr) && include("common.jl")

@test_true
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


# Projected normal vector, direct call
plot_glyphs(z_ellipsoid(; tilt = π / 4), grid_indices((999, 999)); 
    gs = GSVector(;multip = 1000, maxg = 200, dashsize = 4))

# Projected normal vector
args = (; gs = GSVector(;multip = 100, maxg = 100))
grid_fcall_with_background(plot_glyphs!, args, z = z_cylinder_offset(π / 6), Δ = 50)

args = (; gs = GSVector(multip = 500, maxg = 100, color = 0.5 * PALETTE_GRGB[3]))
grid_fcall_with_background(plot_glyphs!, args, z = z_ellipsoid(; tilt = π / 4), Δ = 50)

args = (; gs = GSVector(multip = 100, maxg = 100))
grid_fcall_with_background(plot_glyphs!, args, z = z_paraboloid())

args = (; gs = GSVector(multip = 100, maxg = 100))
grid_fcall_with_background(plot_glyphs!, args, z = z_paraboloid(; a = 500, b = 500))

args = (; gs = GSVector(multip = 100, maxg = 100))
grid_fcall_with_background(plot_glyphs!, args, z = z_sphere(), Δ = 50)

args = (; gs = GSVector(multip = 50))
grid_fcall_with_background(plot_glyphs!, args, z = z_cos(), Δ = 75)

args = (; gs = GSVector(multip = 50))
grid_fcall_with_background(plot_glyphs!, args, z = z_ridge_peak_valleys(), Δ = 25)

args = (; gs = GSVector(multip = 50, color = 0.5 * PALETTE_GRGB[3]))
grid_fcall_with_background(plot_glyphs!, args, z = -z_paraboloid(; a = 500, b = 40))

args = (; gs = GSVector(multip = 50, color = 0.5 * PALETTE_GRGB[3]))
grid_fcall_with_background(plot_glyphs!, args, z = z_exp3())



# Curvature glyphs, direct call

plot_glyphs(z_ellipsoid(; tilt = π / 4), grid_indices((999, 999)), GSTensor(;multip = 30000, directions = 1))
plot_glyphs(z_ellipsoid(; tilt = π / 4), grid_indices((999, 999)), GSTensor(;multip = 30000, directions = 2))
plot_glyphs(z_ellipsoid(; tilt = π / 4), grid_indices((999, 999)), GSTensor(;multip = 30000, directions = 1:2))

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

# Smartly spaced glyphs

@test potential_scattered_placements(CartesianIndices((9, 9)), avgdist = 1.41) == CartesianIndex.([(1, 7),
       (1, 1),
       (3, 7),
       (6, 1),
       (7, 4),
       (7, 7),
       (5, 5),
       (8, 4),
       (5, 6),
       (2, 7)])


x1, y1 = -2, -2
x2, y2 = -1, -1

@test intersection_distance_ellipse(x1, y1, x2, y2, π / 4) ≈ 2√2
@test intersection_distance_ellipse(x1, y1, x2, y2, 3π / 4) ≈ √2
@test intersection_distance_ellipse(x1, y1, x2, y2, 5π / 4) ≈ 2√2
@test intersection_distance_ellipse(x1, y1, x2, y2, 7π / 4) ≈ √2
@test √2 < intersection_distance_ellipse(x1, y1, x2, y2, π / 2) < 2√2
@test √2 < intersection_distance_ellipse(x1, y1, x2, y2, π / 2) < 2√2

x1, y1 = -2, -2
x2, y2 = 1, 1

@test intersection_distance_ellipse(x1, y1, x2, y2, π / 4) ≈ 2√2
@test intersection_distance_ellipse(x1, y1, x2, y2, 3π / 4) ≈ √2
@test intersection_distance_ellipse(x1, y1, x2, y2, 5π / 4) ≈ 2√2
@test intersection_distance_ellipse(x1, y1, x2, y2, 7π / 4) ≈ √2
@test √2 < intersection_distance_ellipse(x1, y1, x2, y2, π / 2) < 2√2
@test √2 < intersection_distance_ellipse(x1, y1, x2, y2, π / 2) < 2√2

x1, y1 = 0, 2
x2, y2 = 1, 0
@test intersection_distance_ellipse(x1, y1, x2, y2, π / 2) ≈ 2
@test intersection_distance_ellipse(x1, y1, x2, y2, 3π / 2) ≈ 2
@test intersection_distance_ellipse(x1, y1, x2, y2, 0) ≈ 1
@test intersection_distance_ellipse(x1, y1, x2, y2, π) ≈ 1
@test 1 < intersection_distance_ellipse(x1, y1, x2, y2, π / 4) < 2
@test 1 < intersection_distance_ellipse(x1, y1, x2, y2, 5π / 4) < 2



# unit tests
gs = GSVector()
# Horizontal vector
v = [0.1, 0.0]
# To the other point 
α = 0.0
intersection_distance_glyph(gs, v, α)






#= This function:
    length(filtered_placements) = 2461
  0.425071 seconds (539.44 k allocations: 62.942 MiB, 5.91% gc time)
    length(filtered_placements) = 655
  0.206523 seconds (230.81 k allocations: 45.922 MiB, 5.69% gc time)
=#

# Packed gradient
pack_glyphs_with_background(z_paraboloid(), GSVector())


# Packed curvature
# Used to pass 5 glyphs
# Then, compare if Any[] slows things down....
img = pack_glyphs_with_background(z_sphere()[100:300, 100:300], GSTensor(multip = 12000)) 
# 0.194754 seconds (52.19 k allocations: 37.286 MiB, 6.93% gc time)
#  0.201228 seconds (52.19 k allocations: 37.287 MiB, 9.67% gc time)
#  0.323076 seconds (52.19 k allocations: 37.287 MiB, 42.70% gc time)
#  0.186682 seconds (52.19 k allocations: 37.287 MiB, 3.51% gc time)
# 0.188074 seconds (60.06 k allocations: 37.507 MiB, 3.32% gc time) TENSORMAP[] -> Any[]
# 0.195374 seconds (52.19 k allocations: 37.287 MiB, 7.49% gc time) typed passed_values
@time  img = pack_glyphs(z_cylinder(0), GSTensor(multip = 12000)) 

#for scatterdist = 201:-10:1
#    @show scatterdist
#    @time img = pack_glyphs(z_cylinder(0), GSTensor(multip = 12000); scatterdist)
#    display(img) 
#end




# More examples
pack_glyphs_with_background(z_cylinder_offset(π / 6), GSTensor(multip = 12000))
gs = GSTensor(multip = 5000, ming = -25)
pack_glyphs_with_background(z_ellipsoid(; tilt = π / 4), gs; scatterdist = 1)
gs = GSTensor(multip = 5000, strength = 10, ming = -25)
pack_glyphs_with_background(z_ellipsoid(; tilt = π / 4, a= 0.3), gs; scatterdist = 1)
pack_glyphs_with_background(z_paraboloid(), GSTensor(multip = 15000))
pack_glyphs_with_background(z_paraboloid(; a = 400, b = 600), GSTensor(multip = 10000 ))
pack_glyphs_with_background(z_sphere(), GSTensor(multip = 10000))
img = pack_glyphs_with_background(z_exp3(), GSTensor(multip = 7000))
pack_glyphs_with_background!(img, z_exp3(), GSTensor(multip = 7000);  seed = MersenneTwister(1))
pack_glyphs_with_background!(img, z_exp3(), GSTensor(multip = 7000);  seed = MersenneTwister(2))
pack_glyphs_with_background(z_cos(; mult = 50), GSTensor(multip = 25000))
img = pack_glyphs_with_background(z_ridge_peak_valleys(), GSTensor(multip = 4000 ), scatterdist = 1)
pack_glyphs!(img, z_ridge_peak_valleys(), GSTensor(multip = 4000), seed = MersenneTwister(2))





# Packed gradient
pack_glyphs_with_background(z_paraboloid()[1:100, 1:100], GSVector(multip = 15000))





# DEV /
#=
using BenchmarkTools 
# 36.886 ms (63 allocations: 12.50 MiB)
# 36.811 ms (74 allocations: 12.50 MiB)
@btime plot_tangent_basis_glyphs!(img, z, [CartesianIndex((315, 215))])

function pro()
    for i = 1:1000
        plot_tangent_basis_glyphs!(img, z, [CartesianIndex((315, 215))])
    end
end
@profview plot_tangent_basis_glyphs!(img, z, [CartesianIndex((315, 215))])

@profview_allocs pro()
=#
# / DEV


# DEV function signatures


b = BidirectionOnGrid(𝐊!, z_cos(; mult = 50)[400:600, 400:600])
img = background(z_cos(; mult = 50)[400:600, 400:600])
#length(passed_placements) = 11
#  0.010389 seconds (90.91 k allocations: 5.896 MiB)
@time pack_glyphs!(img, b, GSTensor(multip = 15000))

pack_glyphs_with_background(z_paraboloid()[1:100, 1:100], GSVector(multip = 15000))
