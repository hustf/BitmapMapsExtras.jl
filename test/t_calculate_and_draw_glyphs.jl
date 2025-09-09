using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: plot_tangent_basis_glyphs, plot_tangent_basis_glyphs!
using BitmapMapsExtras: plot_curvature_glyphs, plot_curvature_glyphs!
using BitmapMapsExtras: plot_𝐧ₚ_glyphs!, PALETTE_GRGB
using BitmapMapsExtras: BidirectionOnGrid, 𝐊!, Domain, TENSORMAP
using BitmapMapsExtras: GSTensor, GSTangentBasis, GSVector, GlyphSpec
using BitmapMapsExtras: is_in_limits, MAG_EPS, VECTOR_REL_HALFWIDTH
using BitmapMapsExtras: is_bidirec_vect_positive
using BitmapMapsExtras.BitmapMaps: mark_at!


import StatsBase
using StatsBase: Weights, sample, norm
import Random
using Random: MersenneTwister, randperm
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

# Smartly spaced glyphs
"""
    potential_scattered_placements(R::CartesianIndices; avgdist = 10, seed = MersenneTwister(123))
    --> Vector{CartesianIndex}

Return a vector of random CartesianIndex scattered within `R` such that the
mean nearest-neighbour distance is approximately `avgdist`. Refer Poisson disk sampling.

# Arguments

- `R` defines the possible values of each returned index. 
- `avgdist` (default keyword value 10.0) is the mean nearest-neighbour distance between points. The 
    resulting set of points will vary slightly from this target.
"""
function potential_scattered_placements(R::CartesianIndices; avgdist = 10.0, seed = MersenneTwister(123))
    area = length(R)
    # Estimate number of points, within bounds.
    n = max(1, min(round(Int, area / (4 * avgdist^2)), area))
    # Draw random unique indices
    chosen = randperm(seed, area)[1:n]
    [R[i] for i in chosen]
end

function potential_scattered_placements(fz, gs::GlyphSpec; scatterdist = gs.dashsize, seed = MersenneTwister(123))
    @assert hasfield(typeof(fz), :z)
    @assert hasfield(typeof(fz), :Ω)
    R = CartesianIndices(fz.z)
    Ω = fz.Ω
    minj = R[1][2] - Ω[1][2]
    maxi = R[end][1] - Ω[end][2]
    maxj = R[end][2] - Ω[end][2]
    mini = R[1][1] - Ω[1][2]
    Ri = CartesianIndices((mini:maxi, minj:maxj))
    potential_scattered_placements(Ri; avgdist = scatterdist, seed)
end



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

function required_radius_for_plotting(gs::GSTensor, K)
     @assert gs.directions == 1:2 "One direction not implemented. "
    # Extract the half-length of primary and secondary principal direction
    # glyphs
    if is_in_limits(gs, K)
        K1 = @view K[:, 1]
        K2 = @view K[:, 2]
        l1 = norm(K1) * gs.multip
        l2 = norm(K2) * gs.multip
    else
        l1 = 0.0
        l2 = 0.0
    end
    max(l1, l2) 
end


"""
    intersection_distance(x1, y1, x2, y2, α)

Approximate distance from origin to ellipse boundary along ray at angle `α`.

- `(x1, y1)` is assumed to lie on the major axis (used for orientation + length).
- `(x2, y2)` only contributes its distance as semi-minor axis length.
- `α` is ray angle in radians (CCW from x-axis).
"""
function intersection_distance(x1, y1, x2, y2, α)
    a = hypot(x1, y1)             # major radius
    b = hypot(x2, y2)             # minor radius (approx)
    @assert a >= b
    θ = atan(y1, x1)              # orientation
    ϕ = α - θ
    denom = (cos(ϕ)^2)/a^2 + (sin(ϕ)^2)/b^2
    return 1 / sqrt(denom)
end

x1, y1 = 2, 2
x2, y2 = 1, 1

@test intersection_distance(x1, y1, x2, y2, π / 4) ≈ 2√2
@test intersection_distance(x1, y1, x2, y2, 3π / 4) ≈ √2
@test intersection_distance(x1, y1, x2, y2, 5π / 4) ≈ 2√2
@test intersection_distance(x1, y1, x2, y2, 7π / 4) ≈ √2
@test √2 < intersection_distance(x1, y1, x2, y2, π / 2) < 2√2
@test √2 < intersection_distance(x1, y1, x2, y2, π / 2) < 2√2

x1, y1 = -2, -2
x2, y2 = -1, -1

@test intersection_distance(x1, y1, x2, y2, π / 4) ≈ 2√2
@test intersection_distance(x1, y1, x2, y2, 3π / 4) ≈ √2
@test intersection_distance(x1, y1, x2, y2, 5π / 4) ≈ 2√2
@test intersection_distance(x1, y1, x2, y2, 7π / 4) ≈ √2
@test √2 < intersection_distance(x1, y1, x2, y2, π / 2) < 2√2
@test √2 < intersection_distance(x1, y1, x2, y2, π / 2) < 2√2

x1, y1 = -2, -2
x2, y2 = 1, 1

@test intersection_distance(x1, y1, x2, y2, π / 4) ≈ 2√2
@test intersection_distance(x1, y1, x2, y2, 3π / 4) ≈ √2
@test intersection_distance(x1, y1, x2, y2, 5π / 4) ≈ 2√2
@test intersection_distance(x1, y1, x2, y2, 7π / 4) ≈ √2
@test √2 < intersection_distance(x1, y1, x2, y2, π / 2) < 2√2
@test √2 < intersection_distance(x1, y1, x2, y2, π / 2) < 2√2

x1, y1 = 0, 2
x2, y2 = 1, 0
@test intersection_distance(x1, y1, x2, y2, π / 2) ≈ 2
@test intersection_distance(x1, y1, x2, y2, 3π / 2) ≈ 2
@test intersection_distance(x1, y1, x2, y2, 0) ≈ 1
@test intersection_distance(x1, y1, x2, y2, π) ≈ 1
@test 1 < intersection_distance(x1, y1, x2, y2, π / 4) < 2
@test 1 < intersection_distance(x1, y1, x2, y2, 5π / 4) < 2


"""
    intersection_distance_glyph(gs::GSTensor, K::TENSORMAP, α)

Approximate distance from origin along ray at angle `α` to eclipsing elliptical boundary .
"""
function intersection_distance_glyph(gs, K, α)
     @assert gs.directions == 1:2 "One direction not implemented. "
    # Extract the half-length of primary and secondary principal direction
    # glyphs.
    # 
    # (x1, y1): Longest vector
    # (x2, y2): Shortest vector
    # The local 'y' is in a right handed, y up, x right coordinate system
    if norm(K[:, 1]) >= norm(K[:, 2])
        x1, x2 = K[1, :]
        y1, y2 = -K[2, :]
    else
        x2, x1 = K[1, :]
        y2, y1 = -K[2, :]
    end
    l1 = hypot(x1, y1)
    l2 = hypot(x2, y2)
    if l2 / l1 < 5 * VECTOR_REL_HALFWIDTH
        # The minor axis is so short compared to the other axis
        # that the determining 'width' is the half-cone axis of the vector glyph.
        # Since we don't want completely neighbouring parallel vectors,
        # we also include a factor of 3.
        # The direction is set to 90° off (x1, y1)
        #x2 = -y1 * VECTOR_REL_HALFWIDTH
        #y2 = x1 * VECTOR_REL_HALFWIDTH
        x2norm = x2 / l2
        y2norm = y2 / l2
        minlen = 5 * l1 * VECTOR_REL_HALFWIDTH
        x2 = x2norm * minlen
        y2 = y2norm * minlen
    end
    idist =  gs.multip * intersection_distance(x1, y1, x2, y2, α)
    # 
end

#= This function:
    length(filtered_placements) = 2461
  0.425071 seconds (539.44 k allocations: 62.942 MiB, 5.91% gc time)
    length(filtered_placements) = 655
  0.206523 seconds (230.81 k allocations: 45.922 MiB, 5.69% gc time)
=#
function placements_and_values(b, gs, ppts)
    passed_placements = CartesianIndex{2}[]
    passed_values = TENSORMAP[]
    passed_radii = Float64[]
    se = Set(ppts)
    while !isempty(se)
        # Set current point, drop it from the set.
        pt = pop!(se)
        K = b(pt)
        r = required_radius_for_plotting(gs, K)
        if r > MAG_EPS
            if isempty(passed_placements)
                # Only this point to consider: There must be room.
                push!(passed_placements, pt)
                push!(passed_values, copy(K))
                push!(passed_radii, r)
            else
                # Coarse filtering:
                # The locations for which glyph enveloping 
                # circles would overlap
                i_candidates = findall(1:length(passed_placements)) do i
                    npt = passed_placements[i]
                    nr = passed_radii[i]
                    r + nr >  norm(npt.I .- pt.I)
                end
                # Finer filtering:
                # We approximate an octahedron around pt and around each npt.
                any_too_close = any(i_candidates) do i
                    npt = passed_placements[i]
                    nr = passed_radii[i]
                    # From pt to npt
                    v = npt.I .- pt.I
                    # The angle is calculated as if we're in a 'y up' c.s.
                    αv = atan(-v[1], v[2]) 
                    # Intersection with pt's ellipse
                    rpt = intersection_distance_glyph(gs, K, αv)
                    # Npt's ellipse. Unfortunately, this slightly heavy calculation is
                    # not stored for the next loop.
                    Kn = b(npt)
                    rnpt = intersection_distance_glyph(gs, Kn, αv + π)
                    # The ellipses overlap when
                    rpt + rnpt >  norm(v)
                end
                # Would passing this glyph
                # overlap already passed glyphs?
                if ! any_too_close
                    push!(passed_placements, pt)
                    push!(passed_values, copy(K))
                    push!(passed_radii, r)
                end
            end
        end
    end
    passed_placements, passed_values
end

function pack_curvature_glyphs!(img, z, gs::GSTensor; scatterdist = gs.dashsize, seed = MersenneTwister(123))
    b = BidirectionOnGrid(𝐊!, z)
    ppts = potential_scattered_placements(b, gs; scatterdist, seed)
    # DEBUG
    # mark_at!(img, ppts, 1, "on_circle")
    filtered_placements, filtered_values = placements_and_values(b, gs, ppts)
    @show length(filtered_placements)
    # DEBUG
    #for (pt, K) in zip(filtered_placements, filtered_values)
    #    r = required_radius_for_plotting(gs, K)
    #    diameter = 1 + 2 * (Int(floor(r)))
    #    mark_at!(img, [pt], diameter, "on_circle")
    #end
    plot_curvature_glyphs!(img, gs, filtered_placements, filtered_values)
    img
end
function pack_curvature_glyphs(z, gs::GSTensor; scatterdist = gs.dashsize, seed = MersenneTwister(123))
    img = background(z)
    pack_curvature_glyphs!(img, z, gs; scatterdist, seed)
end

img = pack_curvature_glyphs(z_sphere()[100:300, 100:300], GSTensor(multip = 12000)) 
@time  img = pack_curvature_glyphs(z_cylinder(0), GSTensor(multip = 12000)) 

#for scatterdist = 201:-10:1
#    @show scatterdist
#    @time img = pack_curvature_glyphs(z_cylinder(0), GSTensor(multip = 12000); scatterdist)
#    display(img) 
#end




# More examples
pack_curvature_glyphs(z_cylinder_offset(π / 6), GSTensor(multip = 12000))

gs = GSTensor(multip = 5000, ming = -25)
pack_curvature_glyphs(z_ellipsoid(; tilt = π / 4), gs; scatterdist = 1)

gs = GSTensor(multip = 5000, strength = 10, ming = -25)
pack_curvature_glyphs(z_ellipsoid(; tilt = π / 4, a= 0.3), gs; scatterdist = 1)

# This would benefit from a more advanced collision detection.
pack_curvature_glyphs(z_paraboloid(), GSTensor(multip = 15000))

pack_curvature_glyphs(z_paraboloid(; a = 400, b = 600), GSTensor(multip = 10000 ))

pack_curvature_glyphs(z_sphere(), GSTensor(multip = 10000))

img = pack_curvature_glyphs(z_exp3(), GSTensor(multip = 7000))
pack_curvature_glyphs!(img, z_exp3(), GSTensor(multip = 7000);  seed = MersenneTwister(1))
pack_curvature_glyphs!(img, z_exp3(), GSTensor(multip = 7000);  seed = MersenneTwister(2))

pack_curvature_glyphs(z_cos(; mult = 50), GSTensor(multip = 25000))

img = pack_curvature_glyphs(z_ridge_peak_valleys(), GSTensor(multip = 4000 ), scatterdist = 1)
pack_curvature_glyphs!(img, z_ridge_peak_valleys(), GSTensor(multip = 4000), seed = MersenneTwister(2))

