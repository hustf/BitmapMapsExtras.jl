using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras.TestMatrices: I0
using BitmapMapsExtras: plot_bidirec_streamlines!, principal_curvature_components
using BitmapMapsExtras: 𝐊! #, DirectionAtXY, is_bidirec_vect_positive, norm
using BitmapMapsExtras: PALETTE_GRGB, get_bidirec_streamlines_points
import BitmapMaps
using BitmapMaps: mark_at!

!@isdefined(hashstr) && include("common.jl")


@test true

# WIP

r = TestMatrices.r
pts = [CartesianIndex(300, 450),
      CartesianIndex(300, 475),
      CartesianIndex(300, 500)]
args = ((𝐊!, true, true), (;rgb = PALETTE_GRGB[4], r = 5f0, dtmax = 1, tstop = 1900))
img = fcall_with_background(plot_bidirec_streamlines!, args, z_ellipsoid(; tilt = -π / 5 ), pts)
mark_at!(img, pts, 9, "in_circle")

# Maybe we ought to have a stop criterion at too large steepness
args = ((𝐊!, false, true), (;rgb = PALETTE_GRGB[4], r = 5f0, dtmax = 1, tstop = 1900))
img = fcall_with_background(plot_bidirec_streamlines!, args, z_ellipsoid(; tilt = -π / 5 ), pts)
mark_at!(img, pts, 9, "in_circle")

args = ((𝐊!, true, false), (;rgb = PALETTE_GRGB[4], r = 5f0, dtmax = 1, tstop = 1900))
img = fcall_with_background(plot_bidirec_streamlines!, args, z_ellipsoid(; tilt = -π / 3 ), pts)
mark_at!(img, pts, 9, "in_circle")


pts = [CartesianIndex(200, 450),
      CartesianIndex(220, 475),
      CartesianIndex(240, 500)]
args = ((𝐊!, false, false), (;rgb = PALETTE_GRGB[4], r = 5f0, dtmax = 1, tstop = 1550))
img = fcall_with_background(plot_bidirec_streamlines!, args, z_cos(), pts)
mark_at!(img, pts, 9, "in_circle")

# Narrow down debug turn
pts = [CartesianIndex(220, 235)]
args = ((𝐊!, false, false), (;rgb = PALETTE_GRGB[4], r = 1f0, dtmax = 0.1, tstop = 100))
img = fcall_with_background(plot_bidirec_streamlines!, args, z_cos(), pts)
mark_at!(img, pts, 9, "in_circle")







background(z_cos())
extrema(z_paraboloid(;a= 5r, b=-5r))



# Second run.
# 1.200752 seconds (3.45 M allocations: 1.399 GiB, 58.50% gc time)
# 0.940869 seconds (2.50 M allocations: 1.322 GiB, 50.24% gc time)
# 1.092556 seconds (2.23 M allocations: 1.310 GiB, 58.66% gc time)
# 0.790105 seconds (1.89 M allocations: 1.294 GiB, 41.84% gc time)
# 0.748579 seconds (1.67 M allocations: 1.286 GiB, 42.00% gc time)
@time let 
    # This format is pretty complicated, but also powerful. 
    vargs = [((𝐊!, false, true),  (rgb = PALETTE_GRGB[3], r = 3f0, strength = 0.4f0)),
            ((𝐊!, false, false), (rgb = PALETTE_GRGB[4], r = 2f0, strength = 0.4f0))]
    grid_fcall_with_background(plot_bidirec_streamlines!, vargs; Δ = 100, offset = (-30, -30))
end

# dtmax fixes errors where a trace loops back
# 0.888065 seconds (9.56 M allocations: 1.749 GiB, 60.76% gc time, 1 lock conflict)
# 3.064361 seconds (18.86 M allocations: 1.600 GiB, 20.63% gc time)
# 2.845686 seconds (13.12 M allocations: 1.343 GiB, 22.42% gc time)
# 2.225794 seconds (6.15 M allocations: 1.031 GiB, 20.39% gc time)
# 2.143907 seconds (1.64 M allocations: 855.565 MiB, 23.44% gc time)
@time let 
    vargs = [((𝐊!, false, true),  (rgb = PALETTE_GRGB[3], r = 1.5f0, strength = 0.4f0, dtmax = 1)),
            ((𝐊!, false, false), (rgb = PALETTE_GRGB[4], r = 1.5f0, strength = 0.4f0, dtmax = 1)), 
            ((𝐊!, true, true),  (rgb = PALETTE_GRGB[1], r = 1.5f0, strength = 0.4f0, dtmax = 1)),
            ((𝐊!, true, false), (rgb = PALETTE_GRGB[2], r = 1.0f0, strength = 0.4f0, dtmax = 1))]
    grid_fcall_with_background(plot_bidirec_streamlines!, vargs; Δ = 200)
end


let 
    vargs = [((𝐊!, false, true),  (rgb = PALETTE_GRGB[3], r = 3f0, strength = 0.4f0)),
            ((𝐊!, false, false), (rgb = PALETTE_GRGB[4], r = 2f0, strength = 0.4f0))]
    grid_fcall_with_background(plot_bidirec_streamlines!, vargs; Δ = 100, offset = (-30, -30), z = z_cylinder(π / 6))
end

# dtmax fixes undulating traces
let 
    vargs = [((𝐊!, false, true),  (rgb = PALETTE_GRGB[3], r = 3f0, strength = 0.4f0, dtmax = 1)),
            ((𝐊!, false, false), (rgb = PALETTE_GRGB[4], r = 2f0, strength = 0.4f0, dtmax = 1))]
    grid_fcall_with_background(plot_bidirec_streamlines!, vargs; Δ = 100, offset = (-30, -30), z = z_cylinder(π / 6))
end


let 
    vargs = [((𝐊!, false, true),  (rgb = PALETTE_GRGB[3], r = 1f0, strength = 0.4f0, dtmax = 1)),
            ((𝐊!, false, false), (rgb = PALETTE_GRGB[4], r = 1f0, strength = 0.4f0, dtmax = 1)), 
            ((𝐊!, true, true),  (rgb = PALETTE_GRGB[1], r = 1f0, strength = 0.4f0, dtmax = 1)),
            ((𝐊!, true, false), (rgb = PALETTE_GRGB[2], r = 0.6f0, strength = 0.4f0, dtmax = 1))]
    grid_fcall_with_background(plot_bidirec_streamlines!, vargs, z = z_ellipsoid())
end

let 
    vargs = [((𝐊!, false, true),  (rgb = PALETTE_GRGB[3], r = 1f0, strength = 0.4f0, dtmax = 1)),
            ((𝐊!, false, false), (rgb = PALETTE_GRGB[4], r = 1f0, strength = 0.4f0, dtmax = 1)), 
            ((𝐊!, true, true),  (rgb = PALETTE_GRGB[1], r = 1f0, strength = 0.4f0, dtmax = 1)),
            ((𝐊!, true, false), (rgb = PALETTE_GRGB[2], r = 0.6f0, strength = 0.4f0, dtmax = 1))]
    grid_fcall_with_background(plot_bidirec_streamlines!, vargs; Δ = 100, offset = (-30, -30), z = z_paraboloid())
end


# This encounters switch primary <-> secondary
let 
    vargs = [((𝐊!, false, true),  (rgb = PALETTE_GRGB[3], r = 3f0, strength = 0.4f0, dtmax = 1, tstop = 3000))]
    grid_fcall_with_background(plot_bidirec_streamlines!, vargs; Δ = 1500, offset = (-20, -20), z = z_ridge_peak_valleys())
end



#=
# DEV Allocs
using BitmapMapsExtras: UnidirectionAtXY, vu0_from_pts, make_tspan, callbacks_streamlines, MVector, ODEProblem, rhs!, update_corners!, VΦ
using BenchmarkTools

f = 𝐊!
primary = true
flip = false
z = z_cylinder(0)
uxy = UnidirectionAtXY(f, z, primary, flip)


@inferred uxy(30.4, 50.2)
# 3.375 μs (12 allocations: 560 bytes)
# 3.400 μs (12 allocations: 496 bytes)
# 3.350 μs (10 allocations: 400 bytes)
# 3.388 μs (10 allocations: 400 bytes)
# 3.388 μs (9 allocations: 352 bytes)
# 3.350 μs (8 allocations: 320 bytes)
# 3.650 μs (19 allocations: 608 bytes) # Change from using update_corners
# 3.625 μs (15 allocations: 352 bytes) # Parametric type also for lpc
# 3.175 μs (2 allocations: 32 bytes)   # use parameters in BidirectionInDomain
# 3.075 μs (0 allocations: 0 bytes)    # use parameters in UniDirectionAtXy
@btime uxy(30.4, 50.2)

# 3.450 μs (11 allocations: 624 bytes)
# 3.375 μs (10 allocations: 608 bytes)
# 3.737 μs (21 allocations: 896 bytes) # Change from using update_corners
# 3.675 μs (17 allocations: 640 bytes)
# 3.225 μs (4 allocations: 240 bytes)   # use parameters in BidirectionInDomain
# 3.138 μs (1 allocation: 48 bytes)     # use parameters in BidirectionAtXy
@btime uxy.baxy(30.4, 50.2)
baxy = uxy.baxy
# 3.625 μs (19 allocations: 608 bytes) # Change from using update_corners
# 3.612 μs (15 allocations: 352 bytes)
# 3.112 μs (2 allocations: 32 bytes)
# 3.075 μs (0 allocations: 0 bytes) # use parameters in BidirectionAtXy
@btime baxy(30.4, 50.2)


bid = uxy.baxy.bid # BidirectionInDomain
# 3.525 μs (11 allocations: 288 bytes)
# 3.100 μs (0 allocations: 0 bytes) # use parameters in BidirectionInDomain
# 3.038 μs (0 allocations: 0 bytes) # re-introduce column swapping
@btime bid(30.4, 50.2)



function pro()
    for i = 0:10e4
        uxy(30.4, 50.2)
    end
end
@profview_allocs pro()


pt = (100, 110)
uxy.baxy.bid.bdog(pt...)

# 1.010 μs (6 allocations: 528 bytes)
# 990.000 ns (6 allocations: 608 bytes)
@btime uxy.baxy.bid.bdog(100, 110)
bdog = uxy.baxy.bid.bdog
# 759.322 ns (0 allocations: 0 bytes)
@btime bdog(100, 110)

=#
