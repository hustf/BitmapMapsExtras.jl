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



# Second run. TODO: Drop allocations in tangent basis (𝐊!)
# 1.200752 seconds (3.45 M allocations: 1.399 GiB, 58.50% gc time)
@time let 
    # This format is pretty complicated, but also powerful. 
    vargs = [((𝐊!, false, true),  (rgb = PALETTE_GRGB[3], r = 3f0, strength = 0.4f0)),
            ((𝐊!, false, false), (rgb = PALETTE_GRGB[4], r = 2f0, strength = 0.4f0))]
    grid_fcall_with_background(plot_bidirec_streamlines!, vargs; Δ = 100, offset = (-30, -30))
end

# dtmax fixes errors where a trace loops back
# 0.888065 seconds (9.56 M allocations: 1.749 GiB, 60.76% gc time, 1 lock conflict)
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




