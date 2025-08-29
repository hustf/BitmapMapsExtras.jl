using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras.TestMatrices: I0
using BitmapMapsExtras: plot_bidirec_streamlines!, principal_curvature_components
using BitmapMapsExtras: 𝐊! #, DirectionAtXY, is_bidirec_vect_positive, norm
using BitmapMapsExtras: PALETTE_GRGB, get_bidirec_streamlines_points
import BitmapMaps
using BitmapMaps: mark_at!

include("common.jl")

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



# Second run. TODO: Droip allocations in tangent basis (𝐊!)
# 0.842606 seconds (5.23 M allocations: 1.470 GiB, 65.65% gc time)
@time let 
    # This format is pretty complicated, but also powerful. 
    vargs = [((𝐊!, false, true),  (rgb = PALETTE_GRGB[3], r = 3f0, strength = 0.4f0)),
            ((𝐊!, false, false), (rgb = PALETTE_GRGB[4], r = 2f0, strength = 0.4f0))]
    grid_fcall_with_background(plot_bidirec_streamlines!, vargs, Δ = 100, centreoffset = -30)
end

# dtmax fixes errors where a trace loops back
# 0.888065 seconds (9.56 M allocations: 1.749 GiB, 60.76% gc time, 1 lock conflict)
@time let 
    vargs = [((𝐊!, false, true),  (rgb = PALETTE_GRGB[3], r = 1.5f0, strength = 0.4f0, dtmax = 1)),
            ((𝐊!, false, false), (rgb = PALETTE_GRGB[4], r = 1.5f0, strength = 0.4f0, dtmax = 1)), 
            ((𝐊!, true, true),  (rgb = PALETTE_GRGB[1], r = 1.5f0, strength = 0.4f0, dtmax = 1)),
            ((𝐊!, true, false), (rgb = PALETTE_GRGB[2], r = 1.0f0, strength = 0.4f0, dtmax = 1))]
    grid_fcall_with_background(plot_bidirec_streamlines!, vargs, Δ = 200)
end


let 
    vargs = [((𝐊!, false, true),  (rgb = PALETTE_GRGB[3], r = 3f0, strength = 0.4f0)),
            ((𝐊!, false, false), (rgb = PALETTE_GRGB[4], r = 2f0, strength = 0.4f0))]
    grid_fcall_with_background(plot_bidirec_streamlines!, vargs, Δ = 100, centreoffset = -30, z = z_cylinder(π / 6))
end

# dtmax fixes undulating traces
let 
    vargs = [((𝐊!, false, true),  (rgb = PALETTE_GRGB[3], r = 3f0, strength = 0.4f0, dtmax = 1)),
            ((𝐊!, false, false), (rgb = PALETTE_GRGB[4], r = 2f0, strength = 0.4f0, dtmax = 1))]
    grid_fcall_with_background(plot_bidirec_streamlines!, vargs, Δ = 100, centreoffset = -30, z = z_cylinder(π / 6))
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
    grid_fcall_with_background(plot_bidirec_streamlines!, vargs, Δ = 100, centreoffset = -30, z = z_paraboloid())
end


# This encounters switch primary <-> secondary
let 
    vargs = [((𝐊!, false, true),  (rgb = PALETTE_GRGB[3], r = 3f0, strength = 0.4f0, dtmax = 1, tstop = 3000))]
    grid_fcall_with_background(plot_bidirec_streamlines!, vargs, Δ = 1500, centreoffset = -20, z = z_ridge_peak_valleys())
end

let 
    vargs = [((𝐊!, false, true),  (rgb = PALETTE_GRGB[3], r = 3f0, strength = 0.4f0, dtmax = 1, tstop = 3000))]
    fcall_with_background(plot_bidirec_streamlines!, vargs, z_ridge_peak_valleys(), !!!)
end



#radius = TestMatrices.r




z = z_cylinder(0.5)
f = (img, z, grpts) -> plot_bidirec_streamlines!(img, 𝐊!, z, grpts, false, false, nsamples = 1000)
grid_fcall_with_background(;f, z, Δ = 500)
f = (img, z, grpts) -> plot_bidirec_streamlines!(img, 𝐊!, z, grpts, false, true, nsamples = 1000)


pts = [CartesianIndex(300, 300)]
@test length(get_bidirec_streamlines_points(𝐊!, z, pts, 1000, true, true)[1]) == 649



# For a paraboloid, curvature varies with x and y.
let 
    f = (img, z, grpts) -> plot_bidirec_streamlines!(img, 𝐊!, z, grpts, false, false)
    grid_fcall_with_background(;f,  Δ = 500, centreoffset = -100)
end

let # For a paraboloid, curvature varies with x and y.
    f = (img, z, grpts) -> begin
        plot_bidirec_streamlines!(img, 𝐊!, z, grpts, false, false)
        plot_bidirec_streamlines!(img, 𝐊!, z, grpts, false, true)
    end
    grid_fcall_with_background(;f,  Δ = 150, centreoffset = -100)
end




let 
    f = [(img, z, grpts) -> plot_bidirec_streamlines!(img, 𝐊!, z, grpts, false, false; rgb = PALETTE_GRGB[2], r = 3f0),
         (img, z, grpts) -> plot_bidirec_streamlines!(img, 𝐊!, z, grpts, false, true; rgb = PALETTE_GRGB[1], r = 8f0)]
    grid_fcall_with_background(;f,  Δ = 150, centreoffset = -100)
end

# TODO: Isolate one successful and one unsuccessful  using fcall_with_background.
# Track down why only the first call calls both.
let 
   f = [(img, z, grpts) -> plot_bidirec_streamlines!(img, 𝐊!, z, grpts, false, false; rgb = PALETTE_GRGB[2], r = 3f0),
         (img, z, grpts) -> plot_bidirec_streamlines!(img, 𝐊!, z, grpts, false, true; rgb = PALETTE_GRGB[1], r = 8f0)]
        # One start points
    Δ = 600
    rngj = 525:Δ:700
    rngi = 225:Δ:925
    startpts = CartesianIndices((rngi, rngj))
    z = z = z_paraboloid()
    fcall_with_background(f, z, startpts)
end






let # For a cylinder, max principal curvature is zero, and min principal is negative
    f = (img, z, grpts) -> plot_bidirec_streamlines!(img, c1a!(), z, grpts )
    grid_fcall_with_background(;f, z = z_cylinder(0.5))
end

let # For a cylinder, max principal curvature is zero, and min principal is negative
    f = (img, z, grpts) -> plot_bidirec_streamlines!(img, c1b!(), z, grpts)
    grid_fcall_with_background(;f, z = z_cylinder(0.5))
end

let # For a cylinder, max principal curvature is zero, and min principal is negative
    f = (img, z, grpts) -> plot_bidirec_streamlines!(img, c2a!(), z, grpts)
    grid_fcall_with_background(;f, z = z_cylinder(0.5))
end

let # For a cylinder, max principal curvature is zero, and min principal is negative
    f = (img, z, grpts) -> plot_bidirec_streamlines!(img, c2b!(), z, grpts)
    grid_fcall_with_background(;f, z = z_cylinder(0.5))
end

let # For a sphere, curvature is theoretically equal in any direction. This result is round-off errors
    # We might want to stop in general when max and min curvatures are too similar.
    f = (img, z, grpts) -> plot_bidirec_streamlines!(img, c1a!(), z, grpts)
    grid_fcall_with_background(;f, z = z_sphere())
end

#######################
# Overcome oscillations
#######################

@time let # For a paraboloid, curvature varies with x and y.
    f = (img, z, grpts) -> plot_bidirec_streamlines!(img, 𝐊!, z, grpts, false)
    grid_fcall_with_background(;f)
end

let # Adjusting parameter minc  
    f = (img, z, grpts) -> plot_bidirec_streamlines!(img, c1a!(; minc = 0.002), z, grpts)
    grid_fcall_with_background(;f)
end

@time let # Adjusting parameter dtfloor
    f = (img, z, grpts) -> plot_bidirec_streamlines!(img, c2a!(), z, grpts; dtfloor = 4.0)
    grid_fcall_with_background(;f)
end

@time let # Adjusting parameter dtmax slows things down
    f = (img, z, grpts) -> plot_bidirec_streamlines!(img, c2a!(), z, grpts; dtmax = 1.0)
    grid_fcall_with_background(;f)
end


#####################
# Difference between 
# c1a! c1b! c2a! c2b!
#####################


let
    f = (img, z, grpts) -> plot_bidirec_streamlines!(img, c1a!(), z, grpts)
    grid_fcall_with_background(;f)
end

let # Draw streamlines in all principal directions from the grid 
    f1 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c1a!(), z, grpts; tstop = 50, rgb = PALETTE_GRGB[1], r = 3f0)
    f2 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c1b!(), z, grpts, tstop = 50, rgb = PALETTE_GRGB[2], r = 3f0)
    f3 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c2a!(), z, grpts, tstop = 30, rgb = PALETTE_GRGB[3], r = 3f0)
    f4 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c2b!(), z, grpts, tstop = 30, rgb = PALETTE_GRGB[4], r = 3f0)
    grid_fcall_with_background(;f = [f1, f2, f3, f4])
end

########################################
# Validate towards known curvature lines
########################################

@time let # For an ellipsoid, lines of curvature are visually known.
    # Here, the accuracy is good enough to complete one ellipisis (red),
    # but not in the other direction (grey). Setting dtmax fixes that issue.
    dtmax = 25.0
    tstop = 1500
    r = 3f0
    f1 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c1a!(), z, grpts; tstop, dtmax, r, rgb = PALETTE_GRGB[1])
    f2 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c1b!(), z, grpts; tstop, r, rgb = PALETTE_GRGB[2])
    f3 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c2a!(), z, grpts; tstop, r, rgb = PALETTE_GRGB[3])
    f4 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c2b!(), z, grpts; tstop, r, rgb = PALETTE_GRGB[4])
    # One start points
    Δ = 5000
    rngj = 325:Δ:700
    rngi = 625:Δ:925
    startpts = CartesianIndices((rngi, rngj))
    fcall_with_background([f1, f2, f3, f4], z_ellipsoid(; a = 0.49), startpts)
end

@time let # For an ellipsoid, lines of curvature are visually known.
    # Here, the accuracy is good enough to complete one ellipisis (red),
    # but not in the other direction (grey). Setting dtmax fixes that issue.
    dtmax = 1
    tstop = 1000
    r = 3f0
    f1 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c1a!(), z, grpts; tstop, dtmax, r, rgb = PALETTE_GRGB[1])
    f2 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c1b!(), z, grpts; tstop, r, rgb = PALETTE_GRGB[2])
    f3 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c2a!(), z, grpts; tstop, r, rgb = PALETTE_GRGB[3])
    f4 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c2b!(), z, grpts; tstop, r, rgb = PALETTE_GRGB[4])
    # Start points
    Δ = 50
    rngj = 450:Δ:450
    rngi = 800:Δ:950
    startpts = CartesianIndices((rngi, rngj))
    fcall_with_background([f1, f2, f3, f4], z_ellipsoid(; a = 0.49), startpts)
end






let # Follow the max curvature c1a and c1b - collects mostly in ridges 
    # Follow the minimum curvature c2 a b - collects mostly in dieders 
    dtmax = 10.0
    f1 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c1a!(), z, grpts; dtmax, rgb = PALETTE_GRGB[3], strength = 0.1f0)
    f2 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c1b!(), z, grpts; dtmax, rgb = PALETTE_GRGB[3], strength = 0.1f0)
    f3 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c2a!(), z, grpts; dtmax, rgb = PALETTE_GRGB[4], strength = 0.1f0)
    f4 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c2b!(), z, grpts; dtmax, rgb = PALETTE_GRGB[4], strength = 0.1f0)
    grid_fcall_with_background(;f = [f1, f2, f3, f4], z = z_ridge_peak_valleys(), Δ = 75)
end

let # Follow the max curvature c1a and c1b - collects mostly in ridges 
    # Follow the minimum curvature c2 a b - collects mostly in dieders 
    dtmax = 1.0
    r = 5f0
    tstop = 2000
    kws = (;dtmax, r, strength = 0.1f0, tstop)
    f1 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c1a!(), z, grpts; rgb = PALETTE_GRGB[3], kws...)
    f2 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c1b!(), z, grpts; rgb = PALETTE_GRGB[3], kws...)
    f3 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c2a!(), z, grpts; rgb = PALETTE_GRGB[4], kws...)
    f4 = (img, z, grpts) -> plot_bidirec_streamlines!(img, c2b!(), z, grpts; rgb = PALETTE_GRGB[4], kws...)
    pts = rand(TestMatrices.R, 1)
    fcall_with_background([f1, f2, f3, f4], z_ridge_peak_valleys(), pts)
end







#=
#####################################
# Unpack the level solve_ensemble....
#####################################
using BitmapMapsExtras: UnidirectionAtXY, ODEProblem, EnsembleProblem, Tsit5, EnsembleThreads, EnsembleSerial
using BitmapMapsExtras: ContinuousCallback, DiscreteCallback, signed_distance_within_domain, affect!
using BitmapMapsExtras: too_flat, affect_flip_bidirection!, condition_flip_bidirection, CallbackSet
using BitmapMapsExtras: rhs!, solve, reset!, remake
vx = [0.0 400.0 800.0; 0.0 400.0 800.0; 0.0 400.0 800.0]
vy = [1000.0 1000.0 1000.0; 600.0 600.0 600.0; 200.0 200.0 200.0]
vu0 = [[x, y] for (x, y) in zip(vx, vy)]
tspan = (0.0, 1000)
z = z_paraboloid();
uxy = UnidirectionAtXY(𝐊!, z, false, true)
vccb = [ContinuousCallback(signed_distance_within_domain, affect!),
        ContinuousCallback(too_flat, affect!)]
vdcb = DiscreteCallback[]
push!(vdcb, DiscreteCallback(condition_flip_bidirection, affect_flip_bidirection!, save_positions=(true,true)))
cbs = CallbackSet(vccb..., vdcb...)
# solve_ensemble...

u0 = first(vu0)
function prob_func(prob, i, repeat)
    reset!(prob.p)
    u0 = vu0[i]
    @assert u0 isa Vector
    println("\n", u0, " flip: $(prob.p.flip[]) ...  ")
    remake(prob, u0 = u0)
end
prob = ODEProblem(rhs!, u0, tspan, uxy, callback = cbs)
ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
sim = solve(ensemble_prob, Tsit5(), EnsembleSerial(), trajectories = length(vu0))
=#
