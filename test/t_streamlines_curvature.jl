# This is not successful at this point

using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras.TestMatrices: I0
using BitmapMapsExtras: plot_streamlines!, principal_curvature_components
using BitmapMapsExtras: 𝐧ₚ!, DirectionAtXY, is_bidirec_vect_positive, norm

include("common.jl")
include("dev_streamlines_curvature.jl")


let # For a cylinder, max principal curvature is zero, and min principal is negative
    f = (img, z, grpts) -> plot_streamlines!(img, c11!(), z, grpts; tstop = 100000.0, r = 1.2)
    grid_fcall_with_background(;f, z = z_cylinder(0.5), Δ = 150)
end

let # For a cylinder, max principal curvature is zero, and min principal is negative
    f = (img, z, grpts) -> plot_streamlines!(img, c12!(), z, grpts; tstop = 100000.0, r = 1.2)
    grid_fcall_with_background(;f, z = z_cylinder(0.5), Δ = 150)
end

let # For a cylinder, max principal curvature is zero, and min principal is negative
    f = (img, z, grpts) -> plot_streamlines!(img, c21!(), z, grpts; tstop = 100000.0, r = 1.2)
    grid_fcall_with_background(;f, z = z_cylinder(0.5), Δ = 150)
end

let # For a cylinder, max principal curvature is zero, and min principal is negative
    f = (img, z, grpts) -> plot_streamlines!(img, c22!(), z, grpts; tstop = 100000.0, r = 1.2)
    grid_fcall_with_background(;f, z = z_cylinder(0.5), Δ = 150)
end


let # For a sphere, curvature is theoretically equal in any direction. This result is round-off errors
    f = (img, z, grpts) -> plot_streamlines!(img, c11!(), z, grpts; tstop = 100000.0, r = 1.2)
    grid_fcall_with_background(;f, z = z_sphere(), Δ = 150)
end

let # For a paraboloid, curvature varies with x and y
    f = (img, z, grpts) -> plot_streamlines!(img, c11!(), z, grpts; tstop = 100000.0, r = 1.2)
    grid_fcall_with_background(;f, z = z_paraboloid(a=0.5r, b=-r), Δ = 150)
end

let # For a paraboloid, curvature varies with x and y
    f = (img, z, grpts) -> plot_streamlines!(img, c12!(), z, grpts; tstop = 100000.0, r = 1.2)
    grid_fcall_with_background(;f, z = z_paraboloid(a=0.5r, b=-r), Δ = 150)
end

let # Draw streamlines in both directions from the grid (but along the same curvature line)
    f1 = (img, z, grpts) -> plot_streamlines!(img, c11!(), z, grpts; tstop = 500000, r = 1.2)
    f2 = (img, z, grpts) -> plot_streamlines!(img, c12!(), z, grpts; tstop = 500000, r = 1.2)
    grid_fcall_with_background(;f = [f1, f2], z = z_paraboloid(a=0.5r, b=-r), Δ = 150)
end

let # Draw streamlines in both directions from the grid (but along the same curvature line)
    f1 = (img, z, grpts) -> plot_streamlines!(img, c21!(), z, grpts; tstop = 500000, r = 1.2)
    f2 = (img, z, grpts) -> plot_streamlines!(img, c22!(), z, grpts; tstop = 500000, r = 1.2)
    grid_fcall_with_background(;f = [f1, f2], z = z_paraboloid(a=0.5r, b=-r), Δ = 150)
end


let # For an ellipsoid, lines of curvature have cusps 
    f1 = (img, z, grpts) -> plot_streamlines!(img, c11!(), z, grpts; tstop = 500000, r = 1.2)
    f2 = (img, z, grpts) -> plot_streamlines!(img, c12!(), z, grpts; tstop = 500000, r = 1.2)
    grid_fcall_with_background(;f = [f1, f2], z = z_ellipsoid(; a = 0.49), Δ = 50)
end

let # For an ellipsoid, lines of curvature have cusps 
    f1 = (img, z, grpts) -> plot_streamlines!(img, c21!(), z, grpts; tstop = 500000, r = 1.2)
    f2 = (img, z, grpts) -> plot_streamlines!(img, c22!(), z, grpts; tstop = 500000, r = 1.2)
    grid_fcall_with_background(;f = [f1, f2], z = z_ellipsoid(; a = 0.49), Δ = 50)
end



let # 
    f1 = (img, z, grpts) -> plot_streamlines!(img, c11!(), z, grpts; tstop = 2000000, r = 1.2)
    f2 = (img, z, grpts) -> plot_streamlines!(img, c12!(), z, grpts; tstop = 2000000, r = 1.2)
    grid_fcall_with_background(;f = [f1, f2], z = z_ridge_peak_valleys(), Δ = 50)
end

let # 
    f1 = (img, z, grpts) -> plot_streamlines!(img, c21!(), z, grpts; tstop = 10000000, r = 1.2)
    f2 = (img, z, grpts) -> plot_streamlines!(img, c22!(), z, grpts; tstop = 10000000, r = 1.2)
    grid_fcall_with_background(;f = [f1, f2], z = z_ridge_peak_valleys(), Δ = 50)
end

# As expected, the lines of curvature do not follow where curvature is greatest.
# In order to trace out lines of maximal curvature, one would have to travel perpendicularly
# to those, but even better would be to make a parameter  (maximal curvature - minimal curvature) 
# and then follow the direction of least change.
# In practice, that's close to what BitmapMaps already does.




##################
# Continuity issue
##################

let # Draw streamlines in both directions from the grid (but along the same curvature line)
    f1 = (img, z, grpts) -> plot_streamlines!(img, c11!(), z, grpts; tstop = 2^20, r = 1.2)
    f2 = (img, z, grpts) -> plot_streamlines!(img, c12!(), z, grpts; tstop = 5, r = 1.2)
    z = z_paraboloid(a = 0.5r, b = -r)
    img = background(z)
    pts = CartesianIndices((300:300, -1:1)) .+ I0
    #f1(img, z, pts)
    pt = pts[1]
    v = [0.0, 0.0]
    Ω = CartesianIndices((-2:2, -2:2))
    c11!()(v, view(z, Ω .+ pts[1]))
    c11!()(v, view(z, Ω .+ pts[2]))
    ### Hmm, maybe we should just stop lines at zero-crossing curvatures?
end;
