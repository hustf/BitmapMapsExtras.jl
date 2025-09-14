using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: paint_curvature_types!

!@isdefined(hashstr) && include("common.jl")

#################################
# Curvature types, graphical test
#################################

@testset "curvature types" begin

r = TestMatrices.r
f = paint_curvature_types!
# All types
args = (; maxcurv_flat = 0.00002)
img = grid_fcall_with_background(f, args; z = z_ridge_peak_valleys(), Δ = 1)
@test hashstr(img) == "4b908a01d21408e01dcf93736fd0df4374031e88"

# All convex (red)
args = (;)
img = grid_fcall_with_background(f, args; z = z_paraboloid(; a = 0.6r, b = 0.5r), Δ = 1)
@test hashstr(img) == "a020b8a5f472f796cfd3fcbefe7e405a2cfd3d9b"

# All concave  (green)
img = grid_fcall_with_background(f, args; z = z_ellipsoid(), Δ = 1)
@test hashstr(img) == "8a968074801f2a0619d5c74a66072177be5c0a01" 

# All convex-concave, saddle (blue)
img = grid_fcall_with_background(f, args; z = z_paraboloid(;a= 0.6r, b = -0.4r), Δ = 1)
@test hashstr(img) == "081a63c8c3177dd948541a10481372941ce6250f"

# A cylinder. Concave (green).
img = grid_fcall_with_background(f, args; z = z_cylinder(1), Δ = 1)
@test hashstr(img) == "89d01a131a2e1bf6e66600850d89b8b2d4ac40c9"

# Lower part of a cylinder. Convex (red)
img = grid_fcall_with_background(f, args; z = -z_cylinder(π / 6), Δ = 1)
@test hashstr(img) == "d52953766db89bf022185a0c76ba19d5f7cd0939"

end