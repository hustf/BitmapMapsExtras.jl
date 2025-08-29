using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: paint_curvature_types!

include("common.jl")

#################################
# Curvature types, graphical test
#################################
# TODO: Switch to new 'grid_fcall_with_background' format.

@testset "curvature types" begin

r = TestMatrices.r

f = (img, z, grpts) -> paint_curvature_types!(img, z, grpts;
    maxcurv_flat = 0.00002)

# All types
img = grid_fcall_with_background(; f, z = z_ridge_peak_valleys(), Δ = 1)
@test hashstr(img) == "98ed0c914b98ebda87efa0ba3cb58f7b53fed303" 

# All convex (red)
img = grid_fcall_with_background(; f, z = z_paraboloid(; a= 0.6r, b = 0.5r), Δ = 1)
@test hashstr(img) == "8f33c5a74801dbe4bbd9b1e8a3f909fa6507aa91"

# All concave  (green)
# Slight difference between calculating in REPL and in test mode
img = grid_fcall_with_background(; f, z = z_ellipsoid(), Δ = 1)
hashstr(img) == "8e144c9b3fd1ed1279f9afec9c9eba21c64bcf41"

# All convex-concave, saddle (blue)
img = grid_fcall_with_background(; f, z = z_paraboloid(;a= 0.6r, b = -0.4r), Δ = 1)
@test hashstr(img) == "4f4ecc0d7f08ce57bdb29e814a41760c60a5ff6e"

# A cylinder. Concave.
img = grid_fcall_with_background(; f, z = z_cylinder(1), Δ = 1)
@test hashstr(img) == "9113d7dcda9f5845e3e924419b7e49cea7a2fc0b"

# Lower part of a cylinder. Convex
img = grid_fcall_with_background(; f, z = -z_cylinder(π / 6), Δ = 1)
@test hashstr(img) == "911d663925b9c0b156f82b432c4018c43b9f6212"

end