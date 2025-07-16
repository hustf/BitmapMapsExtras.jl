using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: paint_curvature_types!

include("common.jl")

#################################
# Curvature types, graphical test
#################################


f = (img, z, grpts) -> paint_curvature_types!(img, z, grpts;
    maxcurv_flat = 0.00002)

# All types
img = grid_fcall_with_background(; f, z = z_ridge_peak_valleys(), Δ = 1)
@test hashstr(img) == "98ed0c914b98ebda87efa0ba3cb58f7b53fed303" 

# All convex (red)
img = grid_fcall_with_background(; f, z = z_paraboloid(; a= 0.6r, b = 0.5r), Δ = 1)
@test hashstr(img) == "8f33c5a74801dbe4bbd9b1e8a3f909fa6507aa91"

# All concave  (green) or flat (gray). Imperfect at edges.
img = grid_fcall_with_background(; f, z = z_ellipsoid(), Δ = 1)
@test let
    h = hashstr(img) 
    h == "a73a45c827d6b4f3036c4bc581e1905e60b93f02" || h == "6844731cbacad0dd802454a3db5058db5f63d243"
end

# All convex-concave, saddle (blue)
img = grid_fcall_with_background(; f, z = z_paraboloid(;a= 0.6r, b = -0.4r), Δ = 1)
@test hashstr(img) == "4f4ecc0d7f08ce57bdb29e814a41760c60a5ff6e"

# A cylinder. Concave.
img = grid_fcall_with_background(; f, z = z_cylinder(1), Δ = 1)
@test hashstr(img) == "9113d7dcda9f5845e3e924419b7e49cea7a2fc0b"

# Lower part of a cylinder. Convex
img = grid_fcall_with_background(; f, z = -z_cylinder(π / 6), Δ = 1)
@test hashstr(img) == "911d663925b9c0b156f82b432c4018c43b9f6212"

