using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: paint_curvature_types!, indices_on_grid, RGB
using BitmapMapsExtras: Random

!@isdefined(is_hash_stored) && include("common.jl")

#################################
# Curvature types, graphical test
#################################

@testset "Curvature types" begin
    vhash = ["217360f5027e03457d1dcdad29ed973d93b8e793", "cbe32ab9d7b405866e8200671817e68eca09e705", "4ee47227650446059979bfff36fa1912f1c7c30c", "37b1b89afdf83f220128cf39eeb8e6de6b2d9a63", "634b828141731d70f784615ec44d71981e3ceaab", "1ef7d86da5b6a5de02cc226221904aacaa43a6ae"]
    COUNT[] = 0
    z = z_ridge_peak_valleys()
    pts = CartesianIndices(z)
    img = RGB.(background(z))
    paint_curvature_types!(img, z, pts)
    @test is_hash_stored(img, vhash)
    # All convex (red)
    r = TestMatrices.r
    z .= z_paraboloid(; a = 0.6r, b = 0.5r)
    img .= RGB.(background(z))
    paint_curvature_types!(img, z, pts)
    is_hash_stored(img, vhash)
    # All concave  (green)
    z .= z_ellipsoid()
    img .= RGB.(background(z))
    paint_curvature_types!(img, z, pts)
    is_hash_stored(img, vhash)
    # All convex-concave, saddle (blue)
    z .= z_paraboloid(;a= 0.6r, b = -0.4r)
    img .= RGB.(background(z))
    is_hash_stored(img, vhash)
    paint_curvature_types!(img, z, pts)
    # A cylinder. Concave (green).
    z .= z_cylinder(1)
    img .= RGB.(background(z))
    paint_curvature_types!(img, z, pts)
    is_hash_stored(img, vhash)
    # Lower part of a cylinder. Convex (red)
    z .= -z_cylinder(π / 6)
    img .= RGB.(background(z))
    paint_curvature_types!(img, z, pts)
    is_hash_stored(img, vhash)
end
