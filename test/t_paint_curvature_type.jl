using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: paint_curvature_types!, indices_on_grid

!@isdefined(is_hash_stored) && include("common.jl")

#################################
# Curvature types, graphical test
#################################

@testset "Curvature types" begin
    vhash = ["9f43a178e18d8b504e09743f22d40babfaafd7ea", "b2a5b2e2f59ffefb2bd2e141bd3ac7cbe168cae1", "6d347d6ef8a5258fa49dac50ea9e94ed7087afbe", "2807b28535986e7bbb42584ff73131c032c10aba", "1785a9c9f6c20885b5f92eaad88d34093820dbd0", "c1d2c54b6cedca7ed64977631d3ed45219e35624"]
    COUNT[] = 0
    z = z_ridge_peak_valleys()
    pts = CartesianIndices(z)
    img = background(z)
    paint_curvature_types!(img, z, pts)
    @test is_hash_stored(img, vhash)

    # All convex (red)
    r = TestMatrices.r
    z .= z_paraboloid(; a = 0.6r, b = 0.5r)
    img .= background(z)
    paint_curvature_types!(img, z, pts)
    @test is_hash_stored(img, vhash)

    # All concave  (green)
    z .= z_ellipsoid()
    img .= background(z)
    paint_curvature_types!(img, z, pts)
    @test is_hash_stored(img, vhash)

    # All convex-concave, saddle (blue)
    z .= z_paraboloid(;a= 0.6r, b = -0.4r)
    img .= background(z)
    paint_curvature_types!(img, z, pts)
    @test is_hash_stored(img, vhash)

    # A cylinder. Concave (green).
    z .= z_cylinder(1)
    img .= background(z)
    paint_curvature_types!(img, z, pts)
    @test is_hash_stored(img, vhash)

    # Lower part of a cylinder. Convex (red)
    z .= -z_cylinder(π / 6)
    img .= background(z)
    paint_curvature_types!(img, z, pts)
    @test is_hash_stored(img, vhash)
end