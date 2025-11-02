using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: paint_curvature_types!, indices_on_grid, RGB

!@isdefined(is_hash_stored) && include("common.jl")

#################################
# Curvature types, graphical test
#################################

@testset "Curvature types" begin
    vhash = ["e0ee520e3e3617086df04e9af1f0dfce36c66700", "38a0cad0dde309084a840aaec412854b7061e932", "f00ffa26ac74848e8b2e496723efe9c19b0ee86d", "d95ae038d180f388c71d5a2355c9d9a904b188d1", "c2e621bc8fa02595c10028ec43107dd85c0423c4", "aa82aa6b54496968c7e4b33319ecc5f8b4b99fc3"]
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
