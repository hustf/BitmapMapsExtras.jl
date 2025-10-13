using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: paint_curvature_types!, indices_on_grid

!@isdefined(hashstr) && include("common.jl")

#################################
# Curvature types, graphical test
#################################

@testset "Curvature types" begin
    z = z_ridge_peak_valleys()
    pts = CartesianIndices(z)
    img = background(z)
    paint_curvature_types!(img, z, pts)
    @test hashstr(img) == "b11ae5ad33f32bd1761ee647682b6d8ead446d6f"

    # All convex (red)
    r = TestMatrices.r
    z .= z_paraboloid(; a = 0.6r, b = 0.5r)
    img .= background(z)
    paint_curvature_types!(img, z, pts)
    @test hashstr(img) == "a020b8a5f472f796cfd3fcbefe7e405a2cfd3d9b"

    # All concave  (green)
    z .= z_ellipsoid()
    img .= background(z)
    paint_curvature_types!(img, z, pts)
    @test hashstr(img) == "8a968074801f2a0619d5c74a66072177be5c0a01"

    # All convex-concave, saddle (blue)
    z .= z_paraboloid(;a= 0.6r, b = -0.4r)
    img .= background(z)
    paint_curvature_types!(img, z, pts)
    @test hashstr(img) == "081a63c8c3177dd948541a10481372941ce6250f"


    # A cylinder. Concave (green).
    z .= z_cylinder(1)
    img .= background(z)
    paint_curvature_types!(img, z, pts)
    @test hashstr(img) == "89d01a131a2e1bf6e66600850d89b8b2d4ac40c9"


    # Lower part of a cylinder. Convex (red)
    z .= -z_cylinder(π / 6)
    img .= background(z)
    paint_curvature_types!(img, z, pts)
    @test hashstr(img) == "d52953766db89bf022185a0c76ba19d5f7cd0939"
end