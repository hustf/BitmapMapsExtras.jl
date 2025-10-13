using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: plot_glyphs!, plot_glyphs, PALETTE_GRGB
using BitmapMapsExtras: GSTensor, GSTangentBasis, GSVector
using BitmapMapsExtras: indices_scattered, indices_on_grid
using BitmapMapsExtras: radial_distance_glyph, pack_glyphs!
using BitmapMapsExtras: indices_on_grid



!@isdefined(hashstr) && include("common.jl")

@testset "Plot tangent basis" begin 
    gs = GSTangentBasis()
    pts = [CartesianIndex((315, 215))]
    img = plot_glyphs(z_paraboloid(), pts, gs)
    @test hashstr(img) == "49d9895e4a2f0babed325677c10adf5c17669f35"
    img = background(z_paraboloid())
    pts = indices_on_grid(size(img))
    plot_glyphs!(img, z_paraboloid(), pts, gs)
    @test hashstr(img) == "44dbb26b32cddae24eed94b18fbce0fc88917081"
end

@testset "Projected normal vector" begin 
    gs = GSVector()
    pts = [CartesianIndex((315, 215))]
    img = plot_glyphs(z_paraboloid(), pts, gs)
    @test hashstr(img) == "0f637ec62c13f4e449c60c6427489a8de62f1cbd"
    img = background(z_paraboloid())
    pts = indices_on_grid(size(img))
    plot_glyphs!(img, z_paraboloid(), pts, gs)
    @test hashstr(img) == "df4f6aee96f0353ce73045994cc5d0f0bf791cb1"
end

@testset "Curvature" begin 
    gs = GSTensor( multip = 20000, strength = 1.0)
    pts = [CartesianIndex((315, 215))]
    img = plot_glyphs(z_paraboloid(), pts, gs)
    @test hashstr(img) == "6b1ca65d6cc8cf2be9091be3b0733489cbb5ba75"
    img = background(z_paraboloid())
    pts = indices_on_grid(size(img))
    plot_glyphs!(img, z_paraboloid(), pts, gs)
    @test hashstr(img) == "90ae945fe278f2662b5b68e4c61bf23d1640d652"
end

@testset "Packed" begin
    z = z_ellipsoid(;tilt = 0.5)[100:400, 100:400]
    #
    img = pack_glyphs_with_background(z,  GSTangentBasis(; halfsize = 15))
    @test hashstr(img) == "82be1509834d02e89b1d762db5c10cf675348fa2"
    # Packed projected normals
    img = pack_glyphs_with_background(z, 
        GSVector(multip = 200, maxg = 200 * 0.2, color = PALETTE_GRGB[4]),
        scatterdist = 4)
    @test hashstr(img) == "38b4d02cd1d97e4125de8dc2fedb4794ece3a96d"
    # Packed curvature
    img = pack_glyphs_with_background(z, 
        GSTensor(multip = 12000, ming = -25, colors=(RGB(0.0,0.049,0.0), RGB(0.0,0.443,1.0))))
    @test hashstr(img) == "35cbb53f00144a777b25994a87a78b9506b045d4"
    # Primary curvature direction (most positive value). 
    img = pack_glyphs_with_background(z, 
        GSTensor(multip = 12000, ming = -25, direction = 1,
        colors=(RGB(0.0,0.049,0.0), RGB(0.0,0.443,1.0)))) 
    @test hashstr(img) == "c2d63dbe0114a0c3a924f6927891fcdd24b950f2"
    # Secondary curvature direction (most positive value). 
    img = pack_glyphs_with_background(z, 
        GSTensor(multip = 12000, ming = -25, direction = 2, 
            colors=(RGB(0.0,0.443,1.0), RGB(0.0,0.443,1.0))))
    @test hashstr(img) == "9d9745bd3b44b8e33627ef5e4df34b1771f63e88"
end

@testset "More examples" begin
    pack_glyphs_with_background(z_cylinder_offset(π / 6), GSTensor(multip = 12000))
    gs = GSTensor(multip = 5000, ming = -25)
    pack_glyphs_with_background(z_ellipsoid(; tilt = π / 4), gs; scatterdist = 1)
    gs = GSTensor(multip = 5000, strength = 10, ming = -25)
    pack_glyphs_with_background(z_ellipsoid(; tilt = π / 4, a= 0.3), gs; scatterdist = 2)
    pack_glyphs_with_background(z_paraboloid(), GSTensor(multip = 15000))
    pack_glyphs_with_background(z_paraboloid(; a = 400, b = 600), GSTensor(multip = 10000 ))
    pack_glyphs_with_background(z_sphere(), GSTensor(multip = 10000))
    pack_glyphs_with_background(z_exp3(), GSTensor(multip = 7000))
    pack_glyphs_with_background(z_cos(; mult = 50), GSTensor(multip = 25000))
    pack_glyphs_with_background(z_ridge_peak_valleys(), GSTensor(multip = 4000 ), scatterdist = 2)
    @test true
end

