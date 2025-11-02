using Test
using BitmapMapsExtras: N0f8, RGBA
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: plot_glyphs!, plot_glyphs, PALETTE_GRGB, RGB
using BitmapMapsExtras: GSTensor, GSTangentBasis, GSVector
using BitmapMapsExtras: Vec2OnGrid, BidirectionOnGrid, 𝐧ₚᵤ!
using BitmapMapsExtras: indices_scattered, indices_on_grid
using BitmapMapsExtras: radial_distance_glyph, pack_glyphs!
using BitmapMapsExtras: indices_on_grid, default_ij_functor, 𝐊ᵤ!, 𝐊!

!@isdefined(hash_image) && include("common.jl")

@testset "Plot tangent basis" begin
    vhash = ["9b85d71b5edfe55a284bbb4907b67a2809ab931e", "04f14ea5f72f7e56cd3179252b5047f3e4011f55"]
    COUNT[] = 0
    gs = GSTangentBasis()
    pts = [CartesianIndex((315, 215))]
    img = plot_glyphs(z_paraboloid(), pts, gs)
    @test is_hash_stored(img, vhash)
    img = background(z_paraboloid())
    pts = indices_on_grid(size(img))
    plot_glyphs!(img, z_paraboloid(), pts, gs)
    @test is_hash_stored(img, vhash)
end


@testset "Descent vector" begin 
    vhash = ["fbfd0a976f81d605a7a67fe9e2456ce691ef3d77", "ca4b35dc0aa0897e34a30e24b91b8b1160a5ff79"]
    COUNT[] = 0
    gs = GSVector()
    pts = [CartesianIndex((315, 215))]
    # This creates the default ij functor based on z and gs: With function 𝐧ₚ!.
    img = plot_glyphs(z_paraboloid(), pts, gs)
    @test is_hash_stored(img, vhash)
    img = background(z_paraboloid())
    pts = indices_on_grid(size(img))
    plot_glyphs!(img, z_paraboloid(), pts, gs) 
    @test is_hash_stored(img, vhash)
end

@testset "Descent unit vector" begin 
    vhash = ["6ddcfde16c64e1736db9a2c65c7981aac4369c3d", "3e238347523804239ee75c33a9409876ef38df6a"]
    COUNT[] = 0
    gs = GSVector()
    pts = [CartesianIndex((315, 215))]
    # 𝐧ₚᵤ! is not the default for a GSVector glyph spec, so we construct 
    # the AbstractIJFunctor here:
    vog = Vec2OnGrid(𝐧ₚᵤ!, z_paraboloid())
    img = plot_glyphs(vog, pts, gs)
    @test is_hash_stored(img, vhash)
    img = background(vog)
    pts = indices_on_grid(vog)
    plot_glyphs!(img, vog, pts, gs)
    @test is_hash_stored(img, vhash)
end


@testset "Curvature" begin
    vhash = ["25461ae74fe83cd2903aadc09368ead0add36d29", "e798766ebd9238e12908aee6f2701411d2986851"]
    COUNT[] = 0
    # We don't make a functor object, because the
    # default for a GSTensor is what we want.
    gs = GSTensor( multip = 12000)
    pts = [CartesianIndex((315, 215))]
    img = plot_glyphs(z_paraboloid(), pts, gs)
    @test is_hash_stored(img, vhash)
    img = background(z_paraboloid())
    pts = indices_on_grid(size(img))
    plot_glyphs!(img, z_paraboloid(), pts, gs)
    @test is_hash_stored(img, vhash)
end

@testset "Unit curvature" begin
    vhash = ["33aa58431f0e539d331d1d800da9755809e653c4", "1905f3f4c98a772c14e4c7f09c0037baea415477"]
    COUNT[] = 0
    gs = GSTensor( multip = 45)
    pts = [CartesianIndex((315, 215))]
    # 𝐊ᵤ! is not the default for a GSTensor glyph spec, so we construct 
    # the AbstractIJFunctor here:
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_paraboloid())
    img = plot_glyphs(bdog, pts, gs)
    @test is_hash_stored(img, vhash)
    img = background(bdog)
    pts = indices_on_grid(bdog)
    plot_glyphs!(img, bdog, pts, gs)
    @test is_hash_stored(img, vhash)
end


@testset "Packed" begin
    vhash = ["43a0f15aec21deae2c85024d7c28ed44221329ea", "19546f7b18ea588a044b65b39e76287b16283e70", "983b300e6bb0dd06142b1f30a780d08dc483efcd", "7a28f4ecdcc9e3ee202c615bff6211e5e453e14a", "7b15c6c74f22881e7da8618b21e7d0a90eeb5478", "fadbbc288db34d0cc6c8c9bc517030d9068b76e1", "1346ecb69b5de921cc40184735c16f82db814956"]
    COUNT[] = 0
    z = z_ellipsoid(;tilt = 0.5)[100:400, 100:400]
    #
    img = background(z)
    pack_glyphs!(img, z,  GSTangentBasis(; halfsize = 15))
    @test is_hash_stored(img, vhash)
    # Packed projected normals
    img = background(z)
    pack_glyphs!(img, z, 
        GSVector(multip = 200, maxg = 200 * 0.2, color = PALETTE_GRGB[4]),
        scatterdist = 4)
    @test is_hash_stored(img, vhash)
    # Packed unit descent
    vog = Vec2OnGrid(𝐧ₚᵤ!, z)
    img = background(vog)
    pack_glyphs!(img, vog, 
        GSVector(multip = 20, maxg = 200 * 0.2, color = PALETTE_GRGB[4]),
        scatterdist = 1)
    @test is_hash_stored(img, vhash)
    # Packed curvature
    img = background(z)
    pack_glyphs!(img, z, 
        GSTensor(multip = 12000, ming = -25, colors=(RGB(0.0,0.049,0.0), RGB(0.0,0.443,1.0))))
    @test is_hash_stored(img, vhash)
    # Packed unit curvature
    bdog = BidirectionOnGrid(𝐊ᵤ!, z)
    img = background(bdog)
    pack_glyphs!(img, bdog, 
        GSTensor(multip = 20, ming = -25, colors=(RGB(0.0,0.049,0.0), RGB(0.0,0.443,1.0))),
        scatterdist = 1)
    @test is_hash_stored(img, vhash)
    # Primary curvature direction (most positive value). 
    img = background(z)
    pack_glyphs!(img, z, 
        GSTensor(multip = 12000, ming = -25, direction = 1,
        colors=(RGB(0.0,0.049,0.0), RGB(0.0,0.443,1.0)))) 
    @test is_hash_stored(img, vhash)
    # Secondary curvature direction (most positive value).
    img = background(z) 
    pack_glyphs!(img, z, 
        GSTensor(multip = 12000, ming = -25, direction = 2, 
            colors=(RGB(0.0,0.443,1.0), RGB(0.0,0.443,1.0))))
    @test is_hash_stored(img, vhash)
end

@testset "More examples curvature" begin
    vhash = ["dda981106fbdba38ccae4d6004cf628d733054e3", "b006ef40610d3be1e13f0dd5a66c016be775860f", "87b953a36bf264576bbf5c4c209ebb040bcac001", "074614f13d7b76f9b8364f3e45ab1de4ef10aa11", "3fe83b80336530f288bb62f1d920fa95c992f57d", "70088b51485d4410647aaf54117aa809928bae9a", "77b148374b5ee437f44aa4ac6213e02cee8b7c87", "2a3322a5b2a43d47fa52b56cf078d7a2fea52ac1"]
    COUNT[] = 0
    bdog = BidirectionOnGrid(𝐊!, z_cylinder_offset(π / 6))
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 12000))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊!, z_ellipsoid(; tilt = π / 4))
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 12000, strength = 10))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊!, z_ellipsoid(; tilt = π / 4, a = 0.3))
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 12000, strength = 10))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊!, z_paraboloid())
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 12000))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊!, z_paraboloid(; a = 400, b = 600))
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 10000))
    @test is_hash_stored(img, vhash)
    # 
    bdog = BidirectionOnGrid(𝐊!, z_exp3())
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 7000, maxg = 100, ming=-100))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊!, z_cos(;mult = 50))
    img = background(bdog; α = 0.3)
    pack_glyphs!(img, bdog, GSTensor(multip = 17000))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊!, z_ridge_peak_valleys())
    img = background(bdog; α = 0.3)
    pack_glyphs!(img, bdog, GSTensor(multip = 5000))
    @test is_hash_stored(img, vhash)
end


@testset "More examples unit curvature" begin
    vhash = ["6662ab9dd72373e6318f77260971902ab9babb5c", "94356cf19e013fc63a1d29b61e98f0ad98f75d7e", "0c3e2e07d1e86952c627ea9ca4d536ed8295ac89", "f8a967162188d9733aabced5ac1fbd6b399294bd", "675b8ef76cee851a263a9c594de9d576c1610dcb", "4ea93a4ea55472c3b490fcefe8ecd01999b5c892", "3917fe66c32271ead12ba7e2ff2422a78e5149f4", "5ee88c4a619397e3d672eb59be610de5c1e6ad0f"]
    COUNT[] = 0
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_cylinder_offset(π / 6))
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor())
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_ellipsoid(; tilt = π / 4))
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 30))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_ellipsoid(; tilt = π / 4, a = 0.3))
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 30))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_paraboloid())
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor())
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_paraboloid(; a = 400, b = 600))
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 10000))
    @test is_hash_stored(img, vhash)
    # 
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_exp3())
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 20, maxg = 100, ming=-100))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_cos(;mult = 50))
    img = background(bdog; α = 0.3)
    pack_glyphs!(img, bdog, GSTensor(multip = 20))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_ridge_peak_valleys())
    img = background(bdog; α = 0.3)
    pack_glyphs!(img, bdog, GSTensor(;multip = 15, strength = 4))
    @test is_hash_stored(img, vhash)
end

# DEV 

# This aims to benchmark the drawing functionality itself.
#=
using BenchmarkTools # test dependency
using BitmapMapsExtras: indices_scattered, MersenneTwister
using BitmapMapsExtras: placements_and_values, plot_glyphs_given_values!
!@isdefined(hash_image) && include("common.jl")

vhash = ["943f92be5048f710722304cbcca17b60e20bc2b2"]
COUNT[] = 0
#
bdog = BidirectionOnGrid(𝐊!, z_ellipsoid(; tilt = π / 4))
img = background(bdog);
display_if_vscode(img)

# We're unpacking this call:
# pack_glyphs!(img, bdog, GSTensor(multip = 12000, strength = 10))
scatterdist = 3.0
seed = MersenneTwister(123)
fij = bdog
gs = GSTensor(multip = 12000, strength = 10)
ppts = indices_scattered(fij; scatterdist, seed)
filtered_placements, filtered_values = placements_and_values(fij, gs, ppts)
# Now benchmark the plotting for later reference
# 33.181 ms (5480 allocations: 7.82 MiB) # Before moving funcs to `DrawAndSpray`
# 27.076 ms (5394 allocations: 15.43 MiB) # After moving funcs to 'DrawAndSpray'
@btime plot_glyphs_given_values!(img, filtered_placements, filtered_values, gs)
@test is_hash_stored(img, vhash)
@profview plot_glyphs_given_values!(img, filtered_placements, filtered_values, gs)
=#