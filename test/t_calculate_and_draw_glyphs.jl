using Test
using BitmapMapsExtras
using BitmapMapsExtras: N0f8, RGBA
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: plot_glyphs!, plot_glyphs, PALETTE_GRGB, RGB
using BitmapMapsExtras: GSTensor, GSTangentBasis, GSVector
using BitmapMapsExtras: Vec2OnGrid, BidirectionOnGrid, 𝐧ₚᵤ!
using BitmapMapsExtras: indices_scattered, indices_on_grid
using BitmapMapsExtras: radial_distance_glyph, pack_glyphs!
using BitmapMapsExtras: indices_on_grid, default_ij_functor, 𝐊ᵤ!, 𝐊!

!@isdefined(hash_image) && include("common.jl")
@testset "Plot tangent basis" begin
    vhash = ["b942849a8d52f1254eda7dea6a5f104f3c62c137", "ccecdb2fe0c8ee28ef03fb0a12d1aa7d57d9d221"]
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
    vhash = ["ae52041380af7ab9d46e34ab9c80a1a8fba62162", "261ec4726d6e7da5477e249d4e81e72dc6625044"]
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
    vhash = ["fdd6ff2992be35aa512590ee846d42c60ebf082c", "146e9cb3f2d5714453450ef16814488ad83a2676"]
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
    vhash = ["4ca534a1680a2170c65d9ccddb42cad1571d8d53", "ca41ecdb038e73476e8b96e813457fec23cf6ba4"]
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
    vhash = ["d23f090b6971f71c1d5a4aa1628f0f773d34b8cf", "da8d0f6f8b045f11914453c10c9dfc6313b4bbed"]
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
    vhash = ["b5084a8c3e562dd5988107ec48b13467faf5391e", "69c5d9e64c56f3b814c119a33e4704b4106cd552", "f4c905d816dcbcfe22ab71dedb803f428dc11708", "86021239fe92bab1a50a16f4c26ada4a62593555", "c5cbf21a2d4431517a44cb424e0541169f927361", "d47ee5cd2575f812898c4346ce5a48dc7e703b66", "a322d3d1222a686ca6352c94b617e60cf4cc5de3"]
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
    vhash = ["dbd797102d14064ccadf4237c5f9df0878974862", "a712a45d997e776bdc6c75200814efd67247b008", "529acd7e863bff4bfb5a56bfa81dcbbd5b18c9fd", "e5019237ebaf0604b025dc66c9adcf003a93da6c", "bda4e88e3ae13f874f0770b0e04edebd5a5aebb2", "34049a2d4851b289cda6440dfeafbf7b5be6dd1c", "5eedea4b2c100491764c0482cd6411d1e5d248ad", "4b88a9f549f9f622859a876c2a246007f6c6e62e"]
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
    vhash = ["4b6a3165f80caf20201c95108363fae64bcc7035", "40e823377ec0ea9c1637689361533a6d610606b8", "dc18466c31c79fcc5250c5bd9ac002b1178ecc3b", "9b66afb0581ba9b82692563a83e5434316e5fa13", "f73319d6b482fbc7ade899eee3f05f0369be5d50", "71c5ddea9ada7190b4f20776527c1b127f685f23", "fc46dd39cb0d61530e966be2165e943b614b4df7", "3abfc61577d52d7b08fdd56aeecc348ae643c3ff"]
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

vhash = String[]
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
=#