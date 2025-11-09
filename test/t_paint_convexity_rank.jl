using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: paint_convexity_rank!, indices_on_grid, RGB
using BitmapMapsExtras: convexity_rank, BidirectionOnGrid, 𝐊!
using BitmapMapsExtras: is_bidirec_vect_positive, N0f8

!@isdefined(is_hash_stored) && include("common.jl")

@testset "Convexity rank" begin
    r = TestMatrices.r
    flatval = 0.0001
    # 1    Convex - convex           Green
    bog = BidirectionOnGrid(𝐊!, z_paraboloid(; a= 0.6r, b = 0.5r))
    @test convexity_rank(bog(400, 400), flatval) == 1
    # 2     Convex - flat            Bleak green
    bog = BidirectionOnGrid(𝐊!, -z_cylinder(0.1))
    @test convexity_rank(bog(400, 400), flatval) == 2
    # 3     Convex - concave         Blue
    bog = BidirectionOnGrid(𝐊!, z_paraboloid(;a = 0.6r, b = -0.4r))
    @test convexity_rank(bog(400, 400), flatval) == 3
    # 4     Flat - flat              Grey
    bog = BidirectionOnGrid(𝐊!, z_plane())
    @test convexity_rank(bog(400, 400), flatval) == 4
    # 5    Flat - concave            Bleak red
    bog = BidirectionOnGrid(𝐊!, z_cylinder(0.1))
    @test convexity_rank(bog(400, 400), flatval) == 5
    # 6    Concave - concave         Red
    bog = BidirectionOnGrid(𝐊!, z_ellipsoid())
    @test convexity_rank(bog(400, 400), flatval) == 6
end



@testset "Edge case, curvature sign" begin
    i, j = 500, 232
    bog = BidirectionOnGrid(𝐊!, z_ellipsoid())
    K = bog(i,j)
    @test !is_bidirec_vect_positive(K[:, 1])
    @test !is_bidirec_vect_positive(K[:, 2])
end

@testset "Curvature types" begin
    vhash = ["dc56b23e5a88e99204b04937551526bc41dc80fc", "25404f3297b8a41cd074246e42341df9548ac600", "13dae453f55d750db5a2321af05bb2293d65dba4", "7945d11d6a4412d0b649e8bbcbd4acc2fbf83a56", "ab99155725d6cb9a1499792a3851510f67b23426", "8344b3ec23f78ac1b5cb9c455c088ef5cb123bd8"]
    COUNT[] = 0
    # We compromise on this value
    flatval = 0.0003
    pts = TestMatrices.R
    r = TestMatrices.r
    function paint_convexity(z, flatval)
        # A bleaker version of the background lets the convexity rank
        # stand out more
        img = RGB.(background(z))
        img .*= 0.5
        img .+= RGB{N0f8}(0.5, 0.5, 0.5)
        paint_convexity_rank!(img, z, pts; flatval)
    end
    # 1    Convex - convex           Green
    img =  paint_convexity(z_paraboloid(; a= 0.6r, b = 0.5r), flatval)
    is_hash_stored(img, vhash)
    # 2     Convex - flat            Bleak green
    img =  paint_convexity(-z_cylinder(0.1), flatval)
    is_hash_stored(img, vhash)
    # 3     Convex - concave         Blue
    img =  paint_convexity(z_paraboloid(;a= 0.6r, b = -0.4r), flatval)
    is_hash_stored(img, vhash)
    # 4     Flat - flat              Grey
    img =  paint_convexity(z_plane(), flatval)
    is_hash_stored(img, vhash)
    # 5    Flat - concave            Bleak red
    img =  paint_convexity(z_cylinder(0.1), flatval)
    is_hash_stored(img, vhash)
    # 6    Concave - concave         Red
    img =  paint_convexity(z_ellipsoid(), flatval)
    is_hash_stored(img, vhash)
end

