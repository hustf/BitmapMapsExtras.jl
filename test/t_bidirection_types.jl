using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: BidirectionOnGrid, BidirectionInDomain, BidirectionAtXY, Domain
using BitmapMapsExtras: DirectionFunctor, UnidirectionAtXY
using BitmapMapsExtras: direction_at_xy!, direction_on_grid!, direction_in_domain!
using BitmapMapsExtras: 𝐊!, MVector, MMatrix, GrayA, norm, reset!

include("common.jl")
radius = TestMatrices.r
#################
# Low-level tests
#################

@testset "Constant curvature" begin
    z = z_cylinder(0.0)[300:500, 100:300];
    @testset "BidirectionOnGrid" begin
        bog = BidirectionOnGrid(𝐊!, z)
        pt = CartesianIndex(100,100)
        @test isapprox(bog((pt)), [0 0; 0   -1/radius],  atol = 2e-4)
        @inferred bog(pt)
    end
    @testset "BidirectionInDomain" begin 
        bid = BidirectionInDomain(𝐊!, z)
        # Exact on grid
        bid(100.0, 100.0)
        @test isapprox(bid(100.0, 100.0), [0 0; 0 -1/radius],  atol = 2e-4)
        @inferred bid(100.0, 100.0)
    end
end
@testset "BidirectionAtXY" begin
    # Smoothly varying curvature, largest at centre lines
    z = z_paraboloid()
    baxy = BidirectionAtXY(𝐊!, z, true) 
    @test baxy(3.0, 3.0) ≈ [0.0003794227888053568, 0.0008893924914720936]
    vbu = [norm(baxy(100.0, y)[:, 1]) for y in 50:0.2:53]
    @test sort(vbu) == vbu
    @test length(unique(vbu)) == length(vbu)
    vbu = [norm(baxy(x, 100.0)[:, 1]) for x in 50:0.2:53]
    @test sort(vbu, rev = true) == vbu
    @test length(unique(vbu)) == length(vbu)
    # Outside domain: zero everywhere, but no error thrown.
    # This is because we find it hard to restrict DiffEq solvers
    # from probing solutions outside the domain.
    @test baxy(2.99, 3.0) ≈ [0, 0]
    @test baxy(997.0, 997.0) ≈ [0.0003794227888056931, 0.0008893924914691132]
    @test baxy(997.01, 997.0) ≈ [0,  0]
    @test baxy(997.00, 997.01) ≈ [0,  0]
    @test baxy(892.0, 960.0) ≈  [0.0002916029305382096, 0.0009477262727801998]
    baxy.primary[] = ! baxy.primary[]
    @test baxy(892.0, 960.0) ≈ [-0.0012225868641178692, -0.0002940188601622608]
end

@testset "UnidirectionAtXY" begin
    # Smoothly varying curvature, largest at centre lines
    z = z_paraboloid()
    uxy = UnidirectionAtXY(𝐊!, z, true, false) 
    @test uxy isa DirectionFunctor
    @test uxy(3.0, 3.0) ≈ [0.0003794227888056931, 0.0008893924914691132]
    vu = [norm(uxy(100.0, y)[:, 1]) for y in 50:0.2:53]
    @test sort(vu) == vu
    @test length(unique(vu)) == length(vu)
    vu = [norm(uxy(x, 100.0)[:, 1]) for x in 50:0.2:53]
    @test sort(vu, rev = true) == vu
    @test length(unique(vu)) == length(vu)
    # Outside domain: zero everywhere, but no error thrown.
    # This is because we find it hard to restrict DiffEq solvers
    # from probing solutions outside the domain.
    @test uxy(2.99, 3.0) ≈ [0, 0]
    @test uxy(997.0, 997.0) ≈ [0.0003794227888056931, 0.0008893924914691132]
    @test uxy(997.01, 997.0) ≈ [0, 0]
    @test uxy(997.00, 997.01) ≈ [0, 0]
    @test uxy(892.0, 960.0) ≈ [0.0002916029305382096, 0.0009477262727801998]
    uxy.flip[] = ! uxy.flip[]
    @test uxy(892.0, 960.0) ≈ [-0.0002916029305382096, -0.0009477262727801998]
    uxy.baxy.primary
    uxy.flip
    @test uxy.flip[]
    reset!(uxy)
    @test ! uxy.flip[]
end
















# 1.690 μs (28 allocations: 1.84 KiB)
#@btime baxy(892.0, 960.0)


# A good part is in principal_curvature_and_direction. 
# Optimization postponed.
#@profview begin 
#      for y = 1.0:0.5:100
#          baxy(100.0, y)
#      end
#end
