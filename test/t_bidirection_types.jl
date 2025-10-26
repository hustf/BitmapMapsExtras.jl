using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: BidirectionOnGrid, BidirectionInDomain, BidirectionAtXY, Domain
using BitmapMapsExtras: AbstractXYFunctor, SelectedVec2AtXY
using BitmapMapsExtras: 𝐊!, norm, reset!


#################
# Low-level tests
#################

@testset "Constant curvature" begin
    radius = TestMatrices.r
    z = z_cylinder(0.0)[300:500, 100:300];
    @testset "BidirectionOnGrid" begin
        bog = BidirectionOnGrid(𝐊!, z)
        pt = (100, 100)
        @test isapprox(bog(pt...), [0 0; 0   -1/radius],  atol = 2e-4)
        @inferred bog(pt...)
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

@testset "SelectedVec2AtXY" begin
    # Smoothly varying curvature, largest at centre lines
    z = z_paraboloid()
    saxy = SelectedVec2AtXY(𝐊!, z, true, false)
    @test saxy isa AbstractXYFunctor
    @test saxy(3.0, 3.0) ≈ [0.0003794227888056931, 0.0008893924914691132]
    vu = [norm(saxy(100.0, y)[:, 1]) for y in 50:0.2:53]
    @test sort(vu) == vu
    @test length(unique(vu)) == length(vu)
    vu = [norm(saxy(x, 100.0)[:, 1]) for x in 50:0.2:53]
    @test sort(vu, rev = true) == vu
    @test length(unique(vu)) == length(vu)
    @test size(saxy) == size(z)
    # Outside domain: zero everywhere, but no error thrown.
    # This is because we find it hard to restrict DiffEq solvers
    # from probing solutions outside the domain.
    @test saxy(2.99, 3.0) ≈ [0, 0]
    @test saxy(997.0, 997.0) ≈ [0.0003794227888056931, 0.0008893924914691132]
    @test saxy(997.01, 997.0) ≈ [0, 0]
    @test saxy(997.00, 997.01) ≈ [0, 0]
    @test saxy(892.0, 960.0) ≈ [0.0002916029305382096, 0.0009477262727801998]
    saxy.flip[] = ! saxy.flip[]
    @test saxy(892.0, 960.0) ≈ [-0.0002916029305382096, -0.0009477262727801998]
    saxy.baxy.primary
    saxy.flip
    @test saxy.flip[]
    reset!(saxy)
    @test ! saxy.flip[]
end
