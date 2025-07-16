using Test
using BitmapMapsExtras
using BitmapMapsExtras: signed_distance_within_domain, NegateY, Domain

@testset "NegateY" begin
    # 6 high, 10 wide matrix
    R = CartesianIndices((6, 10))
    flipy = NegateY(R)
    @test flipy(0) == 7
    @test flipy(1) == 6
    @test flipy(2) == 5
    @test flipy(6) == 1
    @test flipy(6.0) == 1
end

@testset "NegateY reversible" begin
    R = CartesianIndices((6, 10))
    flipy = NegateY(R)
    y = 1.0
    i = flipy(y)
    @test flipy(i) == y
    y = 6.0
    i = flipy(y)
    @test flipy(i) == y
end

@testset "Domain" begin
    # 8 high, 10 wide matrix
    R = CartesianIndices((8, 10))
    d = Domain(R)
    @test d == Domain(1.0, 1.0, 10.0, 8.0)
    #
    Ω = CartesianIndices((-2:2, -2:2))
    d = Domain(R, Ω)
    @test d == Domain(3.0, 3.0, 8.0, 6.0)
    @test !d(0.5, 0.3)
    @test !d(0.5, 3)
    @test d(5, 3)
    @test d(5.5, 3.5)
    @test !d(1000, 32)
    @test !d(30, 3200)
    @test !d(3045, 3200)
    #
    @test signed_distance_within_domain(d, 3.0, 4.0) == 0.0
    @test signed_distance_within_domain(d, 8.0, 4.0) == 0.0
    @test signed_distance_within_domain(d, 4.0, 3.0) == 0.0
    @test signed_distance_within_domain(d, 4.0, 6.0) == 0.0
    #
    @test signed_distance_within_domain(d, 3.5, 4.0) == 0.5
    @test signed_distance_within_domain(d, 7.5, 4.0) == 0.5
    @test signed_distance_within_domain(d, 4.0, 3.5) == 0.5
    @test signed_distance_within_domain(d, 4.0, 5.5) == 0.5
    #
    @test signed_distance_within_domain(d, 2.5, 4.0) == -0.5
    @test signed_distance_within_domain(d, 8.5, 4.0) == -0.5
    @test signed_distance_within_domain(d, 4.0, 2.5) == -0.5
    @test signed_distance_within_domain(d, 4.0, 6.5) == -0.5
end

