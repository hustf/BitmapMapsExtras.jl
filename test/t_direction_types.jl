using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: Vec2OnGrid, Vec2InDomain, Vec2AtXY, Domain
using BitmapMapsExtras: AbstractXYFunctor
using BitmapMapsExtras: 𝐧ₚ!, MVector

#################
# Low-level tests
#################

@testset "Vec2OnGrid" begin
    z = z_cylinder(0.1)[38:end, 106:end]
    vog = Vec2OnGrid(𝐧ₚ!, z)
    pt = CartesianIndex(5,5)
    vog(pt.I...)
    @test vog.lastvalue[1] < 0 # Normal vector points left
    @test vog.lastvalue[2] > 0 # Normal vector points up ("y is up")
    @test vog.lastvalue ≈ [-0.09903836842989576, 0.9865965311601863]
    pt = CartesianIndex(5,6)
    vog(pt.I...)
    @test vog.lastvalue ≈ [-0.09901472874636574, 0.9863907566538036]
end

@testset "Vec2InDomain" begin
    z = z_cylinder(0.1)[38:end, 106:end]
    did = Vec2InDomain(𝐧ₚ!, z)
    @test Vec2InDomain(𝐧ₚ!, z)(5.0, 5.0) ≈ [-0.09903836842989576, 0.9865965311601863]
    v = MVector{2, Float64}([0, 0])
    v .= did(5.0, 5.0)
    @test v ≈ [-0.09903836842989576, 0.9865965311601863]
    v .= did(6.0, 5.0)
    @test v ≈ [-0.09901472874636574, 0.9863907566538036]
end

@testset "Vec2AtXY" begin
    z = z_cylinder(0.1)[38:end, 106:end]
    vaxy = Vec2AtXY(𝐧ₚ!, z)
    @test vaxy isa AbstractXYFunctor
    v = vaxy(3.0, 3.0)
    @test v ≈ [0.09110987888453688, -0.9080563695155759]
    @test size(vaxy) == size(z)
    # Outside domain: zero everywhere, but no error thrown.
    # This is because we find it hard to restrict DiffEq solvers
    # from probing solutions outside the domain.
    @test vaxy(2.99, 3.0) ≈ [0, 0]
    v = vaxy(892.0, 960.0)
    @test v ≈ [-0.0816455901667088, 0.8137317920512829]
    @test vaxy(892.1, 960.0) ≈ [0, 0]
    @test vaxy.d == Domain(3.0, 3.0, 892.0, 960.0)
    # This point up and a little to the left
    𝐧ₚx, 𝐧ₚy = vaxy(962 / 2, 894 - 100)
    @test abs(𝐧ₚx) < abs(𝐧ₚy)
    @test 𝐧ₚx < 0 # Normal vector points left
    @test 𝐧ₚy > 0 # Normal vector points up ("y is up")
end

