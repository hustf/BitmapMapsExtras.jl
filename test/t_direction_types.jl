using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: DirectionOnGrid, DirectionInDomain, DirectionAtXY, Domain
using BitmapMapsExtras: AbstractXYFunctor
using BitmapMapsExtras: 𝐧ₚ!, MVector

!@isdefined(hashstr) && include("common.jl")

#################
# Low-level tests
#################

@testset "DirectionOnGrid" begin
    z = z_cylinder(0.1)[38:end, 106:end]
    dog = DirectionOnGrid(𝐧ₚ!, z)
    pt = CartesianIndex(5,5)
    dog(pt.I...)
    @test dog.lastvalue[1] < 0 # Normal vector points left
    @test dog.lastvalue[2] > 0 # Normal vector points up ("y is up")
    @test dog.lastvalue ≈ [-0.09903836842989576, 0.9865965311601863]
    pt = CartesianIndex(5,6)
    dog(pt.I...)
    @test dog.lastvalue ≈ [-0.09901472874636574, 0.9863907566538036]
end

@testset "DirectionInDomain" begin
    z = z_cylinder(0.1)[38:end, 106:end]
    did = DirectionInDomain(𝐧ₚ!, z)
    @test DirectionInDomain(𝐧ₚ!, z)(5.0, 5.0) ≈ [-0.09903836842989576, 0.9865965311601863]
    v = MVector{2, Float64}([0, 0])
    v .= did(5.0, 5.0)
    @test v ≈ [-0.09903836842989576, 0.9865965311601863]
    v .= did(6.0, 5.0)
    @test v ≈ [-0.09901472874636574, 0.9863907566538036]
end

@testset "DirectionAtXY" begin
    z = z_cylinder(0.1)[38:end, 106:end]
    daxy = DirectionAtXY(𝐧ₚ!, z)
    @test daxy isa AbstractXYFunctor
    v = daxy(3.0, 3.0)
    @test v ≈ [0.09110987888453688, -0.9080563695155759]
    # Outside domain: zero everywhere, but no error thrown.
    # This is because we find it hard to restrict DiffEq solvers
    # from probing solutions outside the domain.
    @test daxy(2.99, 3.0) ≈ [0, 0]
    v = daxy(892.0, 960.0)
    @test v ≈ [-0.0816455901667088, 0.8137317920512829]
    @test daxy(892.1, 960.0) ≈ [0, 0]
    @test daxy.d == Domain(3.0, 3.0, 892.0, 960.0)
    # This point up and a little to the left
    𝐧ₚx, 𝐧ₚy = daxy(962 / 2, 894 - 100)
    @test abs(𝐧ₚx) < abs(𝐧ₚy)
    @test 𝐧ₚx < 0 # Normal vector points left
    @test 𝐧ₚy > 0 # Normal vector points up ("y is up")
end

# DEV
# daxy = DirectionAtXY(𝐧ₚ!, z_cylinder(0.1)[38:end, 106:end])
# 90.292 ns (0 allocations: 0 bytes)
# @btime daxy(892.0, 960.0)