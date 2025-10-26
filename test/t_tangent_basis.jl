using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras.TestMatrices: I0
using BitmapMapsExtras: tangent_unit_2d_vector, dz_over_dy, dz_over_dx
using BitmapMapsExtras: 𝐧ₚ!
using BitmapMapsExtras: MVector, MMatrix, tangent_basis!, tangent_unit_2d_vector!

@testset "tangent_basis 1" begin
    Ω = CartesianIndices((-2:2, -2:2))
    M = z_cylinder(π/6)[I0 + CartesianIndex(100, 100) .+ Ω]
    @test tangent_basis(M) ≈ [0.9900246401335998 0.03340305238927263 0.13687749274229064; 0.0 0.971490433781295 -0.23707875710706655; -0.14089432894313494 0.23471381118824466 0.9617994670975614]
    @test all(tangent_unit_2d_vector(dz_over_dy, M) .≈ [0.9709379516452763, 0.239331347831568])
    @test 𝐧ₚ!([0.0, 0.0], M) ≈ [0.13687749274229064, -0.23707875710706658]
end

@testset "tangent_unit_2d_vector" begin
    Ω = CartesianIndices((-2:2, -2:2))
    pt =  CartesianIndex(100, 100)
    z = z_cylinder(0.01)
    win = view(z, Ω .+ pt)
    v = MVector{2,Float64}(Array{Float64, 1}(undef, 2))
    P = MMatrix{3,3,Float64}(Array{Float64, 2}(undef, (3, 3)))
    @test tangent_unit_2d_vector!(v, dz_over_dx, win) ≈ [0.9999049050052322, 0.013790610808710548]
    @test tangent_unit_2d_vector!(v, dz_over_dy, win) ≈ [0.5870136875366161, -0.80957700723567]
end
#=
# DEV 
using BenchmarkTools 
using BitmapMapsExtras: ⋅, KERN1´, KERN2´

z = z_cylinder(0.01)
Ω = CartesianIndices((-2:2, -2:2))
pt =  CartesianIndex(100, 100)
win = view(z, Ω .+ pt)


# Pre-allocate mutable fixed-sized 2d vector
v = MVector{2,Float64}(Array{Float64, 1}(undef, 2))
# Pre-allocate mutable statically-sized array
P = MMatrix{3,3,Float64}(Array{Float64, 2}(undef, (3, 3)))
# 58.596 ns (0 allocations: 0 bytes)
# 80.228 ns (0 allocations: 0 bytes)
# 59.410 ns (0 allocations: 0 bytes)
# 
@btime tangent_unit_2d_vector!(v, dz_over_dx, win) 
@test v ≈ [0.9999049050052322, 0.013790610808710548]
# 68.916 ns (0 allocations: 0 bytes)
# 97.782 ns (0 allocations: 0 bytes)
@btime tangent_unit_2d_vector!(v, dz_over_dy, win) 
@test v ≈ [0.5870136875366161, -0.80957700723567]
dz_over_dx(win)
@code_warntype dz_over_dy(win)

M = z_cylinder(π/6)[I0 + CartesianIndex(100, 100) .+ Ω] 
#  117.165 ns (0 allocations: 0 bytes)
#  140.909 ns (0 allocations: 0 bytes)
# 123.432 ns (0 allocations: 0 bytes)
@btime tangent_basis!(P, v, M)
@code_typed tangent_basis!(P, v, M)

@btime dz_over_dx(M)

=#