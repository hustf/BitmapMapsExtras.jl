using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras.TestMatrices: I0
using BitmapMapsExtras: tangent_unit_2d_vector, dz_over_dy
using BitmapMapsExtras: 𝐧ₚ!

Ω = CartesianIndices((-2:2, -2:2))
M = z_cylinder(π/6)[I0 + CartesianIndex(100, 100) .+ Ω] 
@test tangent_basis(M) ≈ [0.9900246401335998 0.03340305238927263 0.13687749274229064; 0.0 0.971490433781295 -0.23707875710706655; -0.14089432894313494 0.23471381118824466 0.9617994670975614]
@test all(tangent_unit_2d_vector(dz_over_dy, M) .≈ [0.9709379516452763, 0.239331347831568])
@test 𝐧ₚ!([0.0, 0.0], M) ≈ [0.13687749274229064, -0.23707875710706658]


# DEV 
#=
using BenchmarkTools
# 163.325 ns (5 allocations: 304 bytes)
@btime tangent_basis(M); 

using BitmapMapsExtras: MVector, MMatrix, tangent_basis!
# Pre-allocate mutable fixed-sized 2d vector
v = MVector{2,Float64}(Array{Float64, 1}(undef, 2))
# Pre-allocate mutable statically-sized array
P = MMatrix{3,3,Float64}(Array{Float64, 2}(undef, (3, 3)))
#  117.165 ns (0 allocations: 0 bytes)
@btime tangent_basis!(P, v, M)


@edit tangent_basis(M) 
=#