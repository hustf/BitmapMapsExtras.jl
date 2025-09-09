using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: principal_curvature_components, allocations_curvature
using BitmapMapsExtras: principal_curvature_components!
Ω = CartesianIndices((-2:2, -2:2))
M = z_cylinder(0)[TestMatrices.I0 + CartesianIndex(100, 100) .+ Ω] 
@test true

CM = principal_curvature_components(z_cylinder(0), TestMatrices.I0)
@test isapprox(CM, [-9.372522023933059e-5 0; 0 -0.0020051318513509423])
@test isapprox(CM[1,1], 0, atol = 1e-4)
@test isapprox(CM[2,2], -1 / TestMatrices.r, atol = 1e-5)

#= DEV
z = z_cylinder(0)
pt = TestMatrices.I0
Ri, Ω, v, P, K, vα, vκ, vβ, lpc = allocations_curvature(CartesianIndices(z))
win = view(z, Ω .+ pt)
vϕ = BitmapMapsExtras.VΦ
# 1.050 μs (6 allocations: 304 bytes)
@btime principal_curvature_components!($K, $vα, $vβ, $vκ, $P, $win, $vϕ, $lpc)
=#