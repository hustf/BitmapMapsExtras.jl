using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras.TestMatrices: I0
using BitmapMapsExtras: principal_curvature_components, allocations_curvature

Ω = CartesianIndices((-2:2, -2:2))
M = z_cylinder(0)[I0 + CartesianIndex(100, 100) .+ Ω] 
@test true

CM = principal_curvature_components(z_cylinder(0), I0)
@test isapprox(CM, [-9.372522023933059e-5 0; 0 -0.0020051318513509423])
@test isapprox(CM[1,1], 0, atol = 1e-4)
@test isapprox(CM[2,2], -1/r, atol = 1e-5)

# No directions defined
_, _, _, _, _, _, _, _, f_is_within_limits = allocations_curvature(CartesianIndices((100,100)), [])
@test f_is_within_limits(2) == 2
@test f_is_within_limits(CM[:,1]) == CM[:,1]
@test f_is_within_limits(CM[:,2]) == CM[:,2]

# Principal curvature, can't pass two curvatures
_, _, _, _, _, _, _, _, f_is_within_limits  = allocations_curvature(CartesianIndices((100,100)), [1])
@test_throws MethodError f_is_within_limits(CM)
@test f_is_within_limits(CM[:,1]) == true
@test f_is_within_limits(CM[:,2]) == true

# Principal curvatures, can't pass single curvature
_, _, _, _, _, _, _, _, f_is_within_limits  = allocations_curvature(CartesianIndices((100,100)), [1, 2])
@test f_is_within_limits(CM) == true
@test_throws BoundsError f_is_within_limits(CM[:,1])
@test f_is_within_limits(CM) == true
