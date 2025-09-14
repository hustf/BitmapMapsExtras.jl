using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: principal_curvature_components, allocations_curvature
using BitmapMapsExtras: principal_curvature_components!, VΦ, sample_at_float
using BitmapMapsExtras: sampled_second_derivative_in_direction
Ω = CartesianIndices((-2:2, -2:2))
M = z_cylinder(0)[TestMatrices.I0 + CartesianIndex(100, 100) .+ Ω] 
@test true

CM = principal_curvature_components(z_cylinder(0), TestMatrices.I0)
@test isapprox(CM, [-9.372522023933059e-5 0; 0 -0.0020051318513509423])
@test isapprox(CM[1,1], 0, atol = 1e-4)
@test isapprox(CM[2,2], -1 / TestMatrices.r, atol = 1e-5)

# DEV
#=
using BenchmarkTools

# 12.813 ns (0 allocations: 0 bytes)
# 15.832 ns (0 allocations: 0 bytes)
# 20.963 ns (0 allocations: 0 bytes)
# 9.810 ns (0 allocations: 0 bytes)
@btime sample_at_float($M, 1.5, 1.5)

sample_at_float(M, 1.5, 1.5) ≈ 489.1814488141001
# 50.253 ns (0 allocations: 0 bytes)
# 45.556 ns (0 allocations: 0 bytes)
@btime sampled_second_derivative_in_direction($M, 0, 1.0, 0.0)

z = z_cylinder(0)
pt = TestMatrices.I0
Ri, Ω, v, P, K, vα, vκ, vβ, lpc = allocations_curvature(CartesianIndices(z))
win = view(z, Ω .+ pt)
vϕ = BitmapMapsExtras.VΦ
#  1.030 μs (5 allocations: 240 bytes)
# 1.050 μs (6 allocations: 304 bytes)
# 1.020 μs (5 allocations: 240 bytes)
# 707.857 ns (1 allocation: 48 bytes)
# 627.059 ns (0 allocations: 0 bytes)
@btime principal_curvature_components!($K, $vα, $vβ, $vκ, $P, $win, $vϕ, $lpc)


function pro()
    for i = 1:1000000
        principal_curvature_components!(K, vα, vβ, vκ, P, win, vϕ, lpc)
    end
end

# 1.092 s (5000000 allocations: 228.88 MiB)
#1.112 s (5000000 allocations: 228.88 MiB)
# 729.755 ms (1000000 allocations: 45.78 MiB)
# 669.853 ms (0 allocations: 0 bytes)
@btime pro()


# 705.714 ns (1 allocation: 48 bytes)
# 631.953 ns (0 allocations: 0 bytes)
@btime principal_curvature_components!($K, $vα, $vβ, $vκ, $P, $win, $vϕ, $lpc);

@profview_allocs pro()


@code_warntype principal_curvature_components!(K, vα, vβ, vκ, P, win, vϕ, lpc)



# Unroll principal_curvature_components!
using BitmapMapsExtras: tangent_basis!, angle_tangent_to_xy!, sample_curvature_in_directions!,
    principal_curvature_and_direction
# 0
@btime tangent_basis!(P, vβ, M)
#0
@btime angle_tangent_to_xy!(vα, vϕ, P)
# 700.709 ns (4 allocations: 192 bytes)
@btime sample_curvature_in_directions!($vκ, $M, $P, $vα);
# Dive deeper
function pri()
    for i = 1:10000000
        sample_curvature_in_directions!(vκ, M, P, vα)
    end
end

@btime pri()
@profview pri()

α = vα[1]
αl = α - π / 2 
using BitmapMapsExtras: recenter, RotMatrix, center, WarpedView, Flat, AffineMap, W
tfm = recenter(RotMatrix(αl), center(M))
# 3 allocs, 144 bytes
@btime recenter(RotMatrix(αl), center(M))
# 2 allocs, 112 bytes
@btime recenter(RotMatrix(αl), (3.0, 3.0))
tfm1 = recenter(RotMatrix(αl), (3.0, 3.0))
# can we use LinearMap directly?
@btime LinearMap(RotMatrix(αl))


Mr = WarpedView(M, tfm; fillvalue = Flat())
Mrw = WarpedView(M, tfm; fillvalue = Flat())[W]

sum(KERN´´ .* Mr[W], init = 0.0) / (1 + (-(cos(α) * P[7] + sin(α) * P[8]) / P[9])^2)^(3 / 2)

=#