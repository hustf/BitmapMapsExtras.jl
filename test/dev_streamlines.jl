using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras.TestMatrices: I0
using BitmapMapsExtras: plot_streamlines!, PALETTE_GRGB
using BitmapMapsExtras: SelectedVec2AtXY, 𝐊!, 𝐧ₚ!, Vec2AtXY
using BitmapMapsExtras: Stroke, plot_streamlines, allocations_curvature
using BitmapMapsExtras: is_close_to_perpendicular, is_pointing_roughly_perpendicular
using BenchmarkTools

!@isdefined(hash_image) && include("common.jl")

#r = TestMatrices.r
pts = [CartesianIndex(300, 450),
      CartesianIndex(300, 475),
      CartesianIndex(300, 500)]
z = z_ellipsoid(; tilt = -π / 5 );
img = background(z)
stroke = Stroke(r = 1f0)
saxy = SelectedVec2AtXY(𝐊!, z, true, true);
vaxy = Vec2AtXY(𝐧ₚ!, z);

plot_streamlines!(img, saxy, pts; stroke, dtmax = 1)
plot_streamlines!(img, vaxy, pts; stroke = Stroke(color = PALETTE_GRGB[2]), dtmax = 1)
saxy = SelectedVec2AtXY(𝐊!, z, true, false);
plot_streamlines!(img, saxy, pts; stroke, dtmax = 1)
saxy = SelectedVec2AtXY(𝐊!, z, false, true);
plot_streamlines!(img, saxy, pts; stroke, dtmax = 1)
saxy = SelectedVec2AtXY(𝐊!, z, false, false);
plot_streamlines!(img, saxy, pts; stroke, dtmax = 1)

# Optimize is_close_to_perpendicular / is_pointing_roughly_perpendicular

Ri, Ω, v, P, K, vα, vκ, vβ, lpc = allocations_curvature(TestMatrices.R)

K .= [1 0; 0.1 1]
# 63.367 ns (2 allocations: 64 bytes)
@btime is_close_to_perpendicular(K[:, 1], K[:, 2])
K .= [1 0; 0.1 1] 
# 94.737 ns (2 allocations: 64 bytes)
@btime is_pointing_roughly_perpendicular(saxy, K[:, 1], K[:, 2])

