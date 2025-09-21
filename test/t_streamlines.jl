using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: plot_streamlines!, 𝐧ₚ!, 𝐧ₚᵤ!, PALETTE_GRGB
using BitmapMapsExtras: spray!, apply_color_by_coverage!, RGB
import BitmapMaps
using BitmapMaps: mark_at!, divergence_of_gradients
import StatsBase
using StatsBase: Weights, sample    
import Random
using Random: MersenneTwister

!@isdefined(hashstr) && include("common.jl")


@testset "Streamlines normal projected" begin
    r = TestMatrices.r
    # Non-normalized projection of unit normal vector 
    #  0.217326 seconds (540.22 k allocations: 650.989 MiB)
    #  0.370902 seconds (422.12 k allocations: 653.901 MiB, 48.33% gc time)
    # 0.573795 seconds (646.60 k allocations: 660.699 MiB, 65.42% gc time)
    # 0.369479 seconds (1.27 M allocations: 679.845 MiB, 40.04% gc time)
    @time img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚ!; z = z_cylinder(π / 6))
    @test hashstr(img) == "55ee19054428fd459f6145a604e07891bf9167a6"
    # Normalized gradient 
    # 0.177733 seconds (554.15 k allocations: 651.579 MiB)
    # 0.364872 seconds (465.50 k allocations: 657.474 MiB, 43.03% gc time)
    # 0.350382 seconds (526.51 k allocations: 657.053 MiB, 44.78% gc time)
    # 0.405951 seconds (1.04 M allocations: 672.613 MiB, 42.18% gc time, 4.81% compilation time)
    @time img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!; z = z_cylinder(π / 6))
    @test hashstr(img) == "e7e1bf4cdb1068ebcd26cae723b9103408e2caae"
    # Streamlines for all these geometries end up far into a flat area.
    # This is easily fixed by setting keyword argument dtmax 
    img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!; z = z_sphere())
    @test hashstr(img) == "8ae0a64685db85196fa68fee86c84e376be1ee9f"
    img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!; z = z_paraboloid())
    @test hashstr(img) == "eb48c12bdc4beadb5c57572fa53392511424684b"
    img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!; z = z_paraboloid(a = 0.5r, b = -r))
    @test hashstr(img) == "8a2a9d586303f29c10a8d3b839c2ef16ca0a8cf4"
    img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!, z = z_cos())
    @test hashstr(img) == "e7c9afc7dbdc13dfe1da1489f1ca9c15473e0823"
    # In sharp valley bottoms, there's also some oscillation which could be fixed with dtmax 
    img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!; z = z_ridge_peak_valleys())
    @test hashstr(img) == "a1dd6cba5a52fd25df5179440bc7462d45d2ab3f"
    args = ((𝐧ₚᵤ!, ), (;sol_density = 0.05, dtmax = 5))
    img = grid_fcall_with_background(plot_streamlines!, args; z = z_ridge_peak_valleys())
    @test hashstr(img) == "d1de5a9d82037598d751200f84db3a85195b37cd"
    # Short streamlines down from grid points
    args = ((𝐧ₚᵤ!,), (;dtmax = 5, tstop = 200,  rgb = PALETTE_GRGB[3] ))
    pts = grid_indices(size(img))
    img = fcall_with_background(plot_streamlines!, args, z_ridge_peak_valleys(), pts )
    mark_at!(img, pts, 5, "in_circle")
    @test hashstr(img) == "252c73338cd11459655287a32ac42adf71fe9bd7"
    # Short streamlines up from grid points
    args = ((𝐧ₚᵤ!,), (;dtmax = 5, tstop = -200,  rgb = PALETTE_GRGB[4] ))
    pts = grid_indices(size(img))
    img = fcall_with_background(plot_streamlines!, args, z_ridge_peak_valleys(), pts )
    mark_at!(img, pts, 5, "in_circle")
    @test hashstr(img) == "6b1c48121f1758eb05460cbb1d1b3332c63ae1eb"
    # Short streamlines in both directions from grid points
    args = [((𝐧ₚᵤ!,), (;dtmax = 5, tstop = 100,  rgb = PALETTE_GRGB[3] )),
            ((𝐧ₚᵤ!,), (;dtmax = 5, tstop = -100, rgb = PALETTE_GRGB[4] ))]
    pts = grid_indices(size(img))
    img = fcall_with_background(plot_streamlines!, args, z_ridge_peak_valleys(), pts )
    mark_at!(img, pts, 5, "in_circle")
    @test hashstr(img) == "e964414b7feed44e6d3ea3a464dbe82f8411d83c"
    # Many streamlines in both directions from grid points
    args = [((𝐧ₚᵤ!,), (;dtmax = 5, tstop = 200,  rgb = PALETTE_GRGB[3], strength = 0.06f0)),m n
            ((𝐧ₚᵤ!,), (;dtmax = 5, tstop = -200, rgb = PALETTE_GRGB[4], strength = 0.06f0 ))]
    pts = grid_indices(size(img); Δ = 25)
    img = fcall_with_background(plot_streamlines!, args, z_ridge_peak_valleys(), pts )
    mark_at!(img, pts, 3, "in_circle")
    @test hashstr(img) == "8aeaf33d8279abb03cc93bcbdcd127e24ec887d2"
end

@testset "Streamlines without grid" begin
    args = ((𝐧ₚᵤ!,), (;rgb = PALETTE_GRGB[4], r = 5f0))
    pts = [CartesianIndex(300, 480),
        CartesianIndex(300, 501) ]
    img = fcall_with_background(plot_streamlines!, args, z_cos(), pts)
    mark_at!(img, pts, 15, "in_circle")
    args = ((𝐧ₚᵤ!,), (;rgb = PALETTE_GRGB[3], r = 5f0, tstop = -200))
    fcall(plot_streamlines!, args, img, z_cos(), pts)

    @test hashstr(img) == "909c03d38c6baa3aa054004a52c61a6a1c6ed016"
    #
    seed_dens = clamp.(z_ridge_peak_valleys(), 0.0, 1.0);
    background(seed_dens)
    n = 1000
    pts = sample(MersenneTwister(123), CartesianIndices(seed_dens), Weights(vec(seed_dens)), n);
    args = ((𝐧ₚᵤ!,), (;rgb = PALETTE_GRGB[3], r = 1f0, dtmax = 1))
    img = fcall_with_background(plot_streamlines!, args, z_ridge_peak_valleys(), pts)
    mark_at!(img, pts, 3, "in_circle")
    @test hashstr(img) == "7ea5e7e56304e3c2d249e55ddfe857898acb71a4"
    # Streamlines starting at 'sources', moving towards 'sinks'
    seed_dens = clamp.(divergence_of_gradients(-z_ridge_peak_valleys()), 0.0, 0.03)
    background(seed_dens)
    n = 1000
    pts = sample(MersenneTwister(123), CartesianIndices(seed_dens), Weights(vec(seed_dens)), n);
    args = ((𝐧ₚᵤ!,), (;rgb = PALETTE_GRGB[2], r = 1f0, dtmax = 1))
        args = ((𝐧ₚᵤ!,), (; r = 1f0, dtmax = 1))
    img = fcall_with_background(plot_streamlines!, args, z_ridge_peak_valleys(), pts)
    mark_at!(img, pts, 1, "in_circle")
    @test hashstr(img) == "8a4b9355ae8bcda87b70ff8933fdbb599b93ceeb"
end
# DEV Allocs
#=
using BitmapMapsExtras: DirectionAtXY, vu0_from_pts, make_tspan, callbacks_streamlines, MVector, ODEProblem, rhs!
using BenchmarkTools

f = 𝐧ₚ!
z = z_cylinder(0)
u = (30.4, 50.2)
daxy = DirectionAtXY(f, z)

@inferred daxy(30.4, 50.2)
# 501.042 ns (14 allocations: 640 bytes)
# 413.065 ns (10 allocations: 448 bytes) # {T} types
# 322.707 ns (6 allocations: 192 bytes)  # {F, T} types
# 331.963 ns (6 allocations: 192 bytes)  # {F, T, L} types
# 271.104 ns (0 allocations: 0 bytes)
@btime daxy(30.4, 50.2)

# Dive to the bottom
using BitmapMapsExtras: DirectionOnGrid
dog = DirectionOnGrid(𝐧ₚ!,  z)
@inferred dog(100,200)
# 74.177 ns (0 allocations: 0 bytes)
@btime dog(100,200)

using BitmapMapsExtras: DirectionInDomain
did = DirectionInDomain(𝐧ₚ!,  z)
@inferred did(100.2,200.2)
# 264.181 ns (0 allocations: 0 bytes)
@btime did(100.2,200.2)

@test did(100.0,200.0) == dog(200,100)
did(100.0,200.5)
did(100.0,201.0)


=#