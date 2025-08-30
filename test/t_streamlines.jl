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
    # Non-normalized projection of unit normal vector #  0.217326 seconds (540.22 k allocations: 650.989 MiB)
    @time img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚ!; z = z_cylinder(π / 6))
    @test hashstr(img) == "e862b14e9fd8f77356acd592bf013a1f44279cf9" 
    # Normalized gradient 0.177733 seconds (554.15 k allocations: 651.579 MiB)
    @time img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!; z = z_cylinder(π / 6))
    @test hashstr(img) == "058b424d6d35086c65056f3652c562010d062893"
    # Streamlines for all these geometries end up far into a flat area.
    # This is easily fixed by setting keyword argument dtmax 
    img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!; z = z_sphere())
    @test hashstr(img) == "0beedac7bb2ceb4095ee56801aeb06b86b03f13c"
    img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!; z = z_paraboloid())
    @test hashstr(img) == "695aef015afb1a5baf9581c55b213da4b03ac5e8"
    img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!; z = z_paraboloid(a = 0.5r, b = -r))
    @test hashstr(img) == "78c53562d7a1c4a3a8f73b2ffbfd6f2c795c14bf"
    img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!, z = z_cos())
    @test hashstr(img) == "cc77551c537aabb3bdf9725a2068990e1581c81a"
    # In sharp valley bottoms, there's also some oscillation which could be fixed with dtmax 
    img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!; z = z_ridge_peak_valleys())
    @test hashstr(img) == "9263d03b80a2203092f3631ca7bfabf34315eefc"
    args = ((𝐧ₚᵤ!,), (;dtmax = 5))
    img = grid_fcall_with_background(plot_streamlines!, args; z = z_ridge_peak_valleys())
    @test hashstr(img) == "c30a8feb31672bec286b554745f9ac24c99be261" 
    # Short streamlines down from grid points
    args = ((𝐧ₚᵤ!,), (;dtmax = 5, tstop = 200,  rgb = PALETTE_GRGB[3] ))
    img = grid_fcall_with_background(plot_streamlines!, args; z = z_ridge_peak_valleys())
    @test hashstr(img) == "b9517888bd2c859b7b99b2028763850e2ef657b4"
    # Short streamlines up from grid points
    args = ((𝐧ₚᵤ!,), (;dtmax = 5, tstop = -200,  rgb = PALETTE_GRGB[4] ))
    img = grid_fcall_with_background(plot_streamlines!, args; z = z_ridge_peak_valleys())
    @test hashstr(img) == "70fbe64f770decc9672915afef2342cbd6f0dfb2"
    # Short streamlines in both directions from grid points
    args = [((𝐧ₚᵤ!,), (;dtmax = 5, tstop = 200,  rgb = PALETTE_GRGB[3] )),
            ((𝐧ₚᵤ!,), (;dtmax = 5, tstop = -200, rgb = PALETTE_GRGB[4] ))]
    img = grid_fcall_with_background(plot_streamlines!, args; z = z_ridge_peak_valleys())
    @test hashstr(img) == "4a7a47e04fb21d989d53426160cac2a736f567b0"
    # Many streamlines in both directions from grid points
    args = [((𝐧ₚᵤ!,), (;dtmax = 5, tstop = 200,  rgb = PALETTE_GRGB[3], strength = 0.06f0)),
            ((𝐧ₚᵤ!,), (;dtmax = 5, tstop = -200, rgb = PALETTE_GRGB[4], strength = 0.06f0 ))]
    img = grid_fcall_with_background(plot_streamlines!, args; z = z_ridge_peak_valleys(), Δ = 25)
    @test hashstr(img) == "062909da6970fdc52b42846ade1259577d45a1ec"
end

@testset "Streamlines without grid" begin
    args = ((𝐧ₚᵤ!,), (;rgb = PALETTE_GRGB[4], r = 5f0))
    pts = [CartesianIndex(300, 480),
        CartesianIndex(300, 501) ]
    img = fcall_with_background(plot_streamlines!, args, z_cos(), pts)
    mark_at!(img, pts, 15, "in_circle")
    args = ((𝐧ₚᵤ!,), (;rgb = PALETTE_GRGB[3], r = 5f0, tstop = -200))
    fcall(plot_streamlines!, args, img, z_cos(), pts)
    @test hashstr(img) == "63a051719dce7816fe6ba1fc89034fefe0fe34ed"
    #
    seed_dens = clamp.(z_ridge_peak_valleys(), 0.0, 1.0)
    background(seed_dens)
    n = 1000
    pts = sample(MersenneTwister(123), CartesianIndices(seed_dens), Weights(vec(seed_dens)), n)
    args = ((𝐧ₚᵤ!,), (;rgb = PALETTE_GRGB[3], r = 1f0, dtmax = 1))
    img = fcall_with_background(plot_streamlines!, args, z_ridge_peak_valleys(), pts)
    mark_at!(img, pts, 1, "in_circle")
    @test hashstr(img) == "b3d3e0b5afa59fcec35d09fe0bcf7079c4c9ec1f"
end