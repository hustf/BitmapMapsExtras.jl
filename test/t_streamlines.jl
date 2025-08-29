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

include("common.jl")


#@testset "Streamlines normal projected" begin
    r = TestMatrices.r
    # Non-normalized projection of unit normal vector #  0.217326 seconds (540.22 k allocations: 650.989 MiB)
    @time img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚ!; z = z_cylinder(π / 6))

    @test hashstr(img) == "ac9c072c6646d61ddc9499a45b9d241d51420dfb"
    # Normalized gradient 0.177733 seconds (554.15 k allocations: 651.579 MiB)
    @time img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!; z = z_cylinder(π / 6))
    @test hashstr(img) == "92ae93faa48b0c172498f00d9a69ca7ba5508ff6"
    # Streamlines for all these geometries end up far into a flat area.
    # This is easily fixed by setting keyword argument dtmax 
    img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!; z = z_sphere())
    @test hashstr(img) == "351d579db467f639395b8e65ef999097b12d89ef"
    img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!; z = z_paraboloid())
    @test hashstr(img) == "c600f3d89022799e710add7cb3c571c0d30e699c"
    img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!; z = z_paraboloid(a = 0.5r, b = -r))
    @test hashstr(img) == "0b8e5e0d724f87c755ec80ed1f55f909ff9ca0e6"
    img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!, z = z_cos())
    @test hashstr(img) == "1576d698f3c5d45a8ce7942a4b30017d33491519"
    # In sharp valley bottoms, there's also some oscillation which could be fixed with dtmax 
    img = grid_fcall_with_background(plot_streamlines!, 𝐧ₚᵤ!; z = z_ridge_peak_valleys())
    @test hashstr(img) == "5c053bc7112fb70234a3a5338720fa8580a8c6b7"
    args = ((𝐧ₚᵤ!,), (;dtmax = 5))
    img = grid_fcall_with_background(plot_streamlines!, args; z = z_ridge_peak_valleys())
    @test hashstr(img) == "c4b503d6e6d4e8da4360a839d10668033e6ed795"
    # Short streamlines down from grid points
    args = ((𝐧ₚᵤ!,), (;dtmax = 5, tstop = 200,  rgb = PALETTE_GRGB[3] ))
    img = grid_fcall_with_background(plot_streamlines!, args; z = z_ridge_peak_valleys())
    @test hashstr(img) == "d1668e1b6bfac45070aefb445509c1384dd1903c"
    # Short streamlines up from grid points
    args = ((𝐧ₚᵤ!,), (;dtmax = 5, tstop = -200,  rgb = PALETTE_GRGB[4] ))
    img = grid_fcall_with_background(plot_streamlines!, args; z = z_ridge_peak_valleys())
    @test hashstr(img) == "94e16ecc8570f66a57fe68e17a263bda89e5548a"
    # Short streamlines in both directions from grid points
    args = [((𝐧ₚᵤ!,), (;dtmax = 5, tstop = 200,  rgb = PALETTE_GRGB[3] )),
            ((𝐧ₚᵤ!,), (;dtmax = 5, tstop = -200, rgb = PALETTE_GRGB[4] ))]
    img = grid_fcall_with_background(plot_streamlines!, args; z = z_ridge_peak_valleys())
    @test hashstr(img) == "6eed9576a6b001c2210114cbd4ec91f457bdb6c2"
    # Many streamlines in both directions from grid points
    args = [((𝐧ₚᵤ!,), (;dtmax = 5, tstop = 200,  rgb = PALETTE_GRGB[3], strength = 0.06f0)),
            ((𝐧ₚᵤ!,), (;dtmax = 5, tstop = -200, rgb = PALETTE_GRGB[4], strength = 0.06f0 ))]
    img = grid_fcall_with_background(plot_streamlines!, args; z = z_ridge_peak_valleys(), Δ = 25)
    @test hashstr(img) == "639c9ca44f6a692fd1b4340067845c2325a1f7ce"
#end

#@testset "Streamlines without grid" begin
    args = ((𝐧ₚᵤ!,), (;rgb = PALETTE_GRGB[4], r = 5f0))
    pts = [CartesianIndex(300, 480),
        CartesianIndex(300, 501) ]
    img = fcall_with_background(plot_streamlines!, args, z_cos(), pts)
    mark_at!(img, pts, 15, "in_circle")
    args = ((𝐧ₚᵤ!,), (;rgb = PALETTE_GRGB[3], r = 5f0, tstop = -200))
    fcall(plot_streamlines!, args, img, z_cos(), pts)
    @test hashstr(img) == "ee83342f016540abfcf77b6b69964bd59b4c039f"
    #
    seed_dens = clamp.(z_ridge_peak_valleys(), 0.0, 1.0)
    background(seed_dens)
    n = 1000
    pts = sample(MersenneTwister(123), CartesianIndices(seed_dens), Weights(vec(seed_dens)), n)
    args = ((𝐧ₚᵤ!,), (;rgb = PALETTE_GRGB[3], r = 1f0, dtmax = 1))
    img = fcall_with_background(plot_streamlines!, args, z_ridge_peak_valleys(), pts)
    mark_at!(img, pts, 1, "in_circle")
    @test hashstr(img) == "746bb562bd49022bcfd93c728c35bf77e82c5762"
#end