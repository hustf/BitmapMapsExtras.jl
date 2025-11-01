using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: plot_streamlines!, 𝐧ₚ!, 𝐧ₚᵤ!, PALETTE_GRGB
using BitmapMapsExtras: indices_on_grid, Stroke, Vec2AtXY
using BitmapMapsExtras: mark_at!, display_if_vscode
#using BitmapMaps: divergence_of_gradients
import StatsBase
using StatsBase: Weights, sample
import Random
using Random: MersenneTwister

!@isdefined(is_hash_stored) && include("common.jl")


@testset "Streamlines, 𝐧ₚᵤ! versus 𝐧ₚ!" begin
    vhash = ["53497a3460e8ab75aa87377c1c1cc00c8d1f7716", "67e182cd51d498f03497830c8bc5c37179b9deb0"]
    COUNT[] = 0
    vaxy = Vec2AtXY(𝐧ₚ!, z_cylinder(π / 6));
    pts = indices_on_grid(vaxy)
    img = background(vaxy)
    # @time 0.832305 seconds (1.91 M allocations: 686.920 MiB, 48.47% gc time)
    # Note that the allocations occur while sampling the solution
    plot_streamlines!(img, vaxy, pts; stroke = Stroke(color = PALETTE_GRGB[2]), dtmax = 1)
    @test is_hash_stored(img, vhash)
    # descent_unit!
    vaxy = Vec2AtXY(𝐧ₚᵤ!, z_cylinder(π / 6));
    img = background(vaxy)
    # @time 0.528286 seconds (1.57 M allocations: 675.361 MiB, 32.89% gc time)
    # @time  0.654856 seconds (1.56 M allocations: 675.360 MiB, 44.05% gc time)
    # This is faster, although that doesn't mean much
    plot_streamlines!(img, vaxy, pts; stroke = Stroke(color = PALETTE_GRGB[2]), dtmax = 1)
    @test is_hash_stored(img, vhash)
end

@testset "Solution keyword example" begin
    vhash = ["095dd2d6e0738c196d7799657a6cb2fea1e4cc92", "5d42e47944fe189ed3f4bc25231c5c56ae85587c"]
    COUNT[] = 0
    # Streamlines ending up far into a flat area.
    vaxy = Vec2AtXY(𝐧ₚ!, z_sphere());
    pts = indices_on_grid(vaxy)
    img = background(vaxy)
    plot_streamlines!(img, vaxy, pts)
    @test is_hash_stored(img, vhash)
    # Setting keyword argument dtmax => no flat area solution points.
    img = background(vaxy)
    plot_streamlines!(img, vaxy, pts; dtmax = 1.0)
    @test is_hash_stored(img, vhash)
end

@testset "Streamlines starting on grid" begin
    vhash = ["4d9bd22adb0a69c7f3b42a60015371595faeade4", "3af97fc101ae79685a71f8c67c6a7159c3a77add", "58fdf26029fc652daaf2fbc687d6a57fff790456", "ae79ba6104dd47623a95642f62352e3318159dab", "119030a66f8c00cf9da18b7c42f31e562fbd533d", "49f844e4e80eb4ce9798cf11ba3b190462fa3417", "555f82f483e25ba6d41bdc8ad709752d57899f6f", "5d42e47944fe189ed3f4bc25231c5c56ae85587c"]
    COUNT[] = 0
    vzf = [z_cos,
        () -> z_cylinder(π / 6),
        () -> z_cylinder_offset(π / 3),
        z_ellipsoid,
        z_exp3,
        z_paraboloid,
        z_ridge_peak_valleys,
        z_sphere]
    for fz in vzf
        vaxy = Vec2AtXY(𝐧ₚ!, fz());
        pts = indices_on_grid(vaxy)
        img = background(vaxy)
        # dtmax = 1.0 is often a good idea.
        plot_streamlines!(img, vaxy, pts; dtmax = 1.0)
        @test is_hash_stored(img, vhash)
    end
end


@testset "Uphill and downhill with different appearance" begin
    vhash = ["0638a99a6775677c5f02c0363f4e16be2b68c48a", "3773fd624e96c630c676c4b47dbda48a04c8faff", "f44bece007e5fe0d20e4f5f332242e09fbed25eb"]
    COUNT[] = 0
    # Short streamlines down from grid points
    vaxy = Vec2AtXY(𝐧ₚ!, z_ridge_peak_valleys())
    stroke = Stroke(color = PALETTE_GRGB[3])
    pts = indices_on_grid(vaxy)
    img = background(vaxy)
    plot_streamlines!(img, vaxy, pts; dtmax = 5, tstop = 200, stroke)
    mark_at!(img, pts, 5, "in_circle")
    @test is_hash_stored(img, vhash)
    # Short streamlines up from grid points (by ode keywords, faster than modifying the function)
    vaxy = Vec2AtXY(𝐧ₚ!, z_ridge_peak_valleys())
    stroke = Stroke(color = PALETTE_GRGB[4])
    pts = indices_on_grid(vaxy)
    img = background(vaxy)
    plot_streamlines!(img, vaxy, pts; dtmax = 5, tstop = -200, stroke)
    mark_at!(img, pts, 5, "in_circle")
    @test is_hash_stored(img, vhash)
    # Many streamlines in both directions from grid points
    # Note we use the normalized vector function here.
    vaxy = Vec2AtXY(𝐧ₚᵤ!, z_ridge_peak_valleys()) 
    pts = indices_on_grid(vaxy, Δ = 25)
    img = background(vaxy)
    plot_streamlines!(img, vaxy, pts; dtmax = 1, tstop = 200, 
        stroke = Stroke(color = PALETTE_GRGB[3], strength = 0.3))
    plot_streamlines!(img, vaxy, pts; dtmax = 1, tstop = -200, 
        stroke = Stroke(color = PALETTE_GRGB[4], strength = 0.3))
    mark_at!(img, pts, 3, "in_circle")
    @test is_hash_stored(img, vhash)
end
@testset "Streamlines without grid" begin
    vhash = ["2829699623bcb078046b54de0951b22e9c47e401", "4a0ec1d800241c3834c48ea4ddc45676287b90a0"]
    COUNT[] = 0
    z = z_ridge_peak_valleys()
    # Originating at points where z > 0
    seed_dens = clamp.(z, 0.0, 1.0);
    display_if_vscode(background(seed_dens))
    n = 1000
    pts = sample(MersenneTwister(123), CartesianIndices(seed_dens), Weights(vec(seed_dens)), n);
    vaxy = Vec2AtXY(𝐧ₚᵤ!, z)
    stroke = Stroke(color = PALETTE_GRGB[3], r = 1, strength = 0.3)
    img = background(z)
    mark_at!(img, pts, 3, "in_circle")
    plot_streamlines!(img, vaxy, pts; dtmax = 1, stroke)
    @test is_hash_stored(img, vhash)
    #
    #=
    # Originating at 'sources' (convex terrain)
    vals1 = divergence_of_gradients(-z_ridge_peak_valleys())
    vals2 = clamp.(vals1, 0.003, 1.0) .- 0.003
    seed_dens = clamp.(vals2 .* 100, 0, 1.0)
    # Drop the extreme area near the centre.
    seed_dens[400:600, 300:700] .= 0.0
    n = 1000
    pts = sample(MersenneTwister(123), CartesianIndices(seed_dens), Weights(vec(seed_dens)), n);
    img = background(z)
    mark_at!(img, pts, 3, "in_circle")
    plot_streamlines!(img, vaxy, pts; dtmax = 1)
    @test is_hash_stored(img, vhash)
    =#
end
