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
    vhash = ["e4afae413dec900f7124795dd323ca189b48513e", "8c8eb2d82316d0274ec6a7652d378933487bdaab"]
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
    vhash = ["91cd4ba192a5e2c8593298ac153913f42b6d5817", "6b876dfa7c901ccae0a164f540a0c292db74fc2e"]
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
    vhash = ["aa5122d77b6b4c25c0e5960cf982c81f06aa5d53", "ac181fe1d82fdcacd26bd41b11ac9c49f9b1dbdb", "6cfee7330ccc8450899157b7f70d21a0cfc6c408", "18d4771ce9e5668d24a576436ae69eeaa6befaa8", "60eff50548ed1a5e66e3234443abd68f1849d570", "03e7a8bf15e0774d84843e0d98f7551fca19693e", "73ed530a50f186e71e3de71e4d579bdc9e5299d1", "6b876dfa7c901ccae0a164f540a0c292db74fc2e"]
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
    vhash = ["111d1c783b3dc941f6a7c4a2ea490bdad2c8a1b5", "3d4c6e9ccb17011eec83c51cddbc81994c8fcc3a", "5ee5d2ce26fd22b0daded2ee34f7e1aae36aa036"]
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
    vhash = ["adff000734b6448a50ef8bcd1991a0a1bfd2cb1f"]
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
    # Commented out while BitmapMaps is a troublesome dependency
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
