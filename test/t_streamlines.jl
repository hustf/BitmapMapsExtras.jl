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
    vhash = ["6c3b833da3a0cee35fa076414b1da634247f51d9", "a060f0b8c67c8c6747a052e264794ea2fa9a0ee9"]
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
    vhash =  ["7f7d44aa5486b0247ea3040618522f8ef0fd1238", "3c5b19a9df2cba1c5b58bd08f539452a04e40135"]
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
    vhash = ["a5dd5f4a9ec2058c9c91747232e0501ddb3a8ff3", "5edbfcccee28b502039a6f9be9b21ca3d72f81b8", "c4288308518d82734246b180fdf227b355e79557", "795c4ae66d13e8fdc199c34b1c88d1c43dcb1763", "0805fe7355fdeb26d210899725ed0f61a93ee182", "23321207f7815bb4fb723e83a90e85df1120b5dc", "93d90e879393339722160ed52540898fd4d20cde", "3c5b19a9df2cba1c5b58bd08f539452a04e40135"] 
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
    vhash =  ["7826964a3761926ef6d136820a0a488300f768b8", "45ab184d1308335c349ffdad01a64c8afe29a822", "ac0bd7d1798c752524bb8ad19234fe197f27cfa4"]
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
    vhash =  String[]
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
