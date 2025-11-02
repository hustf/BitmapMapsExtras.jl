using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: plot_streamlines!, 𝐊!, 𝐊ᵤ!, Stroke, indices_on_grid
using BitmapMapsExtras: SelectedVec2AtXY
using BitmapMapsExtras: PALETTE_GRGB
using BitmapMapsExtras: mark_at!, display_if_vscode
#import BitmapMaps
#using BitmapMaps: divergence_of_gradients
import StatsBase
using StatsBase: Weights, sample
import Random
using Random: MersenneTwister


!@isdefined(is_hash_stored) && include("common.jl")

@testset "Streamlines, 𝐊ᵤ! versus less graphically useful 𝐊!" begin
    vhash = ["c068c26f74832c14e932dfbc1f36fd534f2fd34b", "f7ee461e478df960bd471ecb30627dd6801b61af"]
    COUNT[] = 0
    saxy = SelectedVec2AtXY(𝐊!, z_cylinder(π / 6), false, true);
    img = background(saxy)
    pts = indices_on_grid(saxy; Δ = 300)
    mark_at!(img, pts, 5, "in_circle")
    # With non-normalized curvature, 𝐊!, the values du / dt are very small.
    # At time step 1.0, the streamline integration may not have left the initial pixel.
    # That's just a lot of sub-pixel interpolation and extremely in-effective.
    plot_streamlines!(img, saxy, pts; stroke = Stroke(color = PALETTE_GRGB[2]), tstop = 1000000)
    @test is_hash_stored(img, vhash)
    # Normalized curvature. The downside is that we can't get an integrated value of curvature, i.e.
    # local inclination, which might be useful for checks but not for finding the streamline.
    saxy = SelectedVec2AtXY(𝐊ᵤ!, z_cylinder(π / 6), false, true);
    img = background(saxy)
    mark_at!(img, pts, 5, "in_circle")
    plot_streamlines!(img, saxy, pts; stroke = Stroke(color = PALETTE_GRGB[2]), dtmax = 1)
    @test is_hash_stored(img, vhash)
end



@testset "Solution keyword example" begin
    vhash = ["635f10d34bab4966441a911c5ade1d9fea0480d9", "5fbe9a19458055d12f9f020771bb4ccf7382f7eb", "2bfdf353db2fe2f8059abdf4640c1dd9b59482cc"]
    COUNT[] = 0
    saxy = SelectedVec2AtXY(𝐊ᵤ!, z_ellipsoid(), true, true)
    pts = indices_on_grid(saxy)
    img = background(saxy)
    mark_at!(img, pts)
    # This takes too long strides and follow the wrong track
    plot_streamlines!(img, saxy, pts)
    @test is_hash_stored(img, vhash)
    #
    img = background(saxy)
    mark_at!(img, pts)
    # Persistent
    plot_streamlines!(img, saxy, pts; dtmax = 1.0)
    @test is_hash_stored(img, vhash)
    #
    # Secondary direction:
    saxy = SelectedVec2AtXY(𝐊ᵤ!, z_ellipsoid(), false, true)
    # Also persistent
    plot_streamlines!(img, saxy, pts; dtmax = 1, stroke = Stroke(color = PALETTE_GRGB[3]))
    @test is_hash_stored(img, vhash)
end


@testset "Streamlines starting on grid" begin
    vhash = ["aabb60f08d5ac576b5d40646b320c5dd16e2ad41", "688d26635a39906a32892bd6f7c01fe656dd7c8e", "e0c0d90f01c9395e5dec88ff9ecb4ceab682b749", "bde19f8adf00995546d9a17652bff0720fdcb50f", "77daa5c49e14502982cc67bcd229119d378d5d98", "32a336a0b575ff4e3c090f38716a3c31d96f6a99", "cf11016d871e6e0c44671a8c422975bf524c1bae", "1fb02db486fe72d7b21596edafb6abd0aa0ee6ec"]
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
        saxy = SelectedVec2AtXY(𝐊ᵤ!, fz(), true, true);
        pts = indices_on_grid(saxy)
        img = background(saxy)
        # dtmax = 1.0 is 'always' a good idea.
        plot_streamlines!(img, saxy, pts; dtmax = 1.0)
        @test is_hash_stored(img, vhash)
    end
end


@testset "Symmetric directions with different appearance" begin
    vhash = ["a8bfc72a759b56129cbc9c6ef4d459d196cf1f55"]
    COUNT[] = 0
    # Short streamlines from grid points. The direction is the first,
    # which is systematic although hard to predict
    saxy = SelectedVec2AtXY(𝐊ᵤ!, z_ridge_peak_valleys(), true, true)
    pts = indices_on_grid(saxy; Δ = 50)
    img = background(saxy; α = 0.4)
    plot_streamlines!(img, saxy, pts; dtmax = 5, tstop = 50)
    mark_at!(img, pts, 5, "in_circle")
    # Short streamlines from grid points. The direction is the second of two.
    saxy = SelectedVec2AtXY(𝐊ᵤ!, z_ridge_peak_valleys(), true, false)
    stroke = Stroke(color = PALETTE_GRGB[3])
    plot_streamlines!(img, saxy, pts; dtmax = 5, tstop = 50, stroke)
    @test is_hash_stored(img, vhash)
end

@testset "Many streamlines in both directions from grid points" begin
    vhash = ["9ce7899a3307ad5be85819d082167e3b1f187e6b", "a815aa874f50e76336112269e7c0812d9366af16"]
    COUNT[] = 0
    stroke = Stroke(strength = 0.05)
    tstop = 2000 # default is 1000
    dtmax = 1
    # Primary curvature direction. There's a tendency in where they collect,
    # but the tendency is masked by the regular grid
    saxy = SelectedVec2AtXY(𝐊ᵤ!, z_ridge_peak_valleys(), true, true)
    pts = indices_on_grid(saxy, Δ = 50) # 361
    img = background(saxy; α = 0.4)
    #mark_at!(img, pts)
    plot_streamlines!(img, saxy, pts; dtmax, tstop, stroke)
    saxy = SelectedVec2AtXY(𝐊ᵤ!, z_ridge_peak_valleys(), true, false)
    plot_streamlines!(img, saxy, pts; dtmax, tstop, stroke)
    @test is_hash_stored(img, vhash)
    #
    # Primary curvature direction. There's an other tendency in where they collect,
    # but the tendency is masked by the regular grid
    saxy = SelectedVec2AtXY(𝐊ᵤ!, z_ridge_peak_valleys(), false, true)
    img = background(saxy; α = 0.4)
    #mark_at!(img, pts)
    plot_streamlines!(img, saxy, pts; dtmax, tstop, stroke)
    saxy = SelectedVec2AtXY(𝐊ᵤ!, z_ridge_peak_valleys(), false, false)
    plot_streamlines!(img, saxy, pts;  dtmax, tstop, stroke)
    @test is_hash_stored(img, vhash)
end

@testset "Streamlines without grid" begin
    vhash = ["dcf4aa68870019a7a368ef45838e1efe2476c02c", "79e36e9e5e1be84759f36a1e0ad90e59f0008070"]
    COUNT[] = 0
    stroke = Stroke(strength = 0.05)
    tstop = 2000 # default is 1000
    dtmax = 1
    z = z_ridge_peak_valleys()
    # Uniform random distribution
    n = 361
    pts = sample(MersenneTwister(124), CartesianIndices(z), n)
    # Primary curvature. Tends to
    # aggregate on ridges
    img = background(z; α = 0.4)
    plot_streamlines!(img, SelectedVec2AtXY(𝐊ᵤ!, z, true, true), pts; dtmax, tstop, stroke)
    plot_streamlines!(img, SelectedVec2AtXY(𝐊ᵤ!, z, true, false), pts; dtmax, tstop, stroke)
    @test is_hash_stored(img, vhash)
    # Secondary curvature. These do NOT aggregate on ridges,
    # but seems to identify lines of no curvature, as well as valley bottoms
    img = background(z; α = 0.4)
    plot_streamlines!(img, SelectedVec2AtXY(𝐊ᵤ!, z, false, true), pts; dtmax, tstop, stroke)
    plot_streamlines!(img, SelectedVec2AtXY(𝐊ᵤ!, z, false, false), pts; dtmax, tstop, stroke)
    @test is_hash_stored(img, vhash)
end





#=
    # Originating at points where z > 0
    seed_dens = clamp.(z, 0.0, 1.0);
    display_if_vscode(background(seed_dens))
    n = 1000
    pts = sample(MersenneTwister(123), CartesianIndices(seed_dens), Weights(vec(seed_dens)), n);
    saxy = SelectedVec2AtXY(𝐊ᵤ!, z, true, true)
    stroke = Stroke(color = PALETTE_GRGB[3], r = 1, strength = 0.3)
    img = background(z)
    mark_at!(img, pts)
    plot_streamlines!(img, saxy, pts; dtmax = 1, stroke)
    saxy = SelectedVec2AtXY(𝐊ᵤ!, z, true, false)
    plot_streamlines!(img, saxy, pts; dtmax = 1, stroke)
    display_if_vscode(img)
    @test is_hash_stored(img, vhash)
    #
    # Originating at 'sources' (convex terrain)
    vals1 = divergence_of_gradients(-z_ridge_peak_valleys())
    background(vals1)
    vals2 = clamp.(vals1, 0.003, 1.0) .- 0.003
    seed_dens = clamp.(vals2 .* 100, 0, 1.0)
    # Drop the extreme area near the centre.
    seed_dens[400:600, 300:700] .= 0.0
    background(seed_dens)
    n = 1000
    pts = sample(MersenneTwister(123), CartesianIndices(seed_dens), Weights(vec(seed_dens)), n);
    img = background(z)
    mark_at!(img, pts)
    display_if_vscode(img)
    plot_streamlines!(img, saxy, pts; dtmax = 1)
    display_if_vscode(img)
    @test is_hash_stored(img, vhash)
end
=#