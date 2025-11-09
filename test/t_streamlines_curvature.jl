using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: plot_streamlines!, 𝐊!, 𝐊ᵤ!, Stroke, indices_on_grid
using BitmapMapsExtras: SelectedVec2AtXY
using BitmapMapsExtras: PALETTE_GRGB, MersenneTwister
using BitmapMapsExtras: mark_at!, display_if_vscode
import StatsBase
using StatsBase: Weights, sample

!@isdefined(is_hash_stored) && include("common.jl")

@testset "Streamlines, 𝐊ᵤ! versus less graphically useful 𝐊!" begin
    vhash = ["e64b6b2dde7e0904139a838b324c5b3cfe100709", "c126c32ec2e5a06bcd7b15599c0e74b2744d9fbb"]
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
    vhash = ["7248a93f9ec457d8087e33fdca522325ed6cc6fd", "6ec1f7b870cefa772791f3631ff650c5db9207c1", "edd8be4537dfc249656d29e142bfdf3b73c39419"]
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
    vhash = ["ad8e4694820a831cc97ace21b05e2e8873769da0", "c146dd70bf34342e7cbe73e8ed9d60f32f61d3d0", "fcdff1b08a433fbbe09a7feb7349bbddcd9a32d8", "07db5cfc8785612fda59d3fc9587b219322f8324", "9af5944069742124977d530d4378573678d682d7", "f3264a0c84ea20f7d3aec9936e6767c82da59185", "ca3e974448d579743b22e7c195488443a69b77e3", "2fff0a82074bd2e793f08e8ab5fe09a007c8714a"]
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
    vhash = ["bd9ad3dc16da30cab07c0ee3f78e51fd555bff8b"]
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
    vhash = ["15910c6e2dbee78782692f2f9ed8545006320691", "c482dccf5505289d3be7f7f7cd6062869b135bad"]
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
    vhash = ["ecfab2979c884c4db8a71e854889bfe7b7aa3ae4", "614a7695db192fa73192eb1538f03b5cff49c5c2"]
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




#= Commented out while BitmapMaps is a troublesome dependency
    z = z_ridge_peak_valleys()

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