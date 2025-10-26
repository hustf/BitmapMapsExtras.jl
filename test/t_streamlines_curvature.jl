using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: plot_streamlines!, 𝐊!, 𝐊ᵤ!, Stroke, indices_on_grid
using BitmapMapsExtras: SelectedVec2AtXY
using BitmapMapsExtras: PALETTE_GRGB
import BitmapMaps
using BitmapMaps: mark_at!, divergence_of_gradients, display_if_vscode
import StatsBase
using StatsBase: Weights, sample
import Random
using Random: MersenneTwister


!@isdefined(is_hash_stored) && include("common.jl")

@testset "Streamlines, 𝐊ᵤ! versus less graphically useful 𝐊!" begin
    vhash = ["1bc542fd14604a6242058ed4b79807e8c5ece082", "e66a13d3828fb1ee9c5ba86dde93e6d418ec7535"]
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
    vhash = vhash = ["d625d4caa33ddbef09c51072d2f3a05df38b1956", "32ad175e9c45dfce5afbafbe5a35d8467984d6c1", "8a15a49dad200fda49df737f6aed40017a849a49"]
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
    vhash = vhash = ["2edf25f7b1ec183135c03f7fb466f049e5a1c901", "6c435eeba7681eee0244ad9133106b5061e30c86", "db78a7768bd8f5cd0f7844c4cd1f4fa2eada6569", "9c5d765f3654873d7f3c775059cf644f72df05bd", "0310688a14af0a61b2b51525a6a7b14f36f96038", "2ee096ed3b6d37f7541ee38dd1e7f0c743d79046", "3e9a898ee60b93de1c28199c17c64c7cff100791", "68d793cfe9056c56c6881903a4ffe89b211ea071"]
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
    vhash = ["aa6b1645bea451b53748dc8cd1c9f7e044347b67"]
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
    vhash = ["6e95bc2a7aa7bb9e30c8cfede7ee989baeeeb847", "3bf202c19c808e5f7095fa2cb7b96b0d5abe853b"]
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
    vhash = ["231407335af6fa8f3aa37d9c32d7a402c116b59f", "3d28e7a94e0c419462b9d48ef83d22cc17cafa43"]
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