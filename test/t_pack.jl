# Test glyph intersection.
using Test
using BitmapMapsExtras
using BitmapMapsExtras: PALETTE_GRGB, plot_glyph_given_value!, apply_color_by_coverage!
using BitmapMapsExtras: radial_distance_glyph, GSVector, norm, plot_glyphs!, plot_glyphs
using BitmapMapsExtras: placements_and_values, pack_glyphs!, DirectionOnGrid, BidirectionOnGrid
using BitmapMapsExtras: plot_glyphs_given_values!, 𝐧ₚ!, 𝐊!, coarse_radius_for_plotting
using BitmapMapsExtras: GSTensor, draw_bidirectional_vector_glyph_given_value!
using BitmapMapsExtras: coverage_fitting_image, MMatrix
using BitmapMaps: mark_at!
using Random: MersenneTwister
using BitmapMapsExtras.TestMatrices
!@isdefined(hashstr) && include("common.jl")

@testset begin "Mark circumference of a vector glyph"
    img = [1.3f0 *PALETTE_GRGB[1] for i = 1:200, j=1:200]
    cov = [0f0 for i = 1:200, j=1:200]
    gs = GSVector(multip = 1, color = PALETTE_GRGB[3], maxg = 100)
    ptexist = CartesianIndex(100, 100)
    v = [80.0, -20.0] #(x, y)
    plot_glyph_given_value!(cov, ptexist, v, gs)
    apply_color_by_coverage!(img, cov, gs)
    mark_at!(img, ptexist; side = 1)
    rc = coarse_radius_for_plotting(gs, v)
    for α = 0:0.02:2π
        # Point on rough outer circle
        Δi = -Int(round(rc * sin(α)))
        Δj = Int(round(rc * cos(α)))
        npt = ptexist + CartesianIndex((Δi, Δj))
        mark_at!(img, npt; side = 1)
        # 
        # From pt to ptexist
        u = ptexist.I .- npt.I
        # The angle for the ray is calculated as if we're in a 'y up' c.s.
        αu = atan(-u[1], u[2]) 
        # αu ≈ α + π
        rptexist = radial_distance_glyph(gs, v, αu + π)
        ptint = ptexist - CartesianIndex(Int.(round.(u .* rptexist ./ norm(u))))
        mark_at!(img, ptint)
    end
    @test hashstr(img) == "e5fcf38a334e25a24e08a0e6640cd2244cb42cd7"
end

@testset begin "Vector glyph crashless placements and values"
    ppts = [CartesianIndex(100, 100), 
            CartesianIndex(100, 110)]
    gs = GSVector(multip = 1, maxg = 100)
    function b!(v, x)
        v .= [80.0, -20.0]
        v
    end
    d = DirectionOnGrid(b!, zeros(Float64, 199, 199))
    # Visual check
    img = [1.3 *PALETTE_GRGB[1] for i = 1:200, j=1:200]
    plot_glyphs!(img, d, ppts, gs)
    # Only the first glyph is accepted (basic overlap)
    @test placements_and_values(d, gs, ppts) == ([ppts[1]], [[80.0, -20.0]])
    #
    ppts = [CartesianIndex(100, 100), 
            CartesianIndex(100, 130)]
    img = [1.3 *PALETTE_GRGB[1] for i = 1:200, j=1:200]
    plot_glyphs_given_values!(img, ppts, [[80.0, -20.0], [80.0, -20.0]], gs)
    # Only the first glyph is accepted (glancing overlap)
    @test length(placements_and_values(d, gs, ppts)[1]) == 1
    #
    ppts = [CartesianIndex(100, 100), 
            CartesianIndex(100, 140)]
    img = [1.3 *PALETTE_GRGB[1] for i = 1:200, j=1:200]
    plot_glyphs!(img, d, ppts, gs)
    # Both glyphs are accepted (but the order may well be another than in ppts)
    @test length(placements_and_values(d, gs, ppts)[1]) == 2
    @test length(placements_and_values(d, gs, ppts)[1]) == 2
    @test placements_and_values(d, gs, ppts)[2][1] == [80, -20]
end

@testset begin "pack vector glyphs no collision"
    img = [1.0 *PALETTE_GRGB[1] for i = 1:200, j=1:200]
    b(x, y) = [80.0, -20.0]
    @test_throws MethodError pack_glyphs!(img, b, GSVector())
    d = DirectionOnGrid(𝐧ₚ!, z_paraboloid()[400:599, 400:599])
    @test d(3, 99) == [-0.0037305447382605097, -0.36559338434952804]
    gs = GSVector(multip = 30 / norm(d(9,99)), maxg = 100, color = PALETTE_GRGB[2])
    pack_glyphs!(img, d, gs)
    @test hashstr(img) == "628abe11f35d486c112bab65f4178e4ed96c42ff"
end


@testset begin "Mark circumference of a bidirectional vector glyph"
    img = [1.3f0 *PALETTE_GRGB[1] for i = 1:200, j=1:200]
    gs = GSTensor(multip = 0.8, directions = 1, maxg = 200, ming = -200)
    cov = coverage_fitting_image(img, gs)
    @test cov isa Matrix{Float32} # since gs.directions = 1
    ptexist = CartesianIndex(100, 100)
    v = [80.0, 20.0] #(x, y)
    K = MMatrix{2, 2, Float64}(hcat(v, v))
    draw_bidirectional_vector_glyph_given_value!(cov, ptexist, v .* gs.multip, gs.strength)
    apply_color_by_coverage!(img, cov, gs)
    mark_at!(img, ptexist; side = 1)
    rc = coarse_radius_for_plotting(gs, K)
    for α = 0:0.005:2π
        # Point on rough outer circle
        Δi = -Int(round(rc * sin(α)))
        Δj = Int(round(rc * cos(α)))
        npt = ptexist + CartesianIndex((Δi, Δj))
        mark_at!(img, npt; side = 1)
        # 
        # From pt to ptexist
        u = ptexist.I .- npt.I
        # The angle for the ray is calculated as if we're in a 'y up' c.s.
        αu = atan(-u[1], u[2]) 
        # αu ≈ α + π
        rptexist = radial_distance_glyph(gs, v, αu + π)
        ptint = ptexist - CartesianIndex(Int.(round.(u .* rptexist ./ norm(u))))
        mark_at!(img, ptint)
    end
    img
    @test hashstr(img) == "0a38f1e16fc787176f57a8e2fa53b36bbe20171c"
end

@testset begin "Mark circumference of a tensor glyph"
    img = [1.3f0 *PALETTE_GRGB[1] for i = 1:500, j=1:500]
    gs = GSTensor(multip = 2, directions = 1:2, maxg = 500, ming = -500)
    cov = coverage_fitting_image(img, gs)
    @test cov isa Vector{Matrix{Float32}} # since gs.directions = 1:2
    ptexist = CartesianIndex(250, 250)
    K = MMatrix{2, 2, Float64}([80.0 25.0; 10.0 -85.0])
    plot_glyph_given_value!(cov, ptexist, K, gs)
    apply_color_by_coverage!(img, cov, gs)
    mark_at!(img, ptexist; side = 1)
    rc = coarse_radius_for_plotting(gs, K)
    for α = 0:0.005:2π
        # Point on rough outer circle
        Δi = -Int(round(rc * sin(α)))
        Δj = Int(round(rc * cos(α)))
        npt = ptexist + CartesianIndex((Δi, Δj))
        mark_at!(img, npt; side = 1)
        # 
        # From pt to ptexist
        u = ptexist.I .- npt.I
        # The angle for the ray is calculated as if we're in a 'y up' c.s.
        αu = atan(-u[1], u[2]) 
        # αu ≈ α + π
        rptexist = radial_distance_glyph(gs, K, αu + π)
        ptint = ptexist - CartesianIndex(Int.(round.(u .* rptexist ./ norm(u))))
        mark_at!(img, ptint)
    end
    img
    @test hashstr(img) == "0a6fd852795a8ec742c34589c2cc07310ea52886"
end

@testset begin "Tensor glyph crash detection"
    ppts = [CartesianIndex(100, 100), 
            CartesianIndex(100, 110)]
    gs = GSTensor(multip = 0.75, maxg = 100, ming = -100)
    function d!(K, args...)
        K .= [80.0 25.0; 10.0 -85.0]
        K
    end
    d = BidirectionOnGrid(d!, zeros(Float64, 199, 199))
    # Visual check
    img = [1.3f0 *PALETTE_GRGB[1] for i = 1:200, j=1:200]
    plot_glyphs!(img, d, ppts, gs)
    # Only the first glyph is accepted (both directions crashing)
    @test placements_and_values(d, gs, ppts) == ([ppts[1]], [[80.0 25.0; 10.0 -85.0]])
    #
    ppts = [CartesianIndex(100, 65), 
            CartesianIndex(90, 155)]
    img = [1.3 *PALETTE_GRGB[1] for i = 1:200, j=1:200]
    plot_glyphs!(img, d, ppts, gs)
    # Only the first glyph is accepted (primary direction)
    @test length(placements_and_values(d, gs, ppts)[1]) == 1
    #
    img = [1.3f0 *PALETTE_GRGB[1] for i = 1:300, j=1:300]
    d = BidirectionOnGrid(d!, zeros(Float64, 299, 299))
    ppts = [CartesianIndex(80, 125), 
            CartesianIndex(217, 160)]
    plot_glyphs!(img, d, ppts, gs)
    # Only the first glyph is accepted (secondary direction crash)
    @test length(placements_and_values(d, gs, ppts)[1]) == 1
    # 
    img = [1.3f0 *PALETTE_GRGB[1] for i = 1:300, j=1:300]
    ppts = [CartesianIndex(140, 125), 
            CartesianIndex(140, 185)]
    plot_glyphs!(img, d, ppts, gs)
    # Primary axis crossing over secondary
    @test length(placements_and_values(d, gs, ppts)[1]) == 1
end


@testset begin "2d vector glyph crash detection, 1 ends within 2"
    gs = GSVector(multip = 1.16, maxg = 200)
    # Here, we check if the crash test works with the
    # functor machinery, especially getting coordinate systems
    # right
    function fdir!(v, M)
        α = M[3, 3]
        # Output is interpreted as (x, y) oriented
        v .= [70 * cos(α ), 70 * sin(α)]
        v
    end
    ppts = [CartesianIndex(50, 50), 
            CartesianIndex(50, 150)]
    z = zeros(Float64, 199, 199)
    z[ppts[1]] = -36.6 * π / 180
    z[ppts[2]] = -127 * π / 180
    dog = DirectionOnGrid(fdir!, z)
    @test dog(ppts[1].I...) ≈ [56.19722326337802, -41.73574124759311]
    @test dog(ppts[2].I...) ≈ [-42.127051620643385, -55.90448570331049]
    # Visual check
    img = [1.3f0 * PALETTE_GRGB[1] for i = 1:200, j=1:200]
    plot_glyphs!(img, dog, ppts, gs)
    @test hashstr(img) == "d37f369eb6716b1a6e073430645451ec57507065"
    # 
    @test length(placements_and_values(dog, gs, ppts)[1]) == 1
end



#@testset begin "2d vector glyph crash detection, 1 crosses 2"
    gs = GSVector(multip = 1.16, maxg = 200)
    # Here, we check if the crash test works with the
    # functor machinery, especially getting coordinate systems
    # right
    function fdir!(v, M)
        α = M[3, 3]
        # Output is interpreted as (x, y) oriented
        v .= [70 * cos(α ), 70 * sin(α)]
        v
    end
    ppts = [CartesianIndex(70, 70), 
            CartesianIndex(50, 150)]
    z = zeros(Float64, 199, 199)
    z[ppts[1]] = -36.6 * π / 180
    z[ppts[2]] = -127 * π / 180
    dog = DirectionOnGrid(fdir!, z)
    # Visual check
    img = [1.3f0 * PALETTE_GRGB[1] for i = 1:200, j=1:200]
    plot_glyphs!(img, dog, ppts, gs)
    @test hashstr(img) == "71a3e75095473159943ed6d9121dffbb4a09da04"
    # 
    @test length(placements_and_values(dog, gs, ppts)[1]) == 1
#end























@testset begin "Bidirectional vector glyph crash detection"
    gs = GSTensor(multip = 1, maxg = 100, ming = -100, directions = 1)
    function fdir!(K, a, b, c, d, M, e, f)
        α = M[3, 3]
        K .= [-85 * sin(α) -85 * cos(α ); 65 * sin(α + 1.5) 65 * cos(α + 1.5)]
        K
    end
    ppts = [CartesianIndex(65, 65), 
            CartesianIndex(100, 150)]
    z = zeros(Float64, 199, 199)
    z[ppts[1]] = -36.6 * π / 180
    z[ppts[2]] = -127 * π / 180
    bdog = BidirectionOnGrid(fdir!, z)
    # 
    img = [1.3f0 * PALETTE_GRGB[1] for i = 1:200, j=1:200]
    plot_glyphs!(img, bdog, ppts, gs)
    # 
    @test length(placements_and_values(bdog, gs, ppts)[1]) == 1
end

@testset begin "Bidirectional tensor glyph crash detection"
    gs = GSTensor(multip = 1, maxg = 100, ming = -100)
    function fdir!(K, a, b, c, d, M, e, f)
        α = M[3, 3]
        K .= [-85 * sin(α) -85 * cos(α ); 65 * sin(α + 1.5) 65 * cos(α + 1.5)]
        K
    end
    ppts = [CartesianIndex(65, 65), 
            CartesianIndex(100, 150)]
    z = zeros(Float64, 199, 199)
    z[ppts[1]] = -36.6 * π / 180
    z[ppts[2]] = -127 * π / 180
    bdog = BidirectionOnGrid(fdir!, z)
    # 
    img = [1.3f0 * PALETTE_GRGB[1] for i = 1:200, j=1:200]
    plot_glyphs!(img, bdog, ppts, gs)
    # 
    @test length(placements_and_values(bdog, gs, ppts)[1]) == 1
end




#@testset begin "pack tensor glyphs no collision"
    # TODO: TOO FEW?
    n = 399
    img = [RGB{N0f8}(0.9, 0.9, 0.4) for i = 1:n, j=1:n]
    b = BidirectionOnGrid(𝐊!, z_paraboloid()[400:(400 + n), 400:(400 + n)])
    gs = GSTensor(multip = 10 / norm(b(9, n ÷ 4)), maxg = n, ming = -n)
    pack_glyphs!(img, b, gs)
    @test hashstr(img) == "f5a0001e4bf5b8659fa42daf9d59a9ab8c04143b"
#end


using BitmapMapsExtras: potential_scattered_placements
n = 68
img = [RGB{N0f8}(0.9, 0.9, 0.4) for i = 1:n, j=1:n]
b = BidirectionOnGrid(𝐊!, z_paraboloid()[400:(400 + n), 400:(400 + n)])
gs = GSTensor(multip = 10 / norm(b(9, n ÷ 4)), maxg = n, ming = -n)
ppts = potential_scattered_placements(b, gs)
# Too many are accepted... (23) -> ( 18) -> ( 16)
filtered_placements, filtered_values = placements_and_values(b, gs, ppts)
badpts = filter(filtered_placements) do pt
    norm(pt.I .- (52, 22)) < 3 && return true
    norm(pt.I .- (59, 25)) < 3 && return true
    false
end
img = [RGB{N0f8}(0.9, 0.9, 0.4) for i = 1:n, j=1:n]
plot_glyphs!(img, b, badpts, gs)
filtered_p, filtered_v = placements_and_values(b, gs, badpts)









b = BidirectionOnGrid(𝐊!, z_cos(; mult = 50)[400:600, 400:600])
img = background(z_cos(; mult = 50)[400:600, 400:600])
#length(passed_placements) = 11
#  0.010389 seconds (90.91 k allocations: 5.896 MiB)
@time pack_glyphs!(img, b, GSTensor(multip = 15000))






















# DEV
img = background(z_paraboloid())
d = DirectionOnGrid(𝐧ₚ!, z_paraboloid())
gs = GSVector(multip = 100 / norm(d(9,99)), maxg = 100, color = PALETTE_GRGB[3])
#  1.580872 seconds (311.57 k allocations: 79.697 MiB, 1.27% gc time) with 
# 2153 passed
# 1.631167 seconds (311.32 k allocations: 72.143 MiB, 0.62% gc time)
# 1884 passed
# 1.680551 seconds (311.12 k allocations: 66.954 MiB, 0.70% gc time)
# 1695 passed
# 1.831165 seconds (310.79 k allocations: 58.223 MiB, 7.67% gc time)
# 1358 passed
# What's missing is a check along the 'main axis'.
# 1298 passed
#1.804302 seconds (310.73 k allocations: 56.585 MiB, 0.71% gc time)
# 1228 passed
# 1.862029 seconds (310.66 k allocations: 54.603 MiB, 0.28% gc time)
# length(passed_placements) = 1228
#  1.518991 seconds (1.30 k allocations: 15.044 MiB, 0.10% gc time)
# 1.524559 seconds (1.29 k allocations: 15.044 MiB, 0.36% gc time)
@time pack_glyphs!(img, d, gs)
@profview pack_glyphs!(img, d, gs)






# Type stable placements_and_values
img = background(z_paraboloid())
d = DirectionOnGrid(𝐧ₚ!, z_paraboloid())
gs = GSVector(multip = 100 / norm(d(9,99)), maxg = 100, color = PALETTE_GRGB[3])
ppts = [CartesianIndex(100, 100),
        CartesianIndex(110, 100), 
        CartesianIndex(100, 140)]
#length(passed_placements) = 2
#  0.000198 seconds (34 allocations: 1.586 KiB)
#  0.000221 seconds (24 allocations: 1.258 KiB)
@time placements_and_values(d, gs, ppts)
@code_warntype placements_and_values(d, gs, ppts)
@inferred placements_and_values(d, gs, ppts)



# 1.520 s (1293 allocations: 15.04 MiB)
@btime pack_glyphs!(img, d, gs)
gs1 = GSVector(multip = 10 / norm(d(9,99)), maxg = 100, color = PALETTE_GRGB[3])
img = background(z_paraboloid())
@time pack_glyphs!(img, d, gs1)
@profview_allocs  pack_glyphs!(img, d, gs1)