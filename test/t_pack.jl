# Test glyph intersection.
using Test
using BitmapMapsExtras
using BitmapMapsExtras: PALETTE_GRGB, plot_glyph_given_value!, apply_color_by_coverage!
using BitmapMapsExtras: intersection_distance_glyph, GSVector, norm, plot_glyphs!, plot_glyphs
using BitmapMapsExtras: placements_and_values, pack_glyphs!, DirectionOnGrid, 𝐧ₚ!
using BitmapMapsExtras: plot_glyphs_given_values!
using BitmapMaps: mark_at!
using BitmapMapsExtras.TestMatrices
!@isdefined(hashstr) && include("common.jl")

@testset begin "Circumference of a vector"
    img = [1.3 *PALETTE_GRGB[1] for i = 1:200, j=1:200]
    cov = [0f0 for i = 1:200, j=1:200]
    gs = GSVector(multip = 1, color = PALETTE_GRGB[3], maxg = 100)
    ptexist = CartesianIndex(100, 100)
    v = [80.0, -20.0] #(x, y)
    plot_glyph_given_value!(cov, ptexist, v, gs)
    apply_color_by_coverage!(img, cov, gs.color)
    mark_at!(img, ptexist; side = 1)
    for α = 0:0.02:2π
        Δi = -Int(round(99 * sin(α)))
        Δj = Int(round(99 * cos(α)))
        npt = ptexist + CartesianIndex((Δi, Δj))
        mark_at!(img, npt; side = 1)
        # 
        # From pt to ptexist
        u = ptexist.I .- npt.I
        # The angle is calculated as if we're in a 'y up' c.s.
        αu = atan(-u[1], u[2]) 
        rptexist = intersection_distance_glyph(gs, v, αu + π)
        ptint = ptexist - CartesianIndex(Int.(round.(u .* rptexist ./ norm(u))))
        mark_at!(img, ptint)
    end
    @test hashstr(img) == "4541cb301cdd3e0e45585f58144ef0331b1c6d61"
end

@testset begin "glyph placements and values"
    ppts = [CartesianIndex(100, 100), 
            CartesianIndex(100, 110)]
    gs = GSVector(multip = 1, maxg = 100)
    b(x, y) = [80.0, -20.0]
    # Only the first glyph is accepted
    @test placements_and_values(b, gs, ppts) == ([ppts[1]], [[80.0, -20.0]])
    #
    ppts = [CartesianIndex(100, 100), 
            CartesianIndex(100, 130)]
    img = [1.3 *PALETTE_GRGB[1] for i = 1:200, j=1:200]
    plot_glyphs_given_values!(img, ppts, [[80.0, -20.0], [80.0, -20.0]], gs)
    # Only the first glyph is accepted (glancing overlap)
    @test length(placements_and_values(b, gs, ppts)[1]) == 1
    #
    ppts = [CartesianIndex(100, 100), 
            CartesianIndex(100, 140)]
    img = [1.3 *PALETTE_GRGB[1] for i = 1:200, j=1:200]
    plot_glyphs_given_values!(img, ppts, [[80.0, -20.0], [80.0, -20.0]], gs)
    # Both glyphs are accepted (but the order is likely another than in ppts)
    @test length(placements_and_values(b, gs, ppts)[1]) == 2
    @test length(placements_and_values(b, gs, ppts)[1]) == 2
    @test placements_and_values(b, gs, ppts)[2][1] == b(0, 0)
end

@testset begin "pack glyphs no collision"
    img = [1.0 *PALETTE_GRGB[1] for i = 1:200, j=1:200]
    b(x, y) = [80.0, -20.0]
    @test_throws MethodError pack_glyphs!(img, b, gs)
    d = DirectionOnGrid(𝐧ₚ!, z_paraboloid()[400:599, 400:599])
    @test d(3, 99) == [-0.0037305447382605097, -0.36559338434952804]
    gs = GSVector(multip = 30 / norm(d(9,99)), maxg = 100, color = PALETTE_GRGB[2])
    pack_glyphs!(img, d, gs)
    @test hashstr(img) == "628abe11f35d486c112bab65f4178e4ed96c42ff"
end

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