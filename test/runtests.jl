using Test
#using BitmapMapsExtras
#using BitmapMapsExtras.TestMatrices
@testset "differential_geom" begin
    @testset "tangent_basis" begin
        include("t_tangent_basis.jl")
    end
    @testset "curvature" begin
        include("t_curvature.jl")
    end
    @testset "calculate and draw glyphs" begin
        include("t_calculate_and_draw_glyphs.jl")
    end
    @testset "paint curvature  type" begin
        include("t_paint_curvature_type.jl")
    end
end

@testset "differential equations" begin
    @testset "domain" begin
        include("t_domain_types.jl")
    end
    @testset "direction types" begin
        include("t_direction_types.jl")
    end
    @testset "streamlines" begin
        include("t_streamlines.jl")
    end
    # Not included: t_streamlines_curvature.jl
end