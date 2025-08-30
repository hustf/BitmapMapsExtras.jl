import Pkg
if ! haskey(Pkg.project().dependencies, "Random") 
    @warn "For interactive tests: Use a test
    environment outside of this package's /src/ and /test/ folder. 
    The interactive test environment should `Pkg|> develop` BitmapMapsExtras, 
    and `Pkg|> add` test dependencies listed in `test/Project.toml.
    In VSCode, manually change the environment to that environment folder."
end
using Test


@testset "BitmapMapsExtras" begin
@testset "draw direct" begin
    @testset "draw_direct" begin
        include("t_draw_direct.jl")
    end
end
@testset "differential geometry" begin
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
    @testset "bidirection types" begin
        include("t_bidirection_types.jl")
    end
    @testset "streamlines" begin
        include("t_streamlines.jl")
    end
    @testset "streamlines curvature" begin
        include("t_streamlines_curvature.jl")
    end
end
end