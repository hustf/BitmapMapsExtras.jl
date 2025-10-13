using Test
using BitmapMapsExtras
using BitmapMapsExtras: LogMapper, N0f8, spray!, display_if_vscode, RGBA
using BitmapMapsExtras: COLOR_CURVGLYPH, apply_color_by_coverage!

# using UnicodePlots
# using BenchmarkTools
!@isdefined(hashstr) && include("common.jl")


####################
# Non-linear mapping
####################
@testset "Non-linear mapping" begin
    vlin = 0f0:0.5f0:10f0
    # Float32 to Float32
    mapper = LogMapper(1f0/log1p(10f0), 0f0)
    mat = hcat(vlin, mapper.(vlin))
    @test mat ≈ Float32[
        0.0  0.0
        0.5  0.169092
        1.0  0.289065
        1.5  0.382123
        2.0  0.458157
        2.5  0.522443
        3.0  0.57813
        3.5  0.627249
        4.0  0.671188
        4.5  0.710935
        5.0  0.747222
        5.5  0.780602
        6.0  0.811508
        6.5  0.84028
        7.0  0.867194
        7.5  0.892477
        8.0  0.916314
        8.5  0.938862
        9.0  0.960253
        9.5  0.9806
        10.0  1.0
        ]
    # Float32 to N0f8 (mantissa 256)
    mapper = LogMapper()
    @test mapper.(vlin) == N0f8[0.0N0f8, 0.169N0f8, 0.286N0f8, 0.38N0f8, 0.455N0f8, 0.522N0f8, 0.576N0f8, 0.624N0f8, 0.671N0f8, 0.71N0f8, 0.745N0f8, 0.776N0f8, 0.808N0f8, 0.839N0f8, 0.863N0f8, 0.89N0f8, 0.914N0f8, 0.937N0f8, 0.957N0f8, 0.976N0f8, 0.996N0f8]
end
##############
# Spray kernel
##############
@testset "Spray kernel" begin
    cov = zeros(Float32, 11, 11)
    spray!(cov, CartesianIndex((6,6)), 5f0, 1f0)
    @test LogMapper().(cov) == N0f8[
        0.0    0.0    0.0    0.0    0.0    0.329  0.0    0.0    0.0    0.0    0.0
        0.0    0.0    0.329  0.478  0.549  0.569  0.549  0.478  0.329  0.0    0.0
        0.0    0.329  0.525  0.627  0.675  0.69   0.675  0.627  0.525  0.329  0.0
        0.0    0.478  0.627  0.706  0.749  0.761  0.749  0.706  0.627  0.478  0.0
        0.0    0.549  0.675  0.749  0.0    0.0    0.0    0.749  0.675  0.549  0.0
        0.329  0.569  0.69   0.761  0.0    0.0    0.0    0.761  0.69   0.569  0.329
        0.0    0.549  0.675  0.749  0.0    0.0    0.0    0.749  0.675  0.549  0.0
        0.0    0.478  0.627  0.706  0.749  0.761  0.749  0.706  0.627  0.478  0.0
        0.0    0.329  0.525  0.627  0.675  0.69   0.675  0.627  0.525  0.329  0.0
        0.0    0.0    0.329  0.478  0.549  0.569  0.549  0.478  0.329  0.0    0.0
        0.0    0.0    0.0    0.0    0.0    0.329  0.0    0.0    0.0    0.0    0.0
        ]
end

@testset "Rising coverage" begin
    vΣ = Float64[]    
    sumold = 0.0f0
    rn = 0.2f0:0.2f0:8f0
    for r in rn
        ru = Int(ceil(r))
        wh = 2 * ru + 1
        cov = zeros(Float32, wh, wh)
        spray!(cov, CartesianIndex((ru + 1, ru + 1)), r, 1f0)
        α = LogMapper().(cov)
        img = map(x -> RGBA(x,x,0.0N0f8,x), α)
        #display_if_vscode(img)
        Σ = sum(α)
        push!(vΣ, Σ / r)
        #println(rpad("r = $r", 20), rpad("sum = $Σ", 30), 
        #        rpad("sum / r = $(Σ / r)", 30))
        @test Σ >=  sumold
        sumold = Σ
    end
    # UnicodePlots:
    #pl = lineplot(rn, vΣ , width = 200)
    #hline!(pl, 2)
    #vline!(pl, 1:8)
end

############
# Cone lines
############
@testset "Cone lines and optional transpacency layer" begin
    img = zeros(RGBA{N0f8}, 62, 160);
    img .= RGBA{N0f8}(1, 0.494, 0.43, 1);
    cov = zeros(Float32, size(img)...);
    for (j, rr) in zip(10:20:390, 1f0:1f0:8)
        Δr = rr / 60
        for i in 10:1:61
            spray!(cov, CartesianIndex((i, j)), rr, 0.3f0)
            rr -= Δr
        end
    end
    apply_color_by_coverage!(img, cov, RGB{N0f8}(1,1,1))
    @test hashstr(img) == "9a5de89cc254cbe0e59e4080f6915bd790b29200"
    #
    img = zeros(RGB{N0f8}, 62, 160);
    img .= RGB{N0f8}(0, 0, 0);
    # COLOR_CURVGLYPH isa RGB{N0f8}
    apply_color_by_coverage!(img, cov, COLOR_CURVGLYPH)
    @test img isa Matrix{<:RGB}
    @test hashstr(img) == "4fdebbde7a4cd0cfa71af8bd1add54ebc474aa91"
end

