using Test
import ImageFiltering
import ImageFiltering.Kernel
using BitmapMaps
using BitmapMaps: Gray, GrayA, RGBA, scaleminmax00, N0f8, gray, alpha, mapwindow
using BitmapMaps: composite_on_top!, CompositeDestinationOver
using BitmapMaps: line!
using Colors: colormap
import LinearAlgebra
using LinearAlgebra: ⋅, ×, norm
using StaticArrays
using BenchmarkTools
includet("utilties_test_matrix.jl")
using .TestMatrices
using .TestMatrices: r
using OffsetArrays

include("utilties_graphic.jl");

includet("utilties_tangent_5.jl")

# Load 5x5 M, visualize M
Ω = CartesianIndices((-2:2, -2:2))
pt = TestMatrices.I0 + CartesianIndex(-306, 0)
# Take the test array, visualize it
M = z_sphere()[pt .+ Ω]

bluesc(M; mi = 392.0, ma = 397.)

tangent_basis(M);
# 180.857 ns (5 allocations: 304 bytes)
# 189.571 ns (5 allocations: 304 bytes)
# 298.227 ns (17 allocations: 976 bytes)
@btime tangent_basis(M)
# Allocate in a closure...
ftb = let v = MVector{2,Float64}(Array{Float64, 1}(undef, 2)),
    P = MMatrix{3,3,Float64}(Array{Float64, 2}(undef, (3, 3)))
    M -> tangent_basis!(P, v, M)
end
ftb(M)
# 167.492 ms (135633 allocations: 13.03 MiB)
# 208.002 ms (7120179 allocations: 256.64 MiB)
#  224.933 ms (8118181 allocations: 287.10 MiB)
#   226.166 ms (8116713 allocations: 287.05 MiB)
# 358.139 ms (12105785 allocations: 469.66 MiB)
# 342.782 ms (12105785 allocations: 469.66 MiB)
@btime mapwindow(ftb, $(z_sphere()), (5,5));


let v = MVector{2,Float64}(Array{Float64, 1}(undef, 2)),
    P = MMatrix{3,3,Float64}(Array{Float64, 2}(undef, (3, 3)))
    #
    @code_warntype tangent_basis!(P, v, M)
end


function steepness_color(M)
   bluesc(steepness(M); mi = 0.0, ma = 1.0)
end
function northness_color(M)
   bluesc(northness(M); mi = -1.0, ma = 1.0)
end
function eastness_color(M)
   bluesc(eastness(M); mi = -1.0, ma = 1.0)
end

# Takes 15 sec each, but that's due to not doing this the smart way by re-using memory.
mapwindow(steepness_color, z_sphere(), (5,5))
mapwindow(northness_color, z_sphere(), (5,5))
mapwindow(eastness_color, z_sphere(), (5,5))


