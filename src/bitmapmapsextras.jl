module BitmapMapsExtras
import BitmapMaps
using BitmapMaps: scaleminmax00, line!, display_if_vscode, open_as_temp_in_imgedit
import ImageFiltering
import ImageFiltering.Kernel
import ImageCore
using ImageCore: N0f8, GrayA, RGBA, colormap, Gray, RGB, alpha
import StaticArrays
using StaticArrays: SVector, MVector, MMatrix, SMatrix
import LinearAlgebra
using LinearAlgebra: ⋅, norm, normalize!
import ImageTransformations
using ImageTransformations: recenter, center, RotMatrix, WarpedView
import Interpolations
using Interpolations: OnGrid, Flat, BSplineInterpolation
import LinearSolve
using LinearSolve: init, LinearProblem, solve!, NormalCholeskyFactorization
using LinearSolve: OperatorCondition, OperatorAssumptions
import ColorBlendModes
using ColorBlendModes: BlendMode, BlendLighten, BlendHue, blend

export tangent_basis

"""
Kernel for first derivative.
The 5x5 area is probably not perfectly planar. The naive 
way to approach this is picking the horizontal and 
vertical as representative, which would amount to 
an anisotropic bias. Instead, we're differentiating over
wider path to reduce the effect. 

Note that the first, vertical direction is 'pre-negated', 
flipping the y-axis direction: y is screeen-up, x is to right.
"""
const KERN1´ = -Kernel.ando5()[1]
const KERN2´ = Kernel.ando5()[2]

"""
Kernel for second derivative. Applies to a vector with values interpolated from a matrix and a direction.
"""
const KERN´´ = SVector{5, Float64}([0.25, 0.0, -0.5, 0.0, 0.25])

"Indices for second derivative, from a lazily rotated matrix."
const W = CartesianIndex{2}[CartesianIndex(1, 3), CartesianIndex(2, 3), CartesianIndex(3, 3), CartesianIndex(4, 3), CartesianIndex(5, 3)]

"""
Sample angles in the tangent plane.
Values are optimized in 'optimize_sample_angles.jl'.
"""
const VΦ = MVector{4, Float64}([-7, 7, 87, 93] .* (π / 180 ))

"""
Low-effort definition of distinguishable plane colors in roughly similar luminance.
"""
const RED_GREEN_BLUE = SMatrix{3, 3, N0f8}([ 0.95  0.0  0.0
                                        0.0  0.5  0.0
                                        0.0  0.0  0.9])

"""
Low-effort definition of curvature glyph color, used where overlain a colorful picture.
"""
const COLOR_CURVGLYPH = RGB{N0f8}(0.85, 0.5, 0.9)

"""
Low-effort definition of curvature types color, equal luminosity for good effect with BlendHue.
"""
const PALETTE_GRGB = SVector{4, RGB{N0f8}}(
    [RGB(0.467, 0.467, 0.467), RGB(0.957, 0.0, 0.078), RGB(0.0, 0.549, 0.0), RGB(0.0, 0.443, 1.0)])


include("test_matrices.jl")
include("differential_geom/tangent_basis.jl")
include("differential_geom/curvature.jl")
include("differential_geom/bidirectional_quantity.jl")
include("visualization/draw_direct.jl")
include("visualization/draw_bidirectional_glyph.jl")
include("visualization/draw_plane.jl")
include("visualization/draw_vector_glyph.jl")
include("calculate_and_draw_glyphs.jl")
include("calculate_and_paint_curvature_type.jl")
include("domain.jl")
include("hachures.jl")
include("streamlines.jl")

end
