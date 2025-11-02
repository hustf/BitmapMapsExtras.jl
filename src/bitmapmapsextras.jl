module BitmapMapsExtras
import Base
import Base: size
import DrawAndSpray
using DrawAndSpray: line!, spray_along_nested_indices!, draw_bidirectional_vector!
using DrawAndSpray: is_bidirec_vect_positive, spray!, display_if_vscode
using DrawAndSpray: draw_vector!, apply_color_by_any_coverage!, mark_at!
import DrawAndSpray: apply_color_by_coverage!, over!, chromaticity_over!
import ImageFiltering
import ImageFiltering.Kernel
import ImageCore
using ImageCore: N0f8, GrayA, RGBA, colormap, Gray, RGB, Lab, alpha, Colorant
using ImageCore: channelview
import StaticArrays
using StaticArrays: SVector, MVector, MMatrix, SMatrix
import LinearAlgebra
using LinearAlgebra: ⋅, norm, normalize!
import LinearSolve
using LinearSolve: init, LinearProblem, solve!, NormalCholeskyFactorization
using LinearSolve: OperatorCondition, OperatorAssumptions, LinearCache
import ColorSchemes
using ColorSchemes: ColorScheme
import Random
using Random: MersenneTwister, randperm
import Unicode
using Unicode: graphemes

# For streamlines
import OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, ContinuousCallback, DiscreteCallback, CallbackSet
using OrdinaryDiffEq: terminate!, Tsit5, solve, get_proposed_dt, init, solve!
using OrdinaryDiffEq: remake, EnsembleThreads,  EnsembleSerial, EnsembleProblem
import OrdinaryDiffEqCore
using OrdinaryDiffEqCore: ODEIntegrator
import SciMLBase
using SciMLBase: ReturnCode.Success, ReturnCode.Terminated, ReturnCode.DtLessThanMin
#
export tangent_basis # TODO Add more at some later time, or use new public thing.

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

"Kernel for second derivative. Applies to a vector with values interpolated from a matrix and a direction."
const KERN´´ = SVector{5, Float64}([0.25, 0.0, -0.5, 0.0, 0.25])

"Sample angles in the tangent plane. Values were optimized, ref. `optimize_sample_angles.jl`."
const VΦ = MVector{4, Float64}([-7, 7, 87, 93] .* (π / 180 ))


"Low-effort definition of curvature glyph color, used where overlain a colorful picture."
const COLOR_CURVGLYPH = RGB{N0f8}(0.85, 0.5, 0.9)


"Tinted grey, red, green, blue."
const PALETTE_GRGB = SVector{4, RGB{N0f8}}(
    [RGB(0.4785,0.5119,0.5096), RGB(0.957, 0.0, 0.078), RGB(0.0, 0.549, 0.0), RGB(0.0, 0.443, 1.0)])

"Bright colors in roughly similar luminance."
const PALETTE_RGB = (RGB{N0f8}(0.95, 0.0,  0.0),
                     RGB{N0f8}(0.0,  0.5,  0.0),
                     RGB{N0f8}(0.0,  0.0,  0.9))
                                            
"Close to grayscale, maximizes color distance to PALETTE_GRGB"
const PALETTE_BACKGROUND = ColorScheme([
    RGB{N0f8}(0.004, 0.0, 0.008)
    RGB{N0f8}(0.149, 0.075, 0.02)
    RGB{N0f8}(0.325, 0.137, 0.004)
    RGB{N0f8}(0.49, 0.216, 0.004)
    RGB{N0f8}(0.612, 0.29, 0.043)
    RGB{N0f8}(0.769, 0.373, 0.051)
    RGB{N0f8}(0.969, 0.475, 0.016)
    RGB{N0f8}(0.976, 0.518, 0.165)
    RGB{N0f8}(0.976, 0.588, 0.329)
    RGB{N0f8}(0.953, 0.678, 0.518)
        ])

"Limit for meaningful vector normalization."
const MAG_EPS = √(eps(Float64)) 

"Fixed taper for vector glyphs, including bidirectional vectors"
const VECTOR_REL_HALFWIDTH = 0.075

"Type alias"
const TENSORMAP = MMatrix{2, 2, Float64, 4}

"Direction 1 to be shown"
const D1  = 0x01
"Direction 2 to be shown"
const D2  = 0x02
"Bitwise: Specifies whether direction 1 and / or 2 are to be shown"
const D12 = D1 | D2


include("differential_geom/tangent_basis.jl")
include("differential_geom/curvature.jl")
include("differential_geom/bidirectional_quantity.jl")
include("type_definitions/show.jl")
include("type_definitions/domain_types.jl")
include("type_definitions/direction_types.jl")
include("type_definitions/bidirection_types.jl")
include("type_definitions/specification_types.jl")
include("visualization/draw_plane.jl")
include("visualization/plot_glyph_given_value.jl")
include("visualization/pack_glyphs.jl")
include("visualization/calculate_and_draw_glyphs.jl")
include("visualization/calculate_and_draw_streamlines.jl")
include("visualization/calculate_and_paint_curvature_type.jl")
include("visualization/default_ij_functor.jl")
include("visualization/glyph_indices.jl")
include("tensormap_functions.jl")
include("Test_matrices.jl")

# Function aliases
const 𝐊 = principal_curvature_components
const 𝐊! = principal_curvature_components!
const 𝐊ᵤ! = principal_curvature_normalized!
const 𝐧ₚ = descent
const 𝐧ₚ! = descent!
const 𝐧ₚᵤ! = descent_unit!

end
