# Next abstraction after dev_trace_1.j
#    Draw streamlines from gradients and from 
#    principal curvature.
# This drop the closures, dropping the unnecessary LocalInterpolator
# For an ODE, it's better to pack all we need in a struct.
using Test
import ImageFiltering
import ImageFiltering.Kernel
import BitmapMaps
using BitmapMaps: scaleminmax00, line!, display_if_vscode, open_as_temp_in_imgedit
import ImageTransformations
using ImageTransformations: recenter, center, RotMatrix, WarpedView
import Interpolations
using Interpolations: OnGrid, Flat, BSplineInterpolation
import ImageCore
using ImageCore: N0f8, GrayA, RGBA, colormap, Gray, RGB
import ColorBlendModes
using ColorBlendModes: BlendLighten
import ColorTypes
using ColorTypes: alpha
import LinearAlgebra
using LinearAlgebra: ⋅, ×, norm, normalize!
using StaticArrays
using BenchmarkTools
import LinearSolve
using LinearSolve: init, solve!, LinearProblem, NormalCholeskyFactorization
using LinearSolve: OperatorCondition, OperatorAssumptions
import Interpolations
using Interpolations: interpolate, BSpline, Linear
import OffsetArrays
using OffsetArrays: OffsetMatrix
######################
# Utilties
######################
includet("utilties_test_matrix.jl")
using .TestMatrices
using .TestMatrices: r, I0
includet("utilties_graphic_6.jl")
includet("utilties_tangent_8.jl")
includet("utilties_curvature_7.jl")
includet("plot_basis_curvature_8.jl")

"""
Sample angles in the tangent plane.
Values are optimized in 'optimize_sample_angles.jl'.
"""
const VΦ = MVector{4, Float64}([-7, 7, 87, 93] .* (π / 180 ))

# Map image row index y (1.0 at top) to Cartesian-style y (1.0 at bottom)
struct NegateY
    originy::Int
end
NegateY(R::AbstractMatrix) = NegateY(size(R, 1) + 1)
(ny::NegateY)(y) = ny.originy - y

# True if (x, y) lies within rectangular bounds defined by CartesianIndices,
# with y measured from bottom (1) upwards

"""
    Callable: returns true if (x, y) is inside the domain
"""
struct Domain
    minx::Float64
    miny::Float64
    maxx::Float64
    maxy::Float64
end
# Constructor, applies to coordinates (x, y) where y is up
function Domain(R::CartesianIndices)
    flipy = NegateY(R)
    minx = float(R[1][2])
    miny = float(flipy(R[end][1]))
    maxx = float(R[end][2])
    maxy = float(flipy(R[1][1]))
    Domain(minx, miny, maxx, maxy)
end
function Domain(R::CartesianIndices, Ω::CartesianIndices)
    flipy = NegateY(R)
    minx = float(R[1][2] - Ω[1][2])
    miny = float(flipy(R[end][1] - Ω[end][2]))
    maxx = float(R[end][2] - Ω[end][2])
    maxy = float(flipy(R[1][1] - Ω[1][2]))
    Domain(minx, miny, maxx, maxy)
end
# Callable: returns true if (x, y) is inside the domain
(d::Domain)(x, y) = d.minx ≤ x ≤ d.maxx && d.miny ≤ y ≤ d.maxy

"""
    signed_distance_within_domain(d::Domain, x, y)

Returns the distance to the closest border of the domain.
If outside domain, returns a negative number which is the
distance to the closest border.
"""
function signed_distance_within_domain(d::Domain, x, y)
    # Calculate distances to each boundary
    dist_left = x - d.minx
    dist_right = d.maxx - x
    dist_bottom = y - d.miny
    dist_top = d.maxy - y
    # Find the minimum distance to any boundary
    min_dist = min(dist_left, dist_right, dist_bottom, dist_top)
end


struct DirectionOnGrid
    Ω::CartesianIndices # Square neighborhood
    fdir!::Function # e.g 𝐧ₚ!(v, M)
    z::Matrix
end
# Constructor
function DirectionOnGrid(fdir!, z)
    Ω = CartesianIndices((-2:2, -2:2))
    DirectionOnGrid(Ω, fdir!, z)
end
# Direction as function of valid Cartesian 
direction_on_grid!(v, dog::DirectionOnGrid, pt::CartesianIndex) =  dog.fdir!(v, view(dog.z, dog.Ω .+ pt))


struct DirectionInDomain
    dog::DirectionOnGrid
    li::BSplineInterpolation{MVector{2, Float64}, 2, OffsetMatrix{MVector{2, Float64}, Matrix{MVector{2, Float64}}}, BSpline{Linear{Interpolations.Throw{OnGrid}}}, Tuple{Base.Slice{UnitRange{Int64}}, Base.Slice{UnitRange{Int64}}}}
end
# Constructor
function DirectionInDomain(fdir!, z) 
    corners_direction = SMatrix{2,2, MVector{2,Float64}}([ MVector(0.0, 0.0) for i in 1:2, j in 1:2] )
    # The interpolation can be modified in-place through li.coefs
    li = interpolate(corners_direction, BSpline(Linear()))
    DirectionInDomain(DirectionOnGrid(fdir!, z), li)
end
# Callable (allocates)
(did::DirectionInDomain)(x, negy) = direction_in_domain!(MVector{2, Float64}([0, 0]), did, x, negy)
# Direction as function of in-domain float coordinates. 'negy' is 0 at bottom of image.
function direction_in_domain!(v, did::DirectionInDomain, x, negy)
    # Identify (up to four) points on grid
    j1, j2 =  Int(floor(x)), Int(ceil(x))
    i1, i2 =  Int(floor(negy)), Int(ceil(negy))
    # Refer to mutable storage location for the directions at grid points.
    cd = did.li.coefs
    # Refer to 
    dog = did.dog
    # Calculate and store closest grid points
    if i1 == i2 && j1 == j2
        # On grid, no interpolation (but we keep the structure for consistency)
        direction_on_grid!(v, dog, CartesianIndex{2}(i1, j1)) 
        cd[1, 1] .= v
        cd[2, 1] .= v
        cd[1, 2] .= v
        cd[2, 2] .= v
    elseif i1 == i2
        # Line interpolation j1 - j2
        direction_on_grid!(v, dog, CartesianIndex{2}(i1, j1)) 
        cd[1, 1] .= v
        cd[2, 1] .= v
        direction_on_grid!(v, dog, CartesianIndex{2}(i1, j2)) 
        cd[1, 2] .= v
        cd[2, 2] .= v
    elseif j1 == j2
        # Line interpolation i1 - i2
        direction_on_grid!(v, dog, CartesianIndex{2}(i1, j1)) 
        cd[1, 1] .= v
        cd[1, 2] .= v
        direction_on_grid!(v, dog, CartesianIndex{2}(i2, j1))
        cd[2, 1] .= v
        cd[2, 2] .= v
    else
        direction_on_grid!(v, dog, CartesianIndex{2}(i1, j1)) 
        cd[1, 1] .= v
        direction_on_grid!(v, dog, CartesianIndex{2}(i2, j1)) 
        cd[2, 1] .= v
        direction_on_grid!(v, dog, CartesianIndex{2}(i1, j2)) 
        cd[1, 2] .= v
        direction_on_grid!(v, dog, CartesianIndex{2}(i2, j2)) 
        cd[2, 2] .= v
    end
    # Interpolate within (1,1)-(2, 2)
    xn = 1 + x - j1
    yn = 1 + negy - i1
    v .= did.li(yn, xn)
    v
end


"""
    DirectionAtXY

The object type we hand to the differential equation solver. 

Callable: DirectionAtXY(fdir!, z)(x, y) -> 2d directional vector
"""
struct DirectionAtXY
    did::DirectionInDomain
    ny::NegateY
    d::Domain
    v::MVector{2, Float64}
end
# Constructor
function DirectionAtXY(fdir!, z)
    did = DirectionInDomain(fdir!, z)
    ny = NegateY(z)
    R = CartesianIndices(z)
    Ω = CartesianIndices((-2:2, -2:2))
    d = Domain(R, Ω)
    v = MVector{2, Float64}([0, 0])
    DirectionAtXY(did, ny, d, v)
end

# Callable. y is positive up on screen. Returns [0.0,0.0]
# when out-of domain (the callback set doesn't immediately
# stop calculation).
function (daxy::DirectionAtXY)(x, y)
    if daxy.d(x, y)
        direction_at_xy!(daxy.v, daxy.did, x, daxy.ny(y))
    else
        # This occurs fairly often, though isn't presented as part of the solution.
        #printstyled("::DirectionAtXY(x, y) == ($x, $y) ∉ $(daxy.d) \n", color =:176)
        daxy.v .= 0.0
        daxy.v
    end
end
function direction_at_xy!(v, did::DirectionInDomain, x, negy)
    direction_in_domain!(v, did, x, negy)
end

#######################
# Integration as an ODE
#######################
# 
# We're using the heavy dependency and abstraction 
# because it allows us to interpolate within the solution.

using OrdinaryDiffEq
# ODE right-hand side
function rhs!(du, u, p, t)
    v = p(u[1], u[2])
    du[1] = v[1]
    du[2] = v[2]
    du
end

# Exit criteria
function signed_distance_within_domain(u, t, integrator)
    @assert integrator.p isa DirectionAtXY
    d = integrator.p.d
    sd = signed_distance_within_domain(d, u[1], u[2])
    if sd < 0
        #println("u = $(round.(u, digits = 5)) , signed_distance_within_domain = $(round(sd, digits = 5))")
    end
    sd
end
function too_flat(u, t, integrator)
    @assert integrator.p isa DirectionAtXY
    too_flat = max(-1.0, norm(integrator.p(u[1], u[2])) - 0.008)
    if too_flat < 0
        #println("u = $(round.(u, digits = 5)) , too_flat = $(round(too_flat, digits = 5))")
    end
    too_flat
end 

so = let 
    daxy = DirectionAtXY(𝐧ₚ!, z_sphere())
    d = daxy.d
    ny = daxy.ny
    tspan = (0.0, 100.0)
    affect!(integrator) = terminate!(integrator)
    cb = CallbackSet(
        ContinuousCallback(signed_distance_within_domain, affect!),
        ContinuousCallback(too_flat, affect!)
    )
    u0 = [503.0, 453.0]  # start point
    prob = ODEProblem(rhs!, u0, tspan, daxy)
    # We could simply pass many initial conditions, too,
    # then pick unique points.
    sol = solve(prob, Tsit5(), callback = cb, dtmin = 0.02)
end


for t in range(first(so.t), last(so.t), length =1000)
        pt = so(t)
       println(pt)
end



# Visualize trace

function test_trace_plot(fdir!, z, Δ,  dashsize)
    # Background: z-values in color
    img = bluesc(z; mi = -float(r) ) .* 1.0
    # Add simple contour lines, too
    Δc = 20
    wc = Δc / 10
    map!(img, z, img) do zz, pix
        mod(zz, Δc) < wc ? RGBA{N0f8}(0.1, 0.1, 0.1, 1.0) : pix 
    end
    test_trace_plot!(img, fdir!,  z, Δ,  dashsize)
end
function test_trace_plot!(img, fdir!, z, Δ,  dashsize)
    # Black-white buffer
    bbuf = Array{GrayA{N0f8}}(falses( size(img)...))
    test_trace_plot!(bbuf, fdir!, z, Δ,  dashsize)
    # Function converting GrayA{N0f8} to proper color
    f = x -> RGBA{N0f8}(COLOR_CURVGLYPH.r, COLOR_CURVGLYPH.g, COLOR_CURVGLYPH.b, 
                            x.val)
    # Composite bbuf over img
    # Overlay bbuf on img in the proper color
    map!(BlendLighten, img, img, f.(bbuf))
    img
end

function test_trace_plot!(bbuf::Array{GrayA{N0f8}}, fdir!, z, Δ, dashsize)
    daxy = DirectionAtXY(fdir!, z)
    d = daxy.d
    ny = daxy.ny
    tspan = (0.0, 1000.0)
    affect!(integrator) = terminate!(integrator)
    cb = CallbackSet(
        ContinuousCallback(signed_distance_within_domain, affect!),
        ContinuousCallback(too_flat, affect!)
    )
    xs = sort((m=(d.minx+d.maxx)/2; union(m:-Δ:d.minx, m+Δ:Δ:d.maxx)))
    ys = sort((m=(d.miny+d.maxy)/2; union(m:-Δ:d.miny, m+Δ:Δ:d.maxy)))
    for x in xs, y in ys
        u0 = [x, y]  # start point
        prob = ODEProblem(rhs!, u0, tspan, daxy)
        # We could simply pass many initial conditions, too,
        # then pick unique points.
        sol = solve(prob, Tsit5(), callback = cb, dtmin = 0.02)
        if so.retcode != SciMLBase.ReturnCode.Success
            @show so.retcode
            throw("what is this?")
        end
        oldi, oldj = 0, 0
        for t in range(first(sol.t), last(sol.t), length = 1000)
            pt = sol(t) # This is a fast, linear lookup
            x, y = pt
            negy = ny(y)
            i = Int(round(negy))
            j = Int(round(x))
            if i !== oldi || j !== oldj
                # We didn't just visit this pixel before.
                p = CartesianIndex(i, j)
                color_neighbors!(bbuf, CartesianIndices(bbuf), p, dashsize)
                oldi, oldj = i, j
            end
       end
    end
    bbuf
end

test_trace_plot(𝐧ₚ!, z_sphere(), 50, 2)
test_trace_plot(𝐧ₚ!, z_ellipsoid(), 50, 2)
test_trace_plot(𝐧ₚ!, z_paraboloid(; a = r, b= 0.5r), 50, 2)
test_trace_plot(𝐧ₚ!, z_paraboloid(; a = -r, b= -0.5r), 50, 2)
test_trace_plot(𝐧ₚ!, z_paraboloid(; a = -2r, b= 2r), 50, 2)
test_trace_plot(𝐧ₚ!, z_paraboloid(; a = 2r, b= -2r), 50, 2)


