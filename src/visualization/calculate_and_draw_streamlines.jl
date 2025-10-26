# We're using the heavy dependency and abstraction OrdinaryDiffEq
# because it allows us to interpolate within the solution.

##################
# Plot streamlines
##################
"""
    plot_streamlines(fxy::AbstractXYFunctor, pts; 
        stroke = Stroke(), odekws...)
    --> Matrix{RGBA{N0f8}}

# Arguments

- `fxy` determines which curve to follow, and also the size of the output image. 
  See ``
- `pts` is a vector of starting points. See `indices_on_grid`.
- `stroke` specifies the appearance of each streamline.
- `odekws` is a list of keywords for the ordinary differential equation. `dtmax = 1.0` is an example: A smaller timestep can solve many problems.

# Example, 2d vector field

```
julia> using BitmapMapsExtras: TestMatrices

julia> vaxy = Vec2AtXY(𝐧ₚ!, z_cylinder(π / 6));

julia> img = background(vaxy)

julia> img = plot_streamlines!(img, vaxy, pts; stroke = Stroke(color = PALETTE_GRGB[2]), dtmax = 1.0)
```

# Example, 2d vector selected from tensor map



"""
function plot_streamlines(fxy::AbstractXYFunctor, pts; 
        stroke = Stroke(), odekws...)
    # Allocate an empty color image (since user didn't supply one)
    img = zeros(RGBA{N0f8}, size(fxy)...)
    # Modify the image
    plot_streamlines!(img, fxy, pts; stroke, odekws...)
end

"""
    plot_streamlines!(img, fxy::AbstractXYFunctor, pts; 
        stroke = Stroke(), odekws...)

See `plot_streamlines`. Use this method when you already have an image.
"""
function plot_streamlines!(img, fxy::AbstractXYFunctor, pts; 
    stroke = Stroke(), odekws...)
    # Coverage buffer
    cov = zeros(Float32, size(img)...)
    # Modify coverage
    plot_streamlines!(cov, fxy, pts; stroke, odekws...)
    # Apply color by coverage
    apply_color_by_coverage!(img, cov, stroke.color)
    img
end



"""
    plot_streamlines!(cov::Array{Float32}, fxy::AbstractXYFunctor, pts; 
        stroke = Stroke(), odekws...)
"""
function plot_streamlines!(cov::Array{Float32}, fxy::AbstractXYFunctor, pts; 
        stroke = Stroke(), odekws...)
    # Indexed points on the streamlines
    spts = get_streamlines_points(fxy, pts, stroke.tdensity; odekws...)
    # Plot points
    draw_streamlines_points!(cov, spts, stroke.r, stroke.strength)
    cov
end


############################
# Callees for OrdinaryDiffEq
############################

# ODE right-hand side
function rhs!(du, u, fxy!::AbstractXYFunctor, t)
    # Calling fxy mutates fxy.v and returns the value
    du .= fxy!(u[1], u[2])
end

# Exit criterion 
function signed_distance_within_domain(u, t, integrator::ODEIntegrator)
    @assert integrator.p isa AbstractXYFunctor
    signed_distance_within_domain(integrator.p, u[1], u[2])
end
signed_distance_within_domain(saxy::SelectedVec2AtXY, x, y) =  signed_distance_within_domain(saxy.baxy.d, x, y)
signed_distance_within_domain(fxy::AbstractXYFunctor, x, y) =  signed_distance_within_domain(fxy.d, x, y)

# Exit criterion. The hard coded criterion here is not
# intended for fine tuning. Instead, set a threshold 
# when defining the integrator functions (directional
# vector set to zero when under that threshold)
function too_flat(u, t, integrator::ODEIntegrator)
    out = max(-1.0, norm(u[1], u[2]) - 0.008)
    #if out < 0f0
        #@debug "tooflat $out   "
    #end
    out
end 

# Our own function for potential debugging termination causes.
affect!(integrator) = terminate!(integrator)


function condition_flip_bidirection(u, t, integrator::ODEIntegrator)
    @assert t > integrator.tprev
    is_pointing_roughly_opposite(integrator.p, integrator.uprev, u, t)
end
function is_pointing_roughly_opposite(saxy::SelectedVec2AtXY, u0, u1, t)
    dprod = dot_product_with_previous(saxy, u0, u1)
    # Direction change is 
    # 138° < direction change <  228°
    if dprod < - 0.668
        return true
    end
    false
end

function affect_flip_bidirection!(integrator)
    integrator.p.flip[] = ! integrator.p.flip[] 
    integrator.p.v .=- integrator.p.v
    nothing
end

function condition_swap_primary_secondary(u, t, integrator::ODEIntegrator)
    @assert t > integrator.tprev
    # TODO: This repeats some evaluations. Perhaps better share one discrete callback?
    is_pointing_roughly_perpendicular(integrator.p, integrator.uprev, u)
end
function affect_swap_primary_secondary!(integrator)
    saxy = integrator.p
    saxy.baxy.primary[] = ! saxy.baxy.primary[]
    nothing
end


# This needs some prettying up.
# See 'dot_product_with_previous' and 'is_close_to_perpendicular'.
function is_pointing_roughly_perpendicular(saxy::SelectedVec2AtXY, u0, u1)
    dprod = dot_product_with_previous(saxy, u0, u1)
    # Direction change is 
    # +/-48° < direction change <  +/- 138° 
    if -0.668 < dprod < 0.668
        return true
    end
    false
end

"""
    dot_product_with_previous(saxy::SelectedVec2AtXY, u0, u1)

Returns 1.0 if the dot product can't be evaluated.
"""
function dot_product_with_previous(saxy::SelectedVec2AtXY, u0, u1)
    # Vector from last solution point
    Δu = (u1 - u0)
    magΔ = norm(Δu)
    magΔ < MAG_EPS && return 1.0
    # ....normalized to unit length
    Δu ./= magΔ
    # Normalized differential du1 
    du = saxy.v
    norm(du) < 0.98 && return 1.0
    # Dot product    
    Δu[1] * du[1] + Δu[2] * du[2]
end

#################################################################
# Prepare and solve the differential equation for each streamline
#################################################################

function get_streamlines_points(fxy::AbstractXYFunctor, pts, sol_density;
     odekws...)
    #
    @assert eltype(pts) <: CartesianIndex{2}
    # Cache, function, etc. f or integration
    #saxy = SelectedVec2AtXY(f, z, primary, flip)
    # Find solutions, i.e. streamlines
    sols = get_streamlines_xy(fxy, pts; odekws...)
    # NegateY (function for i --> y and for y --> i)
    negy = fxy.negy
    # Extract indexed points from streamlines
    map(sol -> extract_discrete_points_on_streamline(sol, negy, sol_density), sols)
end



function get_streamlines_xy(fxy::T, pts; odekws...) where T<:AbstractXYFunctor
    @assert eltype(pts) <: CartesianIndex{2}
    # Start coordinates in (x, y)
    vu0 = vu0_from_pts(fxy, pts)
    # Get the streamlines in (x, y)
    get_solution_xy(fxy, vu0; odekws...)
end


"""
    vu0_from_pts(baxy::BidirectionAtXY, pts)
    vu0_from_pts(saxy::SelectedVec2AtXY, pts)
    ---> 
"""
vu0_from_pts(fxy::T, pts) where T<:AbstractXYFunctor = vu0_from_pts(fxy.negy, pts)
vu0_from_pts(saxy::SelectedVec2AtXY, pts) = vu0_from_pts(saxy.baxy.negy, pts)
function vu0_from_pts(negy::NegateY, pts)
    @assert eltype(pts) <: CartesianIndex{2}
    vx = map(pt -> float(pt.I[2]), pts)
    vy = map(pt -> negy(float(pt.I[1])), pts)
    [MVector{2, Float64}([x, y]) for (x, y) in zip(vx, vy)]
end



function get_solution_xy(fxy::T, vu0; odekws...) where T<:AbstractXYFunctor
    # We need to define a stopping point. 
    # Take it from keywords if supplied. 
    tspan = make_tspan(;odekws...)
    cbs = callbacks_streamlines(fxy; odekws...)
    # Drop the 'already spent' keywords comprising 'tspan'.
    # The remaining keywords will be passed on to the solver.
    remaining_kws = filter(odekws) do (kw, kwval)
        kw == :tstart && return false
        kw == :tstop && return false
        kw == :dtfloor && return false
        true
    end
    solve_ensemble(fxy, vu0, tspan, cbs; remaining_kws...)
end
function make_tspan(; odekws...)
    if :tspan ∈ keys(odekws)
        @assert :tstop ∉ keys(odekws) "Don't specify both tspan and tstop"
        return odekws[:tspan]
    else
        return (get(odekws, :tstart, 0), get(odekws, :tstop, 1000))
    end
end
function callbacks_streamlines(fxy; odekws...)
    # Early end conditions
    vccb = [ContinuousCallback(signed_distance_within_domain, affect!),
            ContinuousCallback(too_flat, affect!)]
    vdcb = DiscreteCallback[]
    if :dtfloor ∈ keys(odekws)
        push!(vdcb, let dtfloor = odekws[:dtfloor]
                DiscreteCallback(
                   (u, t, integrator ) -> get_proposed_dt(integrator) ≤ dtfloor, 
                       affect!)
            end)
    end
    add_discrete_callbacks!(vdcb, fxy)
    CallbackSet(vccb..., vdcb...)
end
function add_discrete_callbacks!(vdcb, saxy::SelectedVec2AtXY)
    # Flip selected direction callback
    push!(vdcb, DiscreteCallback(condition_flip_bidirection, affect_flip_bidirection!, save_positions=(true,true)))
    # Flip primary <--> secondary
    push!(vdcb, DiscreteCallback(condition_swap_primary_secondary, affect_swap_primary_secondary!, save_positions=(true,true)))
end
add_discrete_callbacks!(vdcb, fxy) = vdcb 

function solve_ensemble(saxy, vu0, tspan, cbs; odekws...)
    u0 = MVector{2, Float64}(first(vu0))
    function prob_func(prob, i, repeat)
        reset!(prob.p)
        u0 = MVector{2, Float64}(vu0[i])
        remake(prob, u0 = u0)
    end
    prob = ODEProblem(rhs!, u0, tspan, saxy, callback = cbs)
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    # EnsembleThreads: Early results indicated tests varying between runs.
    # We changed to EnsembleSerial. This slows down the calculation.
    solve(ensemble_prob, Tsit5(), EnsembleSerial(), trajectories = length(vu0); odekws...)
end

################################################
# Sample streamlines from continuous to discrete  
################################################

function extract_discrete_points_on_streamline(sol, negy::NegateY, sol_density)
    if ! (sol.retcode == Success ||
        sol.retcode == Terminated) # || sol.retcode == DtLessThanMin)
        @show sol.retcode
        throw("What's up?")
    end
    # Pre-allocate
    pts = CartesianIndex{2}[]
    oldi, oldj = 0, 0, 0
    # sol_density determines how far, in solution time, between examined solution points.
    # Examined points are then checked for uniqueness (discrete pixels).
    trng = range(first(sol.t), last(sol.t), step = sol_density * sign( last(sol.t) - first(sol.t)   ))
    for t in trng
        x, y = sol(t) # This is a fast, interpolated lookup 
        ny = negy(y)
        i = Int(round(ny))
        j = Int(round(x))
        if i !== oldi || j !== oldj
            # We didn't just visit this pixel before.
            p = CartesianIndex(i, j)
            push!(pts, p)
            oldi, oldj = i, j
         end
    end
    pts
end
