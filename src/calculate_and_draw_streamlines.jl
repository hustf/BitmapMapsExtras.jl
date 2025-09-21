# We're using the heavy dependency and abstraction OrdinaryDiffEq
# because it allows us to interpolate within the solution.

##################
# Draw streamlines
##################

"""
    plot_streamlines(f, z, pts; 
        sol_density = 0.05, r = 1.2f0, 
        strength = 0.7f0, rgb = COLOR_CURVGLYPH,
        odekws...)
"""
function plot_streamlines(f, z, pts; 
        sol_density = 0.05, r = 1.2f0,
        strength = 0.7f0, rgb = COLOR_CURVGLYPH,
        odekws...)
    # Allocate an empty color image (since user didn't supply one)
    img = zeros(RGBA{N0f8}, size(z)...)
    # Modify the image
    plot_streamlines!(f, z, pts, img; sol_density, r, strength, rgb)
end

"""
    plot_streamlines!(img, z, pts, f; 
        sol_density = 0.05, r = 1.2f0,
        strength = 0.7f0, rgb = COLOR_CURVGLYPH),
        odekws...)
"""
function plot_streamlines!(img, z, pts, f; 
    sol_density = 0.05, r = 1.2f0,
    strength = 0.7f0, rgb = COLOR_CURVGLYPH, 
    odekws...)
    # Coverage buffer
    cov = zeros(Float32, size(img)...)
    # Modify coverage
    plot_streamlines!(cov, z, pts, f; sol_density, r, strength, odekws...)
    # Apply color by coverage
    apply_color_by_coverage!(img, cov, rgb)
    img
end

"""
    plot_streamlines!(cov::Array{Float32}, f, z, pts; 
        sol_density = 0.05, r = 1.2f0, strength = 0.7f0,
        odekws...)

# Keyword arguments

- `sol_density` - How far along a streamline's `t` parameter we move before checking which pixel this point 
    would match. Reduce this to plot all covered pixels. If Δt ≈ 1 pixel width, a value of 10 would result in
    pixels visually spaced apart.
- `r` - This applies after the touched pixels are identified. Around each pixel, apply a 
   sprayed circle with radius r. If r > 3, the centre pixels are not sprayed.
- `strength` The nominal coverage applied in each spray (coverage tapers off). 
   Coverage accumulates 'infinitely' if a pixel is sprayed 'infinitely' many times. After all 
    spraying along streamlines is finished, coverage will be converted non-linearly to color 
    application, see `apply_color_by_coverage`. Hence, small `strength` makes streamlines less 
    visible, unless where many streamlines overlap.
"""
function plot_streamlines!(cov::Array{Float32}, z, pts, f; 
        sol_density = 0.05, r = 1.2f0, strength = 0.7f0,
        odekws...)
    # Indexed points on the streamlines
    spts = get_streamlines_points(f, z, pts, sol_density; odekws...)
    # Plot points
    draw_streamlines_points!(cov, spts, r, strength)
    cov
end


"""
    plot_bidirec_streamlines(f, z, pts, primary::Bool, flip::Bool; 
        sol_density = 0.05, r = 1.2f0, 
        strength = 0.7f0, rgb = COLOR_CURVGLYPH,
        odekws...)
"""
function plot_bidirec_streamlines(f, z, pts, primary::Bool, flip::Bool; 
        sol_density = 0.05, r = 1.2f0,
        strength = 0.7f0, rgb = COLOR_CURVGLYPH,
        odekws...)
    # Allocate an empty color image (since user didn't supply one)
    img = zeros(RGBA{N0f8}, size(z)...)
    # Modify the image
    plot_bidirec_streamlines!(img, z, pts, f, primary, flip; sol_density, r, strength, rgb, odekws...)
end

"""
    plot_bidirec_streamlines!(img, z, pts, f, primary::Bool, flip::Bool; 
        sol_density = 0.05, r = 1.2f0,
        strength = 0.7f0, rgb = COLOR_CURVGLYPH),
        odekws...)
"""
function plot_bidirec_streamlines!(img, z, pts, f, primary::Bool, flip::Bool; 
    sol_density = 0.05, r = 1.2f0,
    strength = 0.7f0, rgb = COLOR_CURVGLYPH, 
    odekws...)
    # Coverage buffer
    cov = zeros(Float32, size(img)...)
    # Modify coverage
    plot_bidirec_streamlines!(cov, f, z, pts, primary, flip; sol_density, r, strength, odekws...)
    # Apply color by coverage
    apply_color_by_coverage!(img, cov, rgb)
    img
end

"""
    plot_bidirec_streamlines!(cov::Array{Float32}, f, z, pts, primary::Bool, flip::Bool; 
        sol_density = 0.05, r = 1.2f0, strength = 0.7f0,
        odekws...)
"""
function plot_bidirec_streamlines!(cov::Array{Float32}, f, z, pts, primary::Bool, flip::Bool; 
        sol_density = 0.05, r = 1.2f0, strength = 0.7f0,
        odekws...)
    # Indexed points on the streamlines
    spts = get_bidirec_streamlines_points(f, z, pts, sol_density, primary, flip; odekws...)
    # Plot points
    draw_streamlines_points!(cov, spts, r, strength)
    cov
end


############################
# Callees for OrdinaryDiffEq
############################



# ODE right-hand side, not bidirectional
function rhs!(du, u, p::DirectionAtXY, t)
    v = p(u[1], u[2])
    du[1] = v[1]
    du[2] = v[2]
    du
end

# ODE right-hand side
function rhs!(du, u, uxy::UnidirectionAtXY, t)
    # Calling uxy mutates uxy.du and returns the value
    mag = norm(uxy(u[1], u[2]))
    # Normalize uxy.du to unit length or zero length
    if mag < MAG_EPS
        #@debug "zero mag"
        uxy.du .= 0.0
    else
        uxy.du ./= mag
    end
    du .= uxy.du
end
# Exit criterion 
function signed_distance_within_domain(u, t, integrator::ODEIntegrator)
    @assert integrator.p isa DirectionFunctor
    signed_distance_within_domain(integrator.p, u[1], u[2])
end
signed_distance_within_domain(uxy::UnidirectionAtXY, x, y) =  signed_distance_within_domain(uxy.baxy.d, x, y)
signed_distance_within_domain(fxy::DirectionFunctor, x, y) =  signed_distance_within_domain(fxy.d, x, y)

# Exit criterion. The hard coded criterion here is not
# intended for fine tuning. Instead, set a threshold 
# when defining the integrator functions (directional
# vector set to zero when under that threshold)
function too_flat(u, t, integrator::ODEIntegrator)
    out = max(-1.0, norm(u[1], u[2]) - 0.008)
    if out < 0f0
        #@debug "tooflat $out   "
    end
    out
end 

# Our own function for potential debugging termination causes.
affect!(integrator) = terminate!(integrator)


function condition_flip_bidirection(u, t, integrator::ODEIntegrator)
    @assert t > integrator.tprev
    is_pointing_roughly_opposite(integrator.p, integrator.uprev, u, t)
end
function is_pointing_roughly_opposite(uxy::UnidirectionAtXY, u0, u1, t)
    dprod = dot_product_with_previous(uxy, u0, u1)
    # Direction change is 
    # 138° < direction change <  228°
    if dprod < - 0.668
        return true
    end
    false
end

function affect_flip_bidirection!(integrator)
    integrator.p.flip[] = ! integrator.p.flip[] 
    integrator.p.du .=- integrator.p.du
    nothing
end

function condition_swap_primary_secondary(u, t, integrator::ODEIntegrator)
    @assert t > integrator.tprev
    # TODO: This repeats some evaluations. Perhaps better share one discrete callback?
    is_pointing_roughly_perpendicular(integrator.p, integrator.uprev, u, t)
end
function affect_swap_primary_secondary!(integrator)
    uxy = integrator.p
    uxy.baxy.primary[] = ! uxy.baxy.primary[]
    # Anyway, we have a 50% chance of picking the right flip direction.
    # Perhaps it would be better to join the conditions into one, and later
    # call an affect with predestined flip and primary fields.



    #= Has no effect on the solution....
    # Modify which vector is taken from BidirectionAtXY into UnidirectionAtXY
    if uxy.baxy.primary[]
        if uxy.flip[]
            uxy.du .=- uxy.baxy.K[:, 1]
        else
            uxy.du .= uxy.baxy.K[:, 1]
        end
    else
        if uxy.flip[]
            uxy.du .=- uxy.baxy.K[:, 2]
        else
            uxy.du .= uxy.baxy.K[:, 2]
        end
    end
    # Normalize uxy.du to unit length or zero length
    mag = norm(uxy.du)
    if mag < MAG_EPS
        uxy.du .= 0.0
    else
        uxy.du ./= mag
    end
    =#
    nothing
end


function is_pointing_roughly_perpendicular(uxy::UnidirectionAtXY, u0, u1, t)
    dprod = dot_product_with_previous(uxy, u0, u1)
    # Direction change is 
    # +/-48° < direction change <  +/- 138° 
    if -0.668 < dprod < 0.668
        #print("prim@$t  ")
        return true
    end
    false
end

"""
    dot_product_with_previous(uxy::UnidirectionAtXY, u0, u1)

Returns 1.0 if the dot product can't be evaluated.
"""
function dot_product_with_previous(uxy::UnidirectionAtXY, u0, u1)
    # Vector from last solution point
    Δu = (u1 - u0)
    magΔ = norm(Δu)
    magΔ < MAG_EPS && return 1.0
    # ....normalized to unit length
    Δu ./= magΔ
    # Normalized differential du1 
    du = uxy.du
    norm(du) < 0.98 && return 1.0
    # Dot product    
    Δu[1] * du[1] + Δu[2] * du[2]
end

#################################################################
# Prepare and solve the differential equation for each streamline
#################################################################

function get_bidirec_streamlines_points(f, z, pts, sol_density, primary::Bool, flip::Bool;
     odekws...)
    #
    @assert eltype(pts) <: CartesianIndex{2}
    # Cache, function, etc. f or integration
    uxy = UnidirectionAtXY(f, z, primary, flip)
    # Find solutions, i.e. streamlines
    sols = get_streamlines_xy(uxy, pts; odekws...)
    # NegateY (function for i --> y and for y --> i)
    negy = uxy.baxy.negy
    # Extract indexed points from streamlines
    map(sol -> extract_discrete_points_on_streamline(sol, negy, sol_density), sols)
end


function get_streamlines_points(f, z, pts, sol_density; 
    odekws...)
    #
    @assert eltype(pts) <: CartesianIndex{2}
    daxy = DirectionAtXY(f, z)
    # NegateY (for i --> y and for y --> i)
    negy = daxy.negy
    # Find solutions, i.e. streamlines
    sols = get_streamlines_xy(daxy, pts; odekws...)
    # Extract indexed points from streamlines
    map(sol -> extract_discrete_points_on_streamline(sol, negy, sol_density), sols)
end

function get_streamlines_xy(fxy::T, pts; odekws...) where T<:DirectionFunctor
    @assert eltype(pts) <: CartesianIndex{2}
    # Start coordinates in (x, y)
    vu0 = vu0_from_pts(fxy, pts)
    # Get the streamlines in (x, y)
    get_solution_xy(fxy, vu0; odekws...)
end


"""
    vu0_from_pts(baxy::BidirectionAtXY, pts)
    vu0_from_pts(uxy::UnidirectionAtXY, pts)
    ---> 
"""
vu0_from_pts(fxy::T, pts) where T<:DirectionFunctor = vu0_from_pts(fxy.negy, pts)
vu0_from_pts(uxy::UnidirectionAtXY, pts) = vu0_from_pts(uxy.baxy.negy, pts)
function vu0_from_pts(negy::NegateY, pts)
    @assert eltype(pts) <: CartesianIndex{2}
    vx = map(pt -> float(pt.I[2]), pts)
    vy = map(pt -> negy(float(pt.I[1])), pts)
    [MVector{2, Float64}([x, y]) for (x, y) in zip(vx, vy)]
end



function get_solution_xy(fxy::T, vu0; odekws...) where T<:DirectionFunctor
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
function add_discrete_callbacks!(vdcb, uxy::UnidirectionAtXY)
    # Flip selected direction callback
    push!(vdcb, DiscreteCallback(condition_flip_bidirection, affect_flip_bidirection!, save_positions=(true,true)))
    # Flip primary <--> secondary
    push!(vdcb, DiscreteCallback(condition_swap_primary_secondary, affect_swap_primary_secondary!, save_positions=(true,true)))
end
add_discrete_callbacks!(vdcb, fxy) = vdcb 

function get_solution_xy(daxy, vu0; odekws...)
    throw("now dead")
    # We need to define a stopping point. 
    # Take it from keywords if supplied. 
    tspan = (0.0, get(odekws, :tstop, 1000))
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
    # Flip selected direction callback
    push!(vdcb, DiscreteCallback(condition_flip_bidirection, affect_flip_bidirection!, save_positions=(true,true)))
    cbs = CallbackSet(vccb..., vdcb...)
    # Drop the 'already spent' keywords 'tstop' and 'dtfloor'.
    # The remaining keywords will be passed on.
    remaining_kws = filter(odekws) do (kw, kwval)
        kw == :tstop && return false
        kw == :dtfloor && return false
        true
    end
    # We're re-using daxy here, for many points. 
    # Unfortunately, daxy may change during solution of one.
    # Stil, we don't define as if the state of daxy is part of the initial condition.

    # TODO
    if isempty(remaining_kws)
        map(vu0) do u0
            println("\n u0 = $u0:")
            println("\n daxy(u0...) = $(daxy(u0...)):")
            solve_in_xy(u0, tspan, daxy, cbs)
        end
    else
        map(vu0) do u0
            println("\n u0 = $u0:")
            solve_in_xy(u0, tspan, daxy, cbs; remaining_kws...)
        end
    end
end


function solve_in_xy(u0::Vector{Float64}, tspan, daxy, cbs; odekws...) 
    prob = ODEProblem(rhs!, u0, tspan, daxy, callback = cbs; odekws...)
    # Initialization evaluates rhs!
    integrator = init(prob, Tsit5())
    solve!(integrator)
end


function solve_ensemble(uxy, vu0, tspan, cbs; odekws...)
    u0 = MVector{2, Float64}(first(vu0))
    function prob_func(prob, i, repeat)
        reset!(prob.p)
        u0 = MVector{2, Float64}(vu0[i])
        remake(prob, u0 = u0)
    end
    prob = ODEProblem(rhs!, u0, tspan, uxy, callback = cbs)
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    # EnsembleThreads: Early results indicated tests varying between runs.
    # Change to EnsembleSerial for testing, slows down some.
    # TEMP 
    sim = solve(ensemble_prob, Tsit5(), EnsembleSerial(), trajectories = length(vu0); odekws...)
    sim
end


# This function takes in the base problem and modifies it to create the new problem that the trajectory actually solves. 


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







