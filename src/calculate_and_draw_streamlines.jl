# We're using the heavy dependency and abstraction OrdinaryDiffEq
# because it allows us to interpolate within the solution.



##################
# Draw streamlines
##################

"""
    plot_streamlines(fdir, z, pts; 
        tstop = 100000.0, nsamples = 10000)
"""
function plot_streamlines(fdir, z, pts; 
    tstop = 100000.0, nsamples = 10000)
    # Allocate an empty color image (since user didn't supply one)
    img = zeros(RGBA{N0f8}, size(z)...)
    # Modify the image
    plot_streamlines!(img, fdir, z, pts; tstop)
end

"""
    plot_streamlines!(img, fdir, z, pts; 
        tstop = 100000.0, nsamples = 10000, r = 5)
"""
function plot_streamlines!(img, fdir, z, pts; 
    tstop = 100000.0, nsamples = 10000, r = 5)
    # Black-white buffer
    bbuf = Array{GrayA{N0f8}}(falses( size(img)...))
    # Modify the buffer
    plot_streamlines_gray!(bbuf, fdir, z, pts; tstop, nsamples, r)
    # Function converting GrayA{N0f8} to color.
    # CONSIDER TODO: Separate color constant for vectors.
    #     Also, use mutable containers for color definitions.
    fcol = x -> RGBA{N0f8}(COLOR_CURVGLYPH.r, COLOR_CURVGLYPH.g, COLOR_CURVGLYPH.b, 
                            x.val)
    # Composite bbuf over img 
    map!(BlendLighten, img, img, fcol.(bbuf))
    img
end

"""
    plot_streamlines_gray!(bbuf, fdir, z, pts; 
        tstop = 100000.0, nsamples = 10000, r = 5)
"""
function plot_streamlines_gray!(bbuf::Array{GrayA{N0f8}}, fdir, z, pts; 
        tstop = 100000.0, nsamples = 10000, r = 2)
    # Indexed points on the streamlines
    spts = get_streamlines_points(fdir, z, pts, tstop, nsamples)
    # Plot points
    plot_streamline_gray!(bbuf, CartesianIndices(z), spts, r)
    bbuf
end

"""
    plot_streamline_gray!(bbuf, R, spts, r)
"""
function plot_streamline_gray!(bbuf, R, spts, r)
    for streamline in spts
        for pt in streamline
            spray_neighbors!(bbuf, R, pt, r)
        end
    end
    bbuf
end

############################
# Callees for OrdinaryDiffEq
############################

# ODE right-hand side
function rhs!(du, u, p, t)
    v = p(u[1], u[2])
    du[1] = v[1]
    du[2] = v[2]
    du
end

# Exit criterion 
function signed_distance_within_domain(u, t, integrator)
    @assert integrator.p isa DirectionAtXY
    d = integrator.p.d
    signed_distance_within_domain(d, u[1], u[2])
end

# Exit criterion 
function too_flat(u, t, integrator)
    @assert integrator.p isa DirectionAtXY
    too_flat = max(-1.0, norm(integrator.p(u[1], u[2])) - 0.008)
    too_flat
end 


#################################################################
# Prepare and solve the differential equation for each streamline
#################################################################

function get_streamlines_points(fdir, z, pts::CartesianIndices, tstop, nsamples)
    # Integrator object 
    daxy = DirectionAtXY(fdir, z)
    # NegateY (for i --> y and for y --> i)
    negy = daxy.ny
    # Find solutions, i.e. streamlines
    sols = get_streamlines_xy(daxy, pts::CartesianIndices, tstop)
    # Extract indexed points from streamlines
    map(sol -> extract_discrete_points_on_streamline(sol, negy, nsamples), sols)
end

function get_streamlines_xy(daxy, pts::CartesianIndices, tstop)
    # Start coordinates
    vx = map(pt -> float(pt.I[2]), pts)
    vy = map(pt -> daxy.ny(float(pt.I[1])), pts)
    # A vector of initial conditions vectors. 
    vu0 = [[x, y] for (x, y) in zip(vx, vy)]
    # Get the streamlines in (x, y)
    get_solution_xy(daxy, vu0, tstop)
end

function get_solution_xy(daxy, vu0, tstop)
    # End of line (if not sooner)
    tspan = (0.0, tstop)
    # Other end conditions
    affect!(integrator) = terminate!(integrator)
    cb = CallbackSet(
        ContinuousCallback(signed_distance_within_domain, affect!),
        ContinuousCallback(too_flat, affect!)
    )
    map(u0 -> solve_in_xy(u0, tspan, daxy, cb), vu0)
end

solve_in_xy(u0::Vector{Float64}, tspan, daxy, cb) = solve(ODEProblem(rhs!, u0, tspan, daxy), Tsit5(), callback = cb)

################################################
# Sample streamlines from continuous to discrete  
################################################

function extract_discrete_points_on_streamline(sol, ny::NegateY, nsamples)
    if ! (sol.retcode == Success || sol.retcode == Terminated)
        @show sol.retcode
        throw("What's up?")
    end
    # Pre-allocate
    pts = CartesianIndex{2}[]
    oldi, oldj = 0, 0, 0
    for t in range(first(sol.t), last(sol.t), length = nsamples)
        x, y = sol(t) # This is a fast, interpolated lookup 
        negy = ny(y)
        i = Int(round(negy))
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







