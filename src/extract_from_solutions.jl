# Callees for 'spray_streamlines`
# concerning extraction of values from diffeq solutions.
#

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