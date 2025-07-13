# This contains NegateY and Domain,
# which helps with index to coordinate transformations.
# These are used from trace.jl, but could be implemented
# for drawing glyps as well. This approach is better suited
# for solving differntial equations, since we have to pretend
# that z is a continuous function instead of a matrix.


"""
Map image row index y (1.0 at top) to Cartesian-style y (1.0 at bottom).
When instantiated, can be called for conversions i -> y and y -> i.
"""
struct NegateY
    originy::Int
end
NegateY(R::AbstractMatrix) = NegateY(size(R, 1) + 1)
# Callable:
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

