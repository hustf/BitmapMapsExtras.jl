 # This builds on domain, is a requisite for
 # streamlines.jl,  and define types for
 # 'non-discretization' of matrices, making them
 # work like continuous functions.
 #
 # Defines types:

 # DirectionOnGrid
 # DirectionInDomain
 # DirectionAtXY 
 #  - the top-level abstraction containing the above as well as Domain and NegateY


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
function DirectionInDomain(fdir!, z::AbstractMatrix{<:Number})
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
    struct DirectionAtXY
        did::DirectionInDomain
        ny::NegateY
        d::Domain
        v::MVector{2, Float64}
    end

The object type we hand to the differential equation solver. 

Callable: DirectionAtXY(fdir!, z)(x, y) -> 2d directional vector

This is a top-level abstraction encompassing the types

   DirectionOnGrid
   DirectionInDomain
   DirectionAtXY 
   Domain 
   NegateY
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
