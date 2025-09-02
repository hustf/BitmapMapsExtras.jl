 # This builds on direction_types, is a requisite for
 # streamlines.jl,  and define types for
 # 'non-discretization' of matrices, making them
 # work like continuous functions.
 #
 # Defines types:
 # BidirectionAtXY
 # BidirectionOnGrid
 # BidirectionInDomain

# Low level. See `allocations_curvature` for comments.
struct BidirectionOnGrid{F, T}
    fdir!::F # e.g 𝐊!. Input signature: (K, vα, vβ, vκ, P, M, vϕ)
    z::Matrix{T} # Image-sized
    Ω::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}
    vα::MVector{4, Float64}
    vβ::MVector{2, Float64}
    vκ::MVector{4, Float64}
    P::MMatrix{3, 3, Float64, 9}
    K::MMatrix{2, 2, Float64, 4}
    lpc::LinearCache
end
# Constructor
function BidirectionOnGrid(fdir!, z)
    _, Ω, _, P, K, vα, vκ, vβ, lpc = allocations_curvature(CartesianIndices(z))
    BidirectionOnGrid(fdir!, z, Ω, vα, vβ, vκ, P, K, lpc)
end
# Callable, returns a tensormap K for the pixel or cell
function (b::BidirectionOnGrid)(pt::CartesianIndex)
    b.fdir!(b.K, b.vα, b.vβ, b.vκ, b.P, view(b.z, b.Ω .+ pt), VΦ, b.lpc)::TENSORMAP
end
@define_show_with_fieldnames BidirectionOnGrid

# Intermediate level
struct BidirectionInDomain
    bdog::BidirectionOnGrid
    corners::SMatrix{2, 2, TENSORMAP, 4}
end
# Constructor
function BidirectionInDomain(fdir!, z::AbstractMatrix{<:Number})
    # Can generate a zero corner value (a 4x4), mutable
    # Establish four mutable corners in a static matrix 
    corners = SMatrix{2, 2, TENSORMAP, 4}([TENSORMAP(zeros(Float64, 2, 2)) for i in 1:2, j in 1:2] )
    # Store as BidirectionInDomain
    BidirectionInDomain(BidirectionOnGrid(fdir!, z), corners)
end
# Callable, returns K, a TENSORMAP for floating point coordinates
function (bid::BidirectionInDomain)(x, negy)
    # Identify (up to four) points on grid
    j1, j2 =  Int(floor(x)), Int(ceil(x))
    i1, i2 =  Int(floor(negy)), Int(ceil(negy))
    update_corners!(bid::BidirectionInDomain, i1, j1, i2, j2)
    interpolate_unit_square(bid.corners, x - j1, negy - i1)
end
@define_show_with_fieldnames BidirectionInDomain

"""
    struct BidirectionAtXY <: DirectionFunctor
        bid::BidirectionInDomain
        negy::NegateY
        d::Domain
        K::TENSORMAP
        primary::Ref{Bool}
        primary_start::Bool
    end

This carries no sense of direction. Calling an instance
with coordinates will return a 180° symmetric vector.
"""
struct BidirectionAtXY
    bid::BidirectionInDomain
    negy::NegateY
    d::Domain
    K::TENSORMAP
    primary::Ref{Bool}
    primary_start::Bool
end
# Constructor
function BidirectionAtXY(fdir!, z, primary::Bool)
    bid = BidirectionInDomain(fdir!, z)
    negy = NegateY(z)
    R = CartesianIndices(z)
    Ω = CartesianIndices((-2:2, -2:2))
    d = Domain(R, Ω)
    K = TENSORMAP(zeros(Float64, 2, 2))
    BidirectionAtXY(bid, negy, d, K, Ref{Bool}(primary), primary)
end
# Callable. The returned 2d vector is ambiguous, i.e. "180°-symmetric". 
function (baxy::BidirectionAtXY)(x::Float64, y::Float64)
    if baxy.d(x, y)
        baxy.K .= baxy.bid(x, baxy.negy(y))
        if baxy.primary[]
            return baxy.K[:, 1]
        else
            return baxy.K[:, 2]
        end
    else
        # (x, y) is out of domain. Solvers do this fairly often.
        baxy.K .= 0.0
        return baxy.K[:, 1]
    end
end
@define_show_with_fieldnames BidirectionAtXY




"""
    struct UnidirectionAtXY <: DirectionFunctor
        baxy::BidirectionAtXY
        du::MVector{2, Float64}
        flip::Ref{Bool}
        flip_start::Bool
    end

Top level for integrating to streamlines.
`flip` carries data on which of two directions we're
currently following. See 'reset!'.
"""
struct UnidirectionAtXY <: DirectionFunctor # DirectionFunctor
    baxy::BidirectionAtXY
    du::MVector{2, Float64}
    flip::Ref{Bool}
    flip_start::Bool
end
# Constructor
function UnidirectionAtXY(fdir!, z, primary::Bool, flip::Bool)
    du = MVector{2, Float64}([0.0, 0.0])
    baxy = BidirectionAtXY(fdir!, z, primary)
    UnidirectionAtXY(baxy, du, Ref{Bool}(flip), flip)
end
# Callable. The returned 2d vector is "360°-symmetric"
function (uxy::UnidirectionAtXY)(x::Float64, y::Float64)
    if uxy.flip[]
        uxy.du .=- uxy.baxy(x, y)
    else
        uxy.du .= uxy.baxy(x, y)
    end
    uxy.du
end
@define_show_with_fieldnames UnidirectionAtXY


# Resetting will not destroy the inherent z-data,
# just prepare for integrating from a new starting point.
function reset!(uxy::UnidirectionAtXY)
    uxy.flip[] = uxy.flip_start    
    uxy.baxy.primary[] = uxy.baxy.primary_start
    uxy.du .= 0.0
    uxy.baxy.K .= 0.0
    uxy
end
reset!(fxy::DirectionFunctor) = fxy