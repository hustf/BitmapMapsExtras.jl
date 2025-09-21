 # This builds on domain, is a requisite for
 # streamlines.jl,  and define types for
 # 'non-discretization' of matrices, making them
 # work like continuous functions.
 #
 # Defines types:

 # DirectionAtXY (top level)
 # DirectionOnGrid
 # DirectionInDomain

struct DirectionOnGrid{F, T}
    Ω::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}# Square neighborhood
    fdir!::F # e.g 𝐧ₚ!(v, M)
    z::Matrix{T} # Image-sized
    lastvalue::MVector{2, Float64}
end
# Constructor
function DirectionOnGrid(fdir!, z)
    Ω = CartesianIndices((-2:2, -2:2))
    lastvalue = zero(MVector{2, Float64})
    DirectionOnGrid(Ω, fdir!, z, lastvalue)
end
# Callable, returns a 2d vector
function (dog::DirectionOnGrid)(i::Int64, j::Int64)
    rngi = dog.Ω.indices[1] .+ i
    rngj = dog.Ω.indices[2] .+ j
    vi = view(dog.z, rngi, rngj)
    dog.fdir!(dog.lastvalue, vi)
    dog.lastvalue
end
@define_show_with_fieldnames DirectionOnGrid


struct DirectionInDomain{F, T}
    dog::DirectionOnGrid{F, T}
    lastvalue::MVector{2, Float64}
end
# Constructor
function DirectionInDomain(fdir!, z::AbstractMatrix{<:AbstractFloat})
    dog = DirectionOnGrid(fdir!, z)
    lastvalue = zero(MVector{2, Float64})
    DirectionInDomain(dog, lastvalue)
end
# Callable, returns a 2d vector. Function of in-domain float coordinates. 'negy' is 0 at bottom of image.
function (did::DirectionInDomain)(x, negy)
    # Identify (up to four) points on grid
    j1, j2 =  Int(floor(x)), Int(ceil(x))
    i1, i2 =  Int(floor(negy)), Int(ceil(negy))
    # Shorthand 
    dog = did.dog
    lastvalue = did.lastvalue
    #
    if i1 != i2 && j1 != j2
        # Both pairs are different
        lastvalue .=  cornerweight(1, 1, x - j1, negy - i1) .* dog(i1, j1)
        lastvalue .+= cornerweight(2, 1, x - j1, negy - i1) .* dog(i2, j1)
        lastvalue .+= cornerweight(1, 2, x - j1, negy - i1) .* dog(i1, j2)
        lastvalue .+= cornerweight(2, 2, x - j1, negy - i1) .* dog(i2, j2)
    elseif i1 == i2 && j1 != j2
        # i1 equals i2, but j1 differs from j2
        lastvalue .=  sideweight(1, x - j1) .* dog(i1, j1)
        lastvalue .+= sideweight(2, x - j1) .* dog(i1, j2)
    elseif j1 == j2 && i1 != i2
        lastvalue .=  sideweight(1, negy - i1) .* dog(i1, j1)
        lastvalue .+= sideweight(2, negy - i1) .* dog(i2, j1)
        # j1 equals j2 
    else
        # i1 equals i2 and j1 equals j2
         lastvalue .= dog(i1, j1)
    end
    lastvalue
end
@define_show_with_fieldnames DirectionInDomain



abstract type DirectionFunctor end

"""
    struct DirectionAtXY <: DirectionFunctor
        did::DirectionInDomain
        negy::NegateY
        d::Domain
        v::MVector{2, Float64}
    end

The object type we hand to the differential equation solver. 

Callable: DirectionAtXY(fdir!, z)(x, y) -> 2d directional vector

This is a top-level abstraction encompassing the types

   DirectionOnGrid
   DirectionInDomain
   Domain 
   NegateY
"""
struct DirectionAtXY{F, T} <: DirectionFunctor
    did::DirectionInDomain{F, T}
    negy::NegateY
    d::Domain
    v::MVector{2, Float64}
end
# Constructor
function DirectionAtXY(fdir!, z)
    did = DirectionInDomain(fdir!, z)
    negy = NegateY(z)
    R = CartesianIndices(z)
    Ω = CartesianIndices((-2:2, -2:2))
    d = Domain(R, Ω)
    v = MVector{2, Float64}([0, 0])
    DirectionAtXY(did, negy, d, v)
end
# Callable. y is positive up on screen. Returns [0.0,0.0]
# when out-of domain (the callback set doesn't immediately
# stop calculation).
function (daxy::DirectionAtXY)(x, y)
    if daxy.d(x, y)
        daxy.v .= daxy.did(x, daxy.negy(y))
    else
        # This occurs fairly often, though isn't presented as part of the solution.
        #printstyled("::DirectionAtXY(x, y) == ($x, $y) ∉ $(daxy.d) \n", color =:176)
        daxy.v .= 0.0
        daxy.v
    end
    daxy.v
end
@define_show_with_fieldnames DirectionAtXY