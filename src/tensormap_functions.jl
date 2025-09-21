# File contains utilty functions for
# the 4x4 mutable matrix ('K') we use as a tensor map
# for the 'world basis'.  

"""
    interpolate_unit_square!(value::TENSORMAP, corners::SMatrix{2, 2, TENSORMAP}, xus::AbstractFloat, yus::AbstractFloat)

Implemented here since we couldn't make Interpolations.jl's return @inferred.
"""
function interpolate_unit_square!(value::TENSORMAP, corners::SMatrix{2, 2, TENSORMAP}, xus::AbstractFloat, yus::AbstractFloat)
    x = clamp(xus, 0.0, 1.0)
    y = clamp(yus, 0.0, 1.0)
    if xus !== x
        throw(" ooo: x")
    end
    if yus !== y
        throw(" ooo: y")
    end    
    f11 = corners[1, 1] # (x, y) = (0, 0) 
    f12 = corners[1, 2] # (x, y) = (0, 1)
    f21 = corners[2, 1] # (x, y) = (1, 0)
    f22 = corners[2, 2] # (x, y) = (1, 1)
    value .= (1 - x) * (1 - y) .* f11 .+ 
                   (1 - x) * y .* f21 .+
                   x * (1 - y) .* f12 .+ 
                         x * y .* f22
    value
end

"""
    cornerweight(i, j, x::Float64, y::Float64)
    --> Float64

For interpolation in a unit square
"""
function cornerweight(i, j, x::Float64, y::Float64)
    i == 1 && j == 1 && return (1 - x) * (1 - y)
    i == 2 && j == 1 && return (1 - x) * y
    i == 1 && j == 2 && return x * (1 - y)
    i == 2 && j == 2 && return x * y
    throw(ErrorException("$i $j"))
end
"""
    sideweight(i, j, x::Float64, y::Float64)
    --> Float64

For interpolation along a line 0..1.
"""
function sideweight(j, x::Float64)
    j == 1 && return (1 - x)
    j == 2 && return x
    throw(ErrorException("$j $x"))
end
"""
    update_corners!(bid::BidirectionInDomain, i1, j1, i2, j2)

Calculate and store tensormap K for closest on-grid coordinates.
Purpose is linear interpolation.

Frimary and secondary directions are not allowed to swap 
column order between corners. This is to ensure
that such a swap is detectable, even with interpolation.
"""
function update_corners!(bid::BidirectionInDomain, i1, j1, i2, j2)
    # Alias
    mc = bid.corners
    bdog = bid.bdog
    # Re-use calculated K when indices are equal.
    # E.g. for x = 100.0, then j1 == j2
    if i1 != i2 && j1 != j2
        # Both pairs are different
        mc[1, 1] .= bdog(i1, j1)
        mc[2, 1] .= bdog(i2, j1)
        mc[1, 2] .= bdog(i1, j2)
        mc[2, 2] .= bdog(i2, j2)
    elseif i1 == i2 && j1 != j2
        # i1 equals i2, but j1 differs from j2
        mc[1, 1] .= bdog(i1, j1)
        mc[2, 1] .= mc[1, 1]
        mc[1, 2] .= bdog(i1, j2)
        mc[2, 2] .= mc[1, 2]
    elseif j1 == j2 && i1 != i2
        # j1 equals j2 
        mc[1, 1] .= bdog(i1, j1)
        mc[2, 1] .= bdog(i2, j1)
        mc[1, 2] .= mc[1, 1]
        mc[2, 2] .= mc[2, 1]
    else
        # i1 equals i2 and j1 equals j2
        mc[1, 1] .= bdog(i1, j1)
        mc[2, 1] .= mc[1, 1]
        mc[1, 2] .= mc[1, 1]
        mc[2, 2] .= mc[1, 1]
    end
    # 
    if is_close_to_perpendicular(mc[1,1][:, 1], mc[2,1][:, 1])
        swap_columns!(mc, 2, 1)
    end
    if is_close_to_perpendicular(mc[1,1][:, 1], mc[1,2][:, 1])
        swap_columns!(mc, 1, 2)
    end
    if is_close_to_perpendicular(mc[1,1][:, 1], mc[2,2][:, 1])
        swap_columns!(mc, 2, 2)
    end
    bid
end

function swap_columns!(mc, i, j)
    # This isn't called very often
    firstc = mc[i, j][:, 1]
    mc[i, j][:, 1] .= mc[i, j][:, 2]
    mc[i, j][:, 2] .= firstc
    mc
end

"""
    is_close_to_perpendicular(v1::T, v2::T) where T <: MVector{2, Float64}

78°–102° returns true, if v1 or v2 are zero, returns false.

This interprets the two first elements as a 2d vector. Ignores other elements.

TODO: Prettify this. Don't return false. @btime this to optimize. norm or like this? 
Combine code `dot_product_with....`
"""
@inline function is_close_to_perpendicular(v1::T, v2::T) where T <: MVector{2, Float64}
    d = v1[1] * v2[1] + v1[2] * v2[2]
    n1 = v1[1]^2 + v1[2]^2
    n2 = v2[1]^2 + v2[2]^2
    # acos(0.2) * 180/π = 78.46
    n1 == 0 || n2 == 0 ? false : abs(d) < 0.2 * sqrt(n1 * n2)
end

