 # This defines types for
 # specifying glyph properties.
 # NOT currently implemented 

abstract type GlyphSpec end 

struct GSTensor <: GlyphSpec
    multip::Float64
    ming::Float64
    maxg::Float64
    dashsize::Float32
    strength::Float32 # Graphical spray strength
    color::RGB{N0f8}
    directions::UnitRange{Int64}
end
function GSTensor(;multip = 50, ming = -50, maxg = 50,
    dashsize = Float32(maxg / 10), strength = 0.7f0, color = COLOR_CURVGLYPH,
    directions = 1:2)
    #
    GSTensor(float(multip), float(ming), float(maxg), Float32(dashsize), Float32(strength),
        color, to_range(directions))
end
to_range(x::Int64) = x:x
to_range(x::UnitRange{Int64}) = x 
@define_show_with_fieldnames GSTensor
function is_in_limits(gs::GSTensor, K::TENSORMAP)
    if 1 ∈ gs.directions && 2 ∈ gs.directions
        # Checks a glyph for max and min principal directions. 
        K1 = @view K[:, 1]
        K2 = @view K[:, 2]
        limit1 = is_bidirec_vect_positive(K1) ? gs.maxg : abs(gs.ming)
        limit2 = is_bidirec_vect_positive(K2) ? gs.maxg : abs(gs.ming)
        norm(K1) * gs.multip ≤ limit1 || return false
        norm(K2) * gs.multip ≤ limit2 || return false
        return true
    elseif 2 ∈ gs.directions || 1 ∈ gs.directions
        # Checks a glyph for one principal direction. 
        v = view(K, :, first(gs.directions))
        limit = is_bidirec_vect_positive(v) ? gs.maxg : abs(gs.ming)
        return norm(v) * gs.multip ≤ limit
    else
        throw(ErrorException("Bad: $(gs.directions)"))
    end
end


struct GSVector <: GlyphSpec
    multip::Float64
    ming::Float64
    maxg::Float64
    dashsize::Float32
    strength::Float32 # Graphical spray strength
    color::RGB{N0f8}
end
function GSVector(; multip = 50, ming = -50, maxg = 50,
    dashsize = Float32(maxg / 10), strength = 0.7f0, color = COLOR_CURVGLYPH)
    #
    GSVector(float(multip), float(ming), float(maxg), Float32(dashsize), Float32(strength), color)
end
@define_show_with_fieldnames GSVector
is_in_limits(gs::GSVector, v) = gs.ming ≤ norm(v) ≤ gs.maxg


struct GSTangentBasis <: GlyphSpec
    halfsize::Int64
end
GSTangentBasis(;halfsize = 30) = GSTangentBasis(halfsize)
@define_show_with_fieldnames GSTangentBasis
