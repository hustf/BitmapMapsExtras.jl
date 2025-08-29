 # This defines types for
 # specifying glyph properties.
 # NOT currently implemented 

abstract type GlyphSpec end 
struct GlyphSpecTensor <: GlyphSpec
    multip::Int64
    maxg::Int64
    ming::Int64
    dashsize::Int64
    strength::Float32
    directions
end
function GlyphSpecTensor(;multip = 50, maxg = 50, ming = -50, 
    dashsize = maxg ÷ 10, strength = 0.7f0, directions = 1:2)
    #
    GlyphSpecTensor(multip, maxg, ming, dashsize, Float32(strength),
        directions)
end


struct GlyphSpecVector <: GlyphSpec
    multip::Int64
    maxg::Int64
    ming::Int64
    dashsize::Int64
    strength::Float32
    color
end
function GlyphSpecVector(; multip = 50, maxg = 50, ming = -50, 
    dashsize = maxg ÷ 10, strength = 0.7f0, color = COLOR_CURVGLYPH)
    #
    GlyphSpecVector(multip, maxg, ming, dashsize, Float32(strength), color)
end


struct GlyphSpecTangentBasis <: GlyphSpec
    halfsize
end
GlyphSpecTangentBasis(;halfsize = 30) = GlyphSpecTangentBasis(halfsize)


