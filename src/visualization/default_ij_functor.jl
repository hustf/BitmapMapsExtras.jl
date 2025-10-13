# This creates the default choice of AbstractIJFunctor,
# given the type of glyph and elevation map z.

default_ij_functor(z, gs::AbstractGlyphSpec) = throw(ArgumentError("No default functor for $gs"))
default_ij_functor(z, gs::GSTensor) = BidirectionOnGrid(𝐊!, z)
default_ij_functor(z, gs::GSVector) = DirectionOnGrid(𝐧ₚ!, z)