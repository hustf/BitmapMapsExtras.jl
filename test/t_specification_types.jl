using Test
using BitmapMapsExtras
using BitmapMapsExtras: GSTensor, GSVector, GSTangentBasis, Stroke
using BitmapMapsExtras: AbstractGlyphSpec, AbstractCurveSpec


@test GSTensor() isa AbstractGlyphSpec
@test GSVector() isa AbstractGlyphSpec
@test GSTangentBasis() isa AbstractGlyphSpec
@test Stroke() isa AbstractCurveSpec