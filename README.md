# BitmapMapsExtra

High dependency, dubious gain: Experiments on top of the BitmapMaps pipeline.

This package is intended for modifying Bitmaps graphics.

## Drawing

On its own, BitmapMaps can draw single-color lines, triangles etc. top of a bitmap already. Text is provided in a vector graphics - svg - layer.

BitmapMapsExtras can:

- Draw colorful and transparent glyphs
- Pack glyphs, respecting the size of each
- Make simple grids of glyphs, or on directly specified locations.
- Draw streamline plots. Most typically, uphill or downhill streamlines, but also
  lines of curvature. Such plots are actually ensembles of (very fast) differential equation solutions.
- Paint regions. Typically convex, concave, convex-concave or plane regions. 

## Types of glyphs 

- Vector glyph 
  
  Similar to an arrow-style glyph, shows a positive value, and a direction. The default glyph is calculated from surface normals.

- Tangent base glyph

  Shows the surface-aligned tangent basis, i.e. a 3d coordinate system which is used to calculate curvature. 

- Bidirectional vector glyph
  
  Shows a signed value and a line of action, but no single direction along that line. The property shown is 180° symmetrical, such as tensile force in a rope. The default glyph is calculated from one of two principal curvatures of a surface. 
  
- Principal directions glyph

  This is a pair of bidirectional vector glyphs. When viewed directly along a surface normal, they are perpendicular. The default shows principal surface curvatures - the signed values of curvature in the direction of maximum and minimum curvature. 

## Differential geometry functions

These functions have discrete domains. They are fast estimators for properties at cells of a gridded surface, like a grey-scale image or elevation map.

- Normal vector projection, 𝐧ₚ

Where BitmapMaps estimates based on a 3x3 grid, this uses a 5x5 grid around each pixel.

`𝐧ₚ!` is the default function for making a `DirectionOnGrid` functor or object. `DirectionOnGrid` is the default `AbstractIJFunctor` for `GSVector` glyph specifications.


- Tangent basis, `tangent_basis`

  Contrary to textbook practice, we choose a tangent x-axis in the elevation surface's yz-plane. This samples a 5x5 window.

- Principal curvatures components, 𝐊

  Samples a 5x5 window around each pixel or cell. This estimator is optimized against known curvature functions for spheroids, cylinders, etc. Still, it is inaccurate for steep terrain (> 60°).

## Streamlines

The familiar example is tracing a curve directly downhill from a point. Doing this accurately requires that we make continuous domains. We do this internally with `DirectionAtXY` or `BidirectionAtXY` types.

`DirectionAtXY` enables integration along lines of descent or ascent. Visually, that
traces out a streamline. If we are only interested in the streamline shape for graphical purposes, we can integrate `𝐧ₚᵤ` (the pure direction of descent) instead of `𝐧ₚ`. The result of integration along such a path is simple the elevation difference between start and finish. 

`BidirectionAtXY` enables streamline integration along lines of curvature. Since curvature is a matrix rather than a vector, we must choose between following lines of maximal or minimal curvature. The result of integration along such a line is change in steepness angle. However, since we're dealing with 2.5D height maps, this restricts the range to <-π, π>. By intuition, any streamline ought to be cyclic, or end in a plane like a lake. 

For both types of streamline, they are more interesting in groups. For example, downhill streamlines collect in valleys, uphill streamlines collect in ridges and local summits. 

A collection of streamlines starting points is an ensemble. An ensemble might be a grid, uniformly random, or based on e.g. divergence (as in fluid flow) or previously drawn [hachures](https://andywoodruff.com/blog/hachures-and-sketchy-relief-maps/) in an iterative approach. 


We're using [OrdinaryDiffEq](https://docs.sciml.ai/OrdinaryDiffEq/stable/) to make these streamlines. Since our functors, or objects, are memory allocation free, we can find streamlines very effectively, even without advanced parallelism.


## Current progress

Version 0.0.14
Fully functional public API for all glyphs:
- `plot_glyphs!`
- `pack_glyphs!` 
- `indices_on_grid`

Glyph packing passes all tests.

- `plot_streamlines!` and `plot_bidirec_streamlines!` might also benefit from some similar, e.g. AbstractStreamlineSpec. Tests are currently broken, as we are dropping 
the temporary interface `grid_fcall` in 'common.jl' and will be working on a new streamlines interface.

Version 0.0.13
Added more crash tests for 'pack_glyphs'. Mostly working, but static dispatch on tensor glyph directions would be useful. Currently, packing 2d vectors will trigger an error. We plan to revise the GSTensor type.

Version 0.0.12
Add 'pack_glyphs' and restructure plotting function hierarchy. So far working for vectors, but more restructuring needed.

Version 0.0.11 Add missed file

Version 0.0.10 drops dependency ImageTransformations. Custom sub-pixel interpolations uses less memory allocation.

Version 0.0.9 outlines (currently in `test/t_calculate_and_draw_glyphs.jl`) a method to densely pack tensormap glyphs.

Version 0.0.8 adds a secondary color to the tensor basis glyph spec. Now, a two-color glyph can be drawn in one operation and faster.

Version 0.0.7 pre-allocates the linear solver for curvature directions, and makes other minor speed improvements.

We are revising and improving on a user interface, i.e. the high-level functions. The somewhat unsatisfactory interfaces reside in `test/common.jl`. Some test files are not updated to the latest version of the interface. One idea for a high-level interface is to pass types like glyph specification into more general functions. `grid_fcall_with_background` etc. is not a nice way to do this. 

`tangent_basis` has some minor allocations we should get rid of, but it may add more fields to some types.

Glyph stacking is currently implemented, but in an inefficient manner. A new algorithm for stacking glyps is planned, but a requirement for that algorithm is to know
the user's glyph specification in the stacking algorithm. Next up is to implement `glyphspec_types.jl` in `test/t_calculate_and_draw_glyphs`.

