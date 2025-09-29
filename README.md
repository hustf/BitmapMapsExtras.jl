# BitmapMapsExtra

High dependency, dubious gain: Experiments on top of the BitmapMaps pipeline.

## Drawing

BitmapMaps can draw single-color lines, triangles, squares on top of a bitmap already (text is provided
in a vector graphics - svg - layer.

Here we add other stuff. So far:

- Draw vector glyps
  
  For gradients and other vector fields

- Draw bidirectional vector glyphs
  
  For curvature, which is 180° symmetrical
  
- Draw principal directions glyphs
  
  For principal curvatures. In the tangent plane, these are orthornormal, but not in the drawing plane.

- Glyph distribution

  Grids of glyps, but also glyph-size based placement. The latter is work in progress.
  
- Draw tangent basis planes
  
  For checking that the basis planes are reasonable.
  
- Paint surface type

  For visualizing valleys and similar curvature-based features. May be dropped.

- Streamlines

  Traces for gradients, but also for lines of curvature. The latter is work in progress.



## Differential geometry functions

- Extract tangent basis from height-maps

  Here, the tangent x-axis lies in the yz-plane. Samples a 5x5 window in a non-allocating way. BitmapMaps already finds the surface normal vector, but for a smaller window. This is more robust, but will also hide a little detail.

- Find principal curvatures from height-maps in a non-allocating way.
  
  The principal curvatures are defined in the tangent plane, so in steep terrain, our estimatition is inaccurate. We optimized the accuracy by comparing with analytically known curvatures for spheroids, cylinders, etc. Samples a 5x5 window around each pixel or cell.

## Integrate to find streamlines

  - Find traces starting at an ensemble of points

  An ensemble might be a grid, uniformly random, or based on e.g. divergence (as in fluid flow) or previously drawn [hachures](https://andywoodruff.com/blog/hachures-and-sketchy-relief-maps/)  in an iterative approach. Work in progress.

  User supplies the vector-field function for the differential equation, typically from differential geometry. We solve the equation in non-discrete coordinates, which requires interpolation between the neighbouring points. We do this 'lazily', i.e. we only interpolate where needed and try to avoid memory allocations.

  We're using [OrdinaryDiffEq](https://docs.sciml.ai/OrdinaryDiffEq/stable/). The solutions are stored very effectively, so we finish finding the solutions before we start drawing traces.


## Current progress

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

