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
  
- Draw tangent basis planes
  
  For checking that the basis planes are reasonable.
  
- Draw traces
  
  For streamplots of 2d vector-fields. Each streamline can add opacity, in order to visualize divergence


## Differential geometry functions

- Extract tangent basis from height-maps

  Here, the tangent x-axis lies in the yz-plane. Samples a 5x5 window in a non-allocating way. BitmapMaps already finds the surface normal vector, but for a smaller window. This is more robust, but will also hide a little detail.

- Find principal curvatures from height-maps in a non-allocating way.
  
  The principal curvatures are defined in the tangent plane, so in steep terrain, our estimatition is inaccurate. We optimized the accuracy by comparing with analytically known curvatures for spheroids, cylinders, etc. Samples a 5x5 window around each pixel or cell.

## Integrate to find streamlines

  - Find traces starting at an ensemble of points

  An ensemble might be a grid, uniformly random, or based on e.g. divergence (as in fluid flow) or previously drawn [hachures](https://andywoodruff.com/blog/hachures-and-sketchy-relief-maps/)  in an iterative approach. 

  User supplies the vector-field function for the differential equation, typically from differential geometry. We solve the equation in non-discrete coordinates, which requires interpolation between the neighbouring points. We do this 'lazily', i.e. we only interpolate where needed and try to avoid memory allocations.

  We're actually just wrapping some functionality from [OrdinaryDiffEq](https://docs.sciml.ai/OrdinaryDiffEq/stable/). The solutions are stored very effectively, so we finish finding the solutions before we start drawing traces.
  
