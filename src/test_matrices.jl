"""
This module is included in module BitmapMapsExtras.

# Example

```
julia> using BitmapMapsExtras

julia> using BitmapMapsExtras.TestMatrices

julia> using BitmapMapsExtras.TestMatrices: I0

julia> varinfo(TestMatrices)
  name                                  size summary
  ––––––––––––––––––––––––––––––– –––––––––– –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
  TestMatrices                    34.326 KiB Module
  principal_curvatures_paraboloid    0 bytes principal_curvatures_paraboloid (generic function with 2 methods)
  r                                  8 bytes Int64
  z_cylinder                         0 bytes z_cylinder (generic function with 1 method)
  z_cylinder_offset                  0 bytes z_cylinder_offset (generic function with 1 method)
  z_ellipsoid                        0 bytes z_ellipsoid (generic function with 1 method)
  z_paraboloid                       0 bytes z_paraboloid (generic function with 1 method)
  z_sphere                           0 bytes z_sphere (generic function with 1 method)
  z_ridge_peak_valleys               0 bytes z_ridge_peak_valleys (generic function with 1 method)
```

"""
module TestMatrices
export r
export z_cylinder, z_cylinder_offset, z_sphere, z_ellipsoid, 
    z_paraboloid, z_ridge_peak_valleys
export principal_curvatures_saddle, principal_curvatures_paraboloid
##############################################
# Define (parameters for) test array functions
##############################################
const w_h = 999
const r = w_h ÷ 2
const R = CartesianIndices((w_h, w_h))
const I0 = CartesianIndex((r + 1, r + 1))

"""
    z_cylinder(α)

A cylinder with its axis in the xy-plane,
axis tilted α to the x-axis.
Principal curvature is (0, -r) everywhere
on the cylinder. (0,0) outside the cylinder, (∞, 0) 
on edges.  
"""
function z_cylinder(α)
    map(R) do I
        y, x = (I - I0).I
        x1 = cos(α) * x - sin(α) * y
        y1 = cos(α) * y + sin(α) * x
        sqrt(max(0, r^2 - y1^2))
    end;
end
"""
    z_cylinder_offset(α)

A cylinder with its axis in the xy-plane,
axis tilted α to the x-axis, centre at a corner.
Principal curvature is (0, -r) everywhere
on the cylinder. (0,0) outside the cylinder, (∞, 0) 
on edges.  
"""
function z_cylinder_offset(α)
    map(R) do I
        y, x = I.I
        x1 = cos(α) * x - sin(α) * y
        y1 = cos(α) * y + sin(α) * x
        sqrt(max(0, r^2 - y1^2))
    end;
end

"""
    z_sphere()

A sphere with its centre in the xy-plane.
Principal curvature is (-r, -r) everywhere
on the sphere. (0,0) outside the circle, (∞, 0) 
on edges.  
"""
function z_sphere()
    map(R) do I
        y, x = (I - I0).I
        sqrt(max(0, r^2 - x^2 - y^2))
    end
end

"""
    z_ellipsoid(; rel_half_axis = 0.5)

An ellipsoid with its center in the xy-plane, major axis along y, and minor axis along x.
The relative half-axis length `rel_half_axis` scales the y-axis (b = rel_half_axis * a).
Principal curvature varies across the surface, becoming (0,0) outside the ellipse and (∞, 0) on edges.

The analytical expressions for curvatures are hard to implement, so not a good 
choice for optimization.
"""
function z_ellipsoid(; rel_half_axis = 0.5)
    map(R) do I
        y, x = (I - I0).I
        a = r  # Major half-axis along y
        b = rel_half_axis * a  # Minor half-axis along x
        discriminant = r^2 - (y^2 / a^2 + x^2 / b^2) * r^2
        sqrt(max(0, discriminant))
    end
end


"""
    z_paraboloid(; a = ;a = -0.5r, b = 0.5r)
    
    z(x, y) = (x^2)/(2a) + (y^2)/(2b)

Elliptic paraboloid, convex everywhere:  a > 0, b > 0
Hyperbolic paraboloid: Default. a and b of opposite sign. Saddle at the orgin.

Principal curvatures at (x, y) = (0, 0):

- κ₁ = 1 / a      (along x-axis)
- κ₂ = 1 / b      (along y-axis)
"""
function z_paraboloid(; a = 1.0, b = 0.5a)
    map(R) do I
        y, x = (I - I0).I
        (x^2)/(2a) + (y^2)/(2b)
    end
end

function zigzag(x)
    sel = mod(x, 2π) 
    y = 2 * mod(x, π) / π - 1
    if sel ≈ 0 # Floating points can be tricky, hence special case
        1.0
    elseif sel ≈ π
        -1.0
    elseif sel < π
        - y
    else
        y
    end
end


function z_ridge_peak_valleys(;
    λ_arc = 0.3182291666666667,   # Wave length fraction of 2π inner of spiral
    λ_rad = 1.875,    # Wave length along radial in terms of half-screen width
    spiral = 0.5
    )
    map(R) do I
        negy, x = (I - I0).I # Centre origin
        y = - negy
        ρ = sqrt(float(x^2 + y^2))
        θ = mod2pi(atan(y, x))
        # Non-linear radial, normalized by half-screen-diagonal
        rnl = max(0.0,  ( ρ / √2r) - λ_rad / 6)
        # Tangential phas
        ωt = 2π * (θ * λ_arc  + spiral * rnl)
        # Phase, radial 
        ωr = ρ * 2π / (λ_rad * r)
        # Pick cos or zigzag 
        shapepicker = (y + r)/(2r) # 1 at top to 0 at bottom
        shapepicker * zigzag(ωt) * zigzag(ωr) + (1 - shapepicker) * cos(ωt) * cos(ωr)
    end
end



"""
    principal_curvatures_paraboloid(x, y; a = 1.0, b = 0.5a)
    principal_curvatures_paraboloid(pt::CartesianIndex; a = 1.0, b = 0.5a)

Returns the principal curvatures (κ₁, κ₂) of the elliptic paraboloid

    z(x,y) = (x²)/(2a) + (y²)/(2b)

where -1 <= x <= 1 
      -1 <= y <= 1

evaluated at point (x,y).
"""
function principal_curvatures_paraboloid(x, y; a = 1.0, b = 0.5a)
    fx = x / a
    fy = y / b
    W2 = 1 + fx^2 + fy^2
    W32 = W2^(3/2)
    H = ((1 + fy^2) * (1/a) + (1 + fx^2) * (1/b)) / (2 * W32)
    K = (1 / (a * b)) / W2^2
    sqrt_discriminant = sqrt(H^2 - K)
    κ₁ = H + sqrt_discriminant
    κ₂ = H - sqrt_discriminant
    sort([κ₁, κ₂], rev = true)
end
function principal_curvatures_paraboloid(pt::CartesianIndex; a = 1.0, b = 0.5a)
    y, x = (pt - I0).I
    principal_curvatures_paraboloid(x, y; a, b)
end

end