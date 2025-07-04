module TestMatrices
export r
export z_cylinder, z_sphere, z_ellipsoid, z_paraboloid
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
on the cylinder. (0,0) outside the circle, (∞, 0) 
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
    z_paraboloid(; a = 1.0, b = 0.5a)
    
    z(x, y) = (x^2)/(2a) + (y^2)/(2b)

Elliptic paraboloid:  a > 0, b > 0
Hyperbolic paraboloid: Switch the sign of a or b. Saddle at the orgin.

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