
"""
# Explanation

## Standard Basis

The standard orthonormal right-handed basis vectors are:

- 𝐢 = [1, 0, 0]  
- 𝐣 = [0, 1, 0]  
- 𝐤 = [0, 0, 1]

## Tangent Basis

At each point on a surface, we define a tangential plane.
The tangent basis vectors are defined in terms of the standard basis:

- 𝐞₁ = a⋅𝐢 + b⋅𝐣 + c⋅𝐤  
- 𝐞₂ = d⋅𝐢 + e⋅𝐣 + f⋅𝐤  
- 𝐞₃ = g⋅𝐢 + h⋅𝐣 + i⋅𝐤

These vectors form a right-handed orthonormal basis where:

‖𝐞₁‖ = ‖𝐞₂‖ = ‖𝐞₃‖ = ‖𝐢‖ = ‖𝐣‖ = ‖𝐤‖ = 1

Both bases share the same origin.

## Vector Representation

Let 𝐩 ∈ ℝ³ be a geometric vector with two coordinate representations:

- In the standard basis:  𝐩 = [x, y, z] = x⋅𝐢 + y⋅𝐣 + z⋅𝐤  
- In the tangent basis:   𝐮 = [u, v, w] = u⋅𝐞₁ + v⋅𝐞₂ + w⋅𝐞₃

Here, 𝐩 and 𝐮 represent the same vector. The rightmost expressions are equal, 
but the bracketed expressions by themselves leave out which base they are defined in.

## Constraints

- 𝐞₁ lies in the xz-plane ⇒ b = 0  
- 𝐞₁ points in the positive x direction ⇒ a ≥ 0

## Transformation from tangent to standard basis

Let P be the 3×3 matrix whose columns are the tangent basis vectors expressed in the standard basis:
         
    P = [𝐞₁ 𝐞₂ 𝐞₃] = 
        [a  d  g;
         0  e  h;
         c  f  i]

Then for any vector 𝐩 ∈ ℝ³ with coordinates 𝐮 = [u, v, w] in the tangent basis, 
the corresponding standard basis coordinates are:

    𝐩 = P * 𝐮 =
        [a   d   g;
         0   e   h;
         c   f   i] * [u;
                      v;
                      w] =
         u𝐞₁ + v𝐞₂ + w𝐞₃


Given 𝐮 = [u, v, w], the standard basis coordinates 𝐩 = [x, y, z] are:

    x = a·u + d·v + g·w
    y = e·v + h·w
    z = c·u + f·v + i·w

## Transformation from standard to tangent basis

As defined above, P is a real, square, orthogonal matrix. Therefore, its inverse equals its transpose:

    P⁻¹ = Pᵀ = 
           [a   0   c;
            d   e   f;
            g   h   i]

This simplifies conversion from standard basis coordinates 𝐩 to tangent basis coordinates 𝐮:

    𝐮 = P⁻¹ * 𝐩 = Pᵀ * 𝐩

Each component of 𝐮 represents how much of 𝐞₁, 𝐞₂ or 𝐞₃ is present in 𝐩. In other words, 𝐮 gives the unique coefficients such that:

    𝐩 = u·𝐞₁ + v·𝐞₂ + w·𝐞₃

Given 𝐩 = [x, y, z], the tangent basis coordinates 𝐮 = [u, v, w] are:

    u = a·x +     c·z
    v = d·x + e·y + f·z
    w = g·x + h·y + i·z

"""
foo



"""
# Reversible projection of an angle to the tangent plane

Let 𝐩 be a vector on the unit circle in the xy-plane, centered at the origin, rotating around the z-axis (𝐤). 
It represents an angle α measured from 𝐢 counterclockwise:

    𝐩 = [cos(α), sin(α), 0] = cos(α)⋅𝐢 + sin(α)⋅𝐣 

Let 𝐪 be a vertical vector (parallel with 𝐤), chosen such that the sum 𝐫 = 𝐩 + 𝐪 lies in the tangent 
plane (𝐞₁, 𝐞₂):

    𝐪 = [0, 0, z], ∣ 𝐫 = 𝐩 + 𝐪 ∈ span(𝐞₁, 𝐞₂)

Expanding:

    𝐫 = cos(α)·𝐢 + sin(α)·𝐣 + z·𝐤


Substituting 𝐢, 𝐣, and 𝐤 with 𝐞₁, 𝐞₂, and 𝐞₃:

    𝐫 = cos(α)(a·𝐞₁ + d·𝐞₂ + g·𝐞₃)
        + sin(α)(e·𝐞₂ + h·𝐞₃)
        + z·(g·𝐞₁ + h·𝐞₂ + i·𝐞₃)

Collecting terms:

    𝐫 = (a·cos(α) + z·g)⋅𝐞₁  +  
        (d·cos(α) + e·sin(α) + z·h)⋅𝐞₂  +
        (g·cos(α) + h·sin(α) + z·i)⋅𝐞₃

Because 𝐫 lies in the uv plane, we solve for z, requiring the 𝐞₃ component to vanish:

    (g·cos(α) + h·sin(α) + z·i) = 0
    ⇒ z = -(g·cos(α) + h·sin(α)) / i
    where  i ≠ 0.

Recall that 𝐞₃ = g·𝐢 + h·𝐣 + i·𝐤 — that is, i is the z-component of the tangent plane normal 
in standard coordinates. If i = 0, then 𝐞₃ lies entirely in the xy-plane, and the tangent plane 
is vertical. In this case, projecting along the z-axis onto the tangent plane is 
undefined. This is consistent with expectations: for a surface defined as a heightfield 
(a "2.5D" surface), the tangent plane is never vertical, so i ≠ 0 is a reasonable and 
typical assumption.

Note: Small values of i (close to zero) may amplify numerical errors in floating-point implementations, affecting the accuracy of z and φ.

With this z, the 𝐞₃ component of 𝐫 vanishes, so 𝐫 lies in the tangent plane.

    z = -(g·cos(α) + h·sin(α)) / i
    ∧
    𝐫 = (a·cos(α) + z·g)⋅𝐞₁ + (d·cos(α) + e·sin(α) + z·h)⋅𝐞₂ 

Then we define ϕ to be the angle from 𝐞₁ to the direction of 𝐫 within the tangent plane. It is computed by:
    
    ϕ = atan(d·cos(α) + e·sin(α) + z·h, a·cos(α) + z·g)

This represents the projected angle of α into the tangent plane along the z-axis. 
This is a coordinate-driven construction, not a universal projection of angle
between planes. It differs from decomposing rotational vectors such as torque and spin.

# Projection of an angle from the tangent plane to the xy plane

This describes the inverse of the earlier projection. Here, we project angle ϕ
along the same axis, z || 𝐤, with α as output. The operation is reversible for 
i ≠ 0, but numerical precision may affect results when i is small.

Let 𝐮 be a vector on the unit circle in the tangent plane, centered at the origin. 
It represents an angle ϕ measured from 𝐞₁ counterclockwise:

    𝐮 = [u, v, w] = [cos(ϕ), sin(ϕ), 0]

We already [deduced](##_transformation_from_tangent_to_standard_basis) how to express 
the same vector in the standard basis:
    𝐩 = [a·u + d·v + g·w, e·v + h·w, c·u + f·v + i·w]

Substituting components of 𝐮:

    𝐩 = [a·cos(ϕ) + d·sin(ϕ), 
         e·sin(ϕ), 
         c·cos(ϕ) + f·sin(ϕ)]

We now project 𝐩 along the z axis onto the xy-plane, which amounts to ignoring the third component.
The projected point lies on a ray from origin, and we define α as the angle this ray makes 
with the x-axis:

    α = atan(e·sin(ϕ), a·cos(ϕ) + d·sin(ϕ))

"""
faa