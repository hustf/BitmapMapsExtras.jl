"""
    components_matrix!(K, κ1, κ2, vβ)

Convert the signed scalar (κ) + angle (β) representation into a column vector per component. This is
believed to require less trigonometric computations in further processing.

Note that due to symmetry, we only need 180° (i.e. π) to represent the direction. But we will
give components in all four quadrants. The additional info is used to determine if the 
curvature is negative or positive. 
    
This is our choice of convention:

| κ11 & κ12| κ21 & κ22| Description                    |
|-----|----|---------------------------------|
| +  | +  | First quadrant, positive, outward arrows  |
| -  | +  | Second quadrant, positive, outward arrows |
| -  | -  | First quadrant, negative, inward arrows   |
| +  | -  | Second quadrant, negative, inward arrows  |

"""
function components_matrix!(K, κ1, κ2, vβ)
    components_vector!(view(K, :, 1), κ1, vβ[1])
    components_vector!(view(K, :, 2), κ2, vβ[2])
end

"""
    components_vector!(v, κ, β)

See `components_matrix!`
"""
function components_vector!(v, κ, β)
    β⁺ = mod2pi(β)
    # Sawtooth 0-π, will map π to π unlike mod(β, π)
    if κ >= 0
        # First or second quadrant
        x = β⁺ > π ? β⁺ - π : β⁺
    else
        # Third or fourth quadrant
        x = β⁺ > π ? β⁺ : β⁺ + π
    end
    v[1] = abs(κ) * cos(x)
    v[2] = abs(κ) * sin(x)
end




"""
    is_bidirec_vect_positive(v::AbstractVector) --> Bool

v is a 2d vector representing a bi-directional quantity, 
i.e. one that is is π-symmetric or 180° symmetric, like uniaxial stress
or curvature.

If v = [x, y] is in the second or fourth quadrant: 
  - its sign is negative.
  - the quantity magnitude is √(x^2 + y^2)
  - it is directed along α = atan(y, x) \n 
    and also along α + π 

# Example
```
julia> is_bidirec_vect_positive([1.0, 0.0])
true

julia> is_bidirec_vect_positive([1.0, -0.0])
false

julia> is_bidirec_vect_positive([1.0, -1.0])
false

julia> is_bidirec_vect_positive([-1.0, -2.0])
false
```
"""
is_bidirec_vect_positive(v::AbstractVector) = !signbit(v[2])



