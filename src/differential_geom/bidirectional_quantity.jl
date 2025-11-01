"""
    components_matrix!(K, κ1, κ2, vβ)

Encode as κ1, κ2 and vβ as K.

K is a 2x2 matrix. It contains a second-order tensor's components in a screen-aligned basis.

Convert the signed scalar (κ) + angle (β) representation into a column vector per component. 
This is believed to require less trigonometric computations in further processing.

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
    K
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




