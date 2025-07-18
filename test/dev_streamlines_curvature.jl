# These functions are not quite successful at present.
# The direction sometimes flip unexpectedly, such as when
# crossing the vertical line at centre of z_paraboloid(a = 0.5r, b = -r)
#
# Curvature streamlines may not be the most useful tool even when
# working fully, but the issue ought to be resolved before adding
# e.g. ridge-following streamlines.


using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras.TestMatrices: I0
using BitmapMapsExtras: allocations_curvature, principal_curvature_components!, VΦ, tangent_basis!


c11!() = func_direction_from_principal_curvature(;max_principal = true, initial_direction = true)
c12!() = func_direction_from_principal_curvature(;max_principal = true, initial_direction = false)
c21!() = func_direction_from_principal_curvature(;max_principal = false, initial_direction = true)
c22!() = func_direction_from_principal_curvature(;max_principal = false, initial_direction = false)


function func_direction_from_principal_curvature(;max_principal::Bool = true, initial_direction::Bool = true)
    # Consider this.
    # Nah, we could make the choice after constructing DirectionAtXY.
    # Or, define a type BiDirectionAtXY if fdir outputs a 2x2 matrix. 
    # The last choice would allow us to dispatch to methods of RHS, parameter p.
    _, _, _, P, K, vα, vκ, vβ, _ = allocations_curvature(CartesianIndices((1:1, 1:1)), [])
    # In-place, mutating function. v mutates to hold result, a 2d vector 
    # representing one component of curvature. M is a 5x5 matrix, presumably a sliding window.
    # v represents the specified component of curvature at M's centre pixel.
    (v, M) -> let P = P, K = K, vα = vα, vκ = vκ, vβ = vβ, vϕ = VΦ
        # Update P in-place, v is temporarily used as storarge 
        tangent_basis!(P, v, M)
        # Find bidirectional quantity. For simplicity,
        # we do calculate both although only one will be used.
        # Update K 
        principal_curvature_components!(K, vα, vβ, vκ, P, M, vϕ)
        # TEMP @show K norm(K[:,1]) norm(K[:,2])
        # Pick max or min principal
        column = max_principal ? 1 : 2
        # Form output directional 2d vector.
        v .= K[:, column]
        if ! is_bidirec_vect_positive(v)
            # TEMP println("positive")
            v .=- v
        else 
            # TEMP println("negative")
        end
        if initial_direction
            v .=- v
        end
        # So far, the norm of v is 'actual' curvature. 
        # But the variation in the norm is large for practical resolution of streamlines.
        # Hence, we 'normalize' the length of vectors by imagining a 3d-vector of unit length. 
        # Thus, output is better suited to effective intergration where only the direction matters.
        # mag = sqrt(1 + v[1]^2 + v[2]^2)
        #v[1] = -dz_x / mag
        #v[2] = -dz_y / mag
        #mag = sqrt(v[1]^2 + v[2]^2)
        #v ./= mag
        v
    end
end