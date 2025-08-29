# Contains `principal_curvature_components` and it's non-allocating cousin, 
# as well as supporting functions.
# Also contains `allocations_curvature`, including glyph limits 
# for visualization. 
# Relies on `tangent_basis.jl` and `bidirectional_quantity.jl`.
# Relies on constants KERN´´, W, VΦ.

"""
    principal_curvature_components(z::Matrix, pt::CartesianIndex)

Note that there is an alias 𝐊 for brevity.
"""
function principal_curvature_components(z::Matrix, pt::CartesianIndex)
    Ri, Ω, v, P, K, vα, vκ, vβ, _ = allocations_curvature(CartesianIndices(z), [])
    if ! (pt ∈ Ri)
        # Too close to edge of z    
        K .= NaN
    else
        # Window 
        win = view(z, Ω .+ pt)
        # Update K etc. 
        principal_curvature_components!(K, vα, vβ, vκ, P, win, VΦ)
    end
    K
end


"""
    principal_curvature_components!(K, vα, vβ, vκ, P, M, vϕ)

M is the 5x5 input
vϕ are sample angles, normally the optimized values VΦ.

K, vα, vβ, vκ and P are pre-allocated.

Note that there is an alias 𝐊! for brevity.
"""
function principal_curvature_components!(K, vα, vβ, vκ, P, M, vϕ)
    @assert size(vα) == (4,)
    @assert size(vβ) == (2,)
    @assert size(vκ) == (4,)
    @assert size(P) == (3, 3)
    @assert size(M) == (5, 5)
    # Find P in-place. vβ is a temporary storage here. 
    tangent_basis!(P, vβ, M)
    # vϕ gives the sampling direction in the tangent plane
    # described by P. Sampling must be done in the screen plane,
    # since our data is an elevation matrix with indices as -y, x
    # 
    # Find 'screen angles' vα:
    angle_tangent_to_xy!(vα, vϕ, P)
    # Sample curvatures vκ
    sample_curvature_in_directions!(vκ, M, P, vα)
    # Calculations in the tangent plane
    κ1, κ2, ϕ1 = principal_curvature_and_direction(vκ, vϕ)
    # Principal angles vβ in the yx (screen) plane. Note that ϕ1 and ϕ2 are 
    # orthonormal in the tangent plane, but not generally in the yx-plane
    angle_tangent_to_xy!(vβ, [ϕ1, ϕ1 + π / 2], P)
    # Express curvatures κ and directions β as 4x4 matrix K.
    # K is a second-order tensor's components in a screen-aligned basis.
    components_matrix!(K, κ1, κ2, vβ)
end

"""
    sample_curvature_in_directions!(vκ, M, P, vα)
"""
function sample_curvature_in_directions!(vκ, M, P, vα)
    # Tangent plane normal's components:
    𝐧ₜₚ₁ = P[7]
    𝐧ₜₚ₂ = P[8]
    𝐧ₜₚ₃ = P[9]
    # This gives us an estimate for the derivative of the surface 
    # in the plane intersected by α
    # ζ´ = -(cos(α) * 𝐧ₜₚ₁ + sin(α) * 𝐧ₜₚ₂) / 𝐧ₜₚ₃
    # Mutates vκ according to vα
    map!(vκ, vα) do α
        # Shorthand
        c = cos(α)
        s = sin(α)
        #
        # Regarding direction, consider the comments on flipping to a 
        # right-handed situation below.
        αl = α - π / 2 
        # Define the transformation
        tfm = recenter(RotMatrix(αl), center(M))
        # Lazy rotated matrix.
        Mr = WarpedView(M, tfm; fillvalue = Flat())
        # Normal curvature and unit tangent vector in the 'matrix coordinate system',
        # where when α = 0 : 
        #           local axis lx points along increasing rows in the unrotated matrix (that is: down in an image).
        #           Sampling from increasing rows can be much faster than sampling from decreasing rows because of how
        #           memory access works.
        #       when α = π / 4: 
        #           local axis lx points along decreasing columns (that is: left in an image)
        # First derivative, taken from the tangent plane
        z´ = -(c * 𝐧ₜₚ₁ + s * 𝐧ₜₚ₂) / 𝐧ₜₚ₃
        # Second derivative
        z´´ = KERN´´ ⋅ Mr[W]
        # Extrinsic curvature in this direction (from Kreyszig AEM)
        κ =  z´´ / (1 + z´^2)^(3 / 2)
        if isnan(κ)
            @show z´ z´´ M 𝐧ₜₚ₁ 𝐧ₜₚ₂ 𝐧ₜₚ₃ P
            @show Mr[W] α αl
            throw(ErrorException("'Spurious' interpolation error "))
        end
        κ
    end
    vκ
end

"""
    angle_xy_to_tangent!(vϕ, vα, P)

This works with tangent bases with the condition that
the tangent plane's primary axis lies in the yz plane, and 
the tangent plane's third axis is flipped in the direction of z.

Projects angles from the x-axis around the z-axis into the tangent
plane, along the z-axis. The projected angle similarly refers to the 
tangent plane's primary axis.

Reversible transformation, see `angle_tangent_to_xy!`
"""
function angle_xy_to_tangent!(vϕ, vα, P)
    a, _, c, d, e, f, g, h, i = P
    @assert i !== 0
    map!(vϕ, vα) do α
        co = cos(α)
        si = sin(α)
        z = -(g·co + h·si) / i
        u = a·co +        c⋅z   # component along 𝐞₁
        v = d·co + e·si + f·z  # component along 𝐞₂
        ϕ = atan(v, u)
    end
end

"""
    angle_tangent_to_xy!(vα, vϕ, P)

This works with tangent bases with the condition that
the tangent plane's primary axis lies in the yz plane, and 
the tangent plane's third axis is flipped in the direction of z.

Projects angles from the tangent plane's primary axis around its third 
axis into the xy plane, along the z-axis. The projected angle similarly
refers to the x axis.

Reversible transformation, see `angle_xy_to_tangent!`
"""
function angle_tangent_to_xy!(vα, vϕ, P)
    # This works with tangent bases with 
    # the conditions we set on orientation.
    map!(vα, vϕ) do ϕ
        # Tangent-space vector [u, v, w]
        u = cos(ϕ)
        v = sin(ϕ)
        # Standard-basis vector 𝐩 = P ⋅ 𝐮
        x = P[1]·u + P[4]·v
        y = P[5]·v
        α = atan(y, x)
    end
end



"""
    principal_curvature_and_direction(vκ::T, vϕ::T) where T <: SVector{4, Float64}
    principal_curvature_and_direction(vκ, vφ) 
    -> (Float64, Float64, Float64)

Compute the principal curvatures κ₁ ≥ κ₂ and the principal direction angle φₚ from three curvature samples.

Given measurements κ(φᵢ) at angles φᵢ (in a common tangent-plane frame), we model

    κ(φ) = κ₁ cos(φ - φₚ)^2 + κ₂ sin(φ - φₚ)^2

where φₚ is the angle (in the same frame) at which the maximum curvature κ₁ occurs.

Expanding via trigonometric identities:
    cos^2(φ - φₚ) = [1 + cos(2φ - 2φₚ)]/2
    sin^2(φ - φₚ) = [1 - cos(2φ - 2φₚ)]/2

leads to the equivalent form

    κ(φ) = a + b cos(2φ) + c sin(2φ)
    where
      a = (κ₁ + κ₂) / 2
      r = (κ₁ - κ₂) / 2  (so b = r cos(2φₚ), c = r sin(2φₚ))

We solve the linear system

    A * [a b c]ᵀ = vκ,

where

    A = [1  cos(2φ₁)  sin(2φ₁);
         1  cos(2φ₂)  sin(2φ₂);
         1  cos(2φ₃)  sin(2φ₃)
         1  cos(2φ₄)  sin(2φ₄)]
      = hcat(ones(Float64, 4), cos.(2vφ), sin.(2vφ))

and solve with an optimized equivalent to:

    a, b, c = A \\ vκ

Then:
    r = √(b² + c²)
    κ₁ = a + r
    κ₂ = a - r
    φₚ = ½ ⋅ atan2(c, b)

Ensuring φₚ is in the correct branch of [0, 2π).

# Parameters
- vκ    Three curvature samples κ(φᵢ)
- vφ    Three angles φᵢ in the same tangent-plane frame

# Returns
- κ₁    Maximum principal curvature
- κ₂    Minimum principal curvature
- φₚ    Principal-direction angle (same units/frame as vφ)


# References
1. https://en.wikipedia.org/wiki/Euler%27s_theorem_(differential_geometry)
"""
function principal_curvature_and_direction(vκ::T, vϕ::T) where T <: SVector{4, Float64}
    # TODO: Non-allocating, mutating version.
    @assert length(vκ) == length(vϕ) == 4
    # Construct design matrix M and solve M * [a; b; c] = vκ
    A = SMatrix{4, 3, Float64, 12}(hcat(ones(Float64, size(vϕ, 1)), cos.(2vϕ), sin.(2vϕ)))
    # Modify vκ to fit the equations
    cache = init(LinearProblem{true}(A, vκ), NormalCholeskyFactorization(), OperatorAssumptions(false, condition=OperatorCondition.WellConditioned))
    a, b, c = solve!(cache)
    # Principal‐curvature values from a, b, c
    r = hypot(b, c)       # faster than √(b^2 + c^2)
    κ1 = a + r            # maximum principal curvature
    κ2 = a - r            # minimum principal curvature
    # First maximum of I) equals κ1 and occurs at (from differentiation)
    ϕp = atan(c, b) / 2
    if isnan(ϕp)
        throw(ErrorException("Unexpected NaN: vκ, vϕ, a, b, c = $vκ $vϕ $a $b $c"))
    end
    κ1, κ2, ϕp
end
principal_curvature_and_direction(vκ, vϕ) = principal_curvature_and_direction(SVector{4, Float64}(vκ), SVector{4, Float64}(vϕ))

"""
    allocations_curvature(R::CartesianIndices, directions; maxg = 50, ming = -50)

Allocate once, re-use at every point! This is re-used in calculations which require fewer 
such containers than curvature.

# Arguments

R             e.g. CartesianIndices(size(z))
directions    e.g. [] or [2] or [1, 2]

# Output

Ri, Ω, v, P, K, vα, vκ, vβ, f_is_within_limits

where

typeof(Ri) = CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}
typeof(Ω) = CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}
typeof(v) = StaticArraysCore.MVector{2, Float64}
typeof(P) = StaticArraysCore.MMatrix{3, 3, Float64, 9}
typeof(K) = StaticArraysCore.MMatrix{2, 2, Float64, 4}
typeof(vα) = StaticArraysCore.MVector{4, Float64}
typeof(vκ) = StaticArraysCore.MVector{4, Float64}
typeof(vβ) = StaticArraysCore.MVector{2, Float64}
f_is_within_limits isa Function, depending on argument directions.
"""
function allocations_curvature(R::CartesianIndices, directions; maxg = 50, ming = -50)
    @assert length(directions) <= 2
    # Define an internal domain Ri, since simple padding options
    # would not yield interesting curvature anyway.
    n = 5
    side = n ÷ 2
    is, js = side + 1, side + 1
    ie, je = size(R,1) - side - 1, size(R,2) - side - 1
    Ri = CartesianIndices((is:ie, js:je))
    # Relative indices, for defining a sliding window around each point
    Ω = CartesianIndices((-side:side, -side:side))
    # Pre-allocate mutable containers for re-using memory. 
    #   fixed-sized 2d vector, multiple purpose
    v = MVector{2, Float64}(Array{Float64, 1}(undef, 2))
    #   statically-sized, tangent basis
    P = MMatrix{3, 3, Float64}(Array{Float64, 2}(undef, (3, 3)))
    #  fixed-sized, curvature components 
    K = MMatrix{2, 2, Float64}(Array{Float64, 2}(undef, (2, 2)))
    #  fixed-sized, curvature sample angles in screen plane 
    vα = MVector{4, Float64}(Array{Float64, 1}(undef, 4))
    #  fixed-sized, sample curvature values  
    vκ = MVector{4, Float64}(Array{Float64, 1}(undef, 4))
    # Principal and secondary principal angles. 
    vβ = MVector{2, Float64}(Array{Float64, 1}(undef, 2))
    #
    f_is_within_limits = func_is_glyph_within_limits(directions, maxg, ming)
    Ri, Ω, v, P, K, vα, vκ, vβ, f_is_within_limits
end

