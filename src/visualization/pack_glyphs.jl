# Smartly spaced glyphs
"""
    potential_scattered_placements(R::CartesianIndices; avgdist = 10, seed = MersenneTwister(123))
    --> Vector{CartesianIndex}

Return a vector of random CartesianIndex scattered within `R` such that the
mean nearest-neighbour distance is approximately `avgdist`. Refer Poisson disk sampling.

# Arguments

- `R` defines the possible values of each returned index. 
- `avgdist` (default keyword value 10.0) is the mean nearest-neighbour distance between points. The 
    resulting set of points will vary slightly from this target.
"""
function potential_scattered_placements(R::CartesianIndices; avgdist = 10.0, seed = MersenneTwister(123))
    area = length(R)
    # Estimate number of points, within bounds.
    n = max(1, min(round(Int, area / (4 * avgdist^2)), area))
    # Draw random unique indices
    chosen = randperm(seed, area)[1:n]
    [R[i] for i in chosen]
end

function potential_scattered_placements(fz::T, gs::U; scatterdist = gs.dashsize, seed = MersenneTwister(123)) where
    {T<: Union{BidirectionOnGrid, DirectionOnGrid}, U<:GlyphSpec}
    #
    R = CartesianIndices(fz.z)
    Ω = fz.Ω
    minj = R[1][2] - Ω[1][2]
    maxi = R[end][1] - Ω[end][2]
    maxj = R[end][2] - Ω[end][2]
    mini = R[1][1] - Ω[1][2]
    Ri = CartesianIndices((mini:maxi, minj:maxj))
    potential_scattered_placements(Ri; avgdist = scatterdist, seed)
end

function coarse_radius_for_plotting(gs::GSTensor, K)
    # Extract the half-length of primary and secondary principal direction
    # glyphs
    if is_in_limits(gs, K)
        if length(gs.directions) == 1
            v = @view K[:, first(gs.directions)]
            l1 = norm(v) * gs.multip
            if is_bidirec_vect_positive(v)
                # Add min radius
                Δl = min(0.5, VECTOR_REL_HALFWIDTH * l1)
            else
                # Add max radius
                Δl = max(0.5, VECTOR_REL_HALFWIDTH * l1)
            end
            dist = l1 + Δl
        else
            v1 = @view K[:, 1]
            v2 = @view K[:, 2]
            l1 = norm(v1) * gs.multip
            l2 = norm(v2) * gs.multip
            if is_bidirec_vect_positive(v1)
                # Add min radius
                Δl1 = min(0.5, VECTOR_REL_HALFWIDTH * l1)
            else
                # Add max radius
                Δl1 = max(0.5, VECTOR_REL_HALFWIDTH * l1)
            end
            if is_bidirec_vect_positive(v2)
                # Add min radius
                Δl2 = min(0.5, VECTOR_REL_HALFWIDTH * l2)
            else
                # Add max radius
                Δl2 = max(0.5, VECTOR_REL_HALFWIDTH * l2)
            end
            dist = max(l1 + Δl1, l2 + Δl2)
        end
    else
        # We return 0, because we do not want
        # this potential glyph to be one of the selected ones.
        dist = 0.0
    end
    # Anchor point to tip
    return dist
end
function coarse_radius_for_plotting(gs::GSVector, v)
    if is_in_limits(gs, v)
        norm(v) * gs.multip
    else
        Float64(gs.dashsize)
    end
end

"""
    intersection_distance_ellipse(x1, y1, x2, y2, α)

Approximate distance from origin to ellipse boundary along ray at angle `α`.

- `(x1, y1)` is assumed to lie on the major axis (used for orientation + length).
- `(x2, y2)` only contributes its distance as semi-minor axis length.
- `α` is ray angle in radians (CCW from x-axis).
"""
function intersection_distance_ellipse(x1, y1, x2, y2, α)
    @warn "Drop this dead" maxlog = 100
    throw("dead")
    a = hypot(x1, y1)             # major radius
    b = hypot(x2, y2)             # minor radius (approx)
    @assert a >= b
    θ = atan(y1, x1)              # orientation
    ϕ = α - θ
    denom = (cos(ϕ)^2)/a^2 + (sin(ϕ)^2)/b^2
    return 1 / sqrt(denom)
end

"""
    mirror_to_right_half(θ)
"""
function mirror_to_right_half(θ)
    #  into [-π, π]
    θ1 = θ > π ?  θ - 2π : θ   
    sign(θ1) * min(abs(θ1), π - abs(θ1))
end

"""
    radial_distance_bidirectional_glyph(::Type{GSVector}, vscaled, α)

See `radial_distance_glyph`.
"""
function radial_distance_bidirectional_glyph(::Type{GSVector}, vscaled, α)
    l = norm(vscaled)
    # Max, min radius
    if is_bidirec_vect_positive(vscaled)
        rA = Float32(max(0.5, VECTOR_REL_HALFWIDTH * l))
        rB = min(0.5f0, rA)
    else
        rB = Float32(max(0.5, VECTOR_REL_HALFWIDTH * l))
        rA = min(0.5f0, rB)
    end
    # The glyph points along vscaled (given in x-y-up directions)
    # Angle of 'pointing axis' rel to x (in x-y-up c.s) 
    β = atan(vscaled[2], vscaled[1])
    # Angle center to tip to side 
    δ = atan(rA - rB, l) 
    # Angle from the glyph direction to α, the line to the other placement.
    # We add 2π where necessary to get positive angles.
    # Since this is a bidirectional, symmetric around γ = π / 2,
    # we're mirroring angles in the 2nd and 3d quadrant
    γ = mod2pi(mirror_to_right_half(α - β))
    # Angle from tip center to tip tangent point
    η = atan(rB, l) 
    @assert rA > 0
    @assert rB > 0
    if 0 < γ <= η
        dist = l + rB
    elseif η < γ <= π / 2 - δ
        # Angle between ray and point on circle rA where the side of the arrow is tangent
        θ = π / 2 - δ - γ
        dist = rA / cos(θ)
    elseif π / 2 - δ < γ <= 3π / 2 + δ
        dist = rA
    elseif 3π / 2 + δ < γ <= 2π - η
        # Angle between ray and point on circle rA where the side of the arrow is tangent
        θ = 3π / 2 + δ - γ
        dist = rA / cos(θ)
    else 
        dist = l + rB
    end
    @assert min(rA, rB) <=  dist
    dist
end
"""
    radial_distance_glyph(gs::GSTensor, K, α)
    radial_distance_glyph(gs::GSVector, v, α)
    radial_distance_glyph(::Type{GSVector}, vscaled, α)

Distance from a glyph's anchor point along a ray at angle `α` to eclipsing glyph boundary. 
The shape of the glyph is strictly convex.

`α` is the counter-clockwise angle from x in 'x-y up' coordinates between screen horizontal (j or x-axis) and a ray to a point on the circumference.
`v` is the 2-d vector.
`vscaled` is `v * gs.multip`
`K` is a TENSORMAP, i.e. two bidirectional and signed vectors.

Here, we assume that glyph values are within limits.
"""
function radial_distance_glyph(gs::GSTensor, K, α)
    if length(gs.directions) == 1
        dirno = first(gs.directions)
        # Scale the single bidirectional vector's value
        vscaled = gs.multip .* view(K, :, dirno)
        return radial_distance_bidirectional_glyph(GSVector, vscaled, α)
    else
        @assert gs.directions == 1:2
        # Scale the two glyph values
        vs1 = gs.multip .* view(K, :, 1)
        vs2 = gs.multip .* view(K, :, 2)
        lv1 = radial_distance_bidirectional_glyph(GSVector, vs1, α)
        lv2 = radial_distance_bidirectional_glyph(GSVector, vs2, α)
        return max(lv1, lv2)
    end
end
function radial_distance_glyph(gs::GSVector, v, α)
    radial_distance_glyph(GSVector, gs.multip .* v, α)
end
function radial_distance_glyph(::Type{GSVector}, vscaled, α)
    l = norm(vscaled)
    # Max, min radius
    rA = Float32(max(0.5, VECTOR_REL_HALFWIDTH * l))
    rB = min(0.5f0, rA)
    # The glyph points along vscaled (given in x-y-up directions)
    # Angle of 'pointing axis' rel to x (in x-y-up c.s) 
    β = atan(vscaled[2], vscaled[1]) 
    # Angle center to tip to side 
    δ = atan(rA - rB, l) 
    # Angle from the glyph direction to α, the line to the other placement.
    # We add 2π where necessary to get positive angles.
    γ = mod(α - β, 2π)
    # Angle from tip center to tip tangent point
    η = atan(rB, l) 
    @assert rA > 0
    @assert rB > 0
    if 0 < γ <= η
        dist = l + rB
    elseif η < γ <= π / 2 - δ
        # Angle between ray and point on circle rA where the side of the arrow is tangent
        θ = π / 2 - δ - γ
        dist = rA / cos(θ)
    elseif π / 2 - δ < γ <= 3π / 2 + δ
        dist = rA
    elseif 3π / 2 + δ < γ <= 2π - η
        # Angle between ray and point on circle rA where the side of the arrow is tangent
        θ = 3π / 2 + δ - γ
        dist = rA / cos(θ)
    else 
        dist = l + rB
    end
    @assert min(rA, rB) <=  dist
    dist
end

"""
    too_close(gs::U, pt1, coarse_r, v1, pt2, coarse_r2, v2) where U<:GlyphSpec

..If the potential glyph at pt1 would intersect with the
already selected glyph at pt2.

A coarse filtering should be done first. Here, we assume
that glyph values are within limits. Input should be pre-filtered.
"""
function too_close(gs::U, pt1, coarse_r, v1, pt2, coarse_r2, v2) where U<:GlyphSpec
    # Vector from pt1 to pt2
    u = pt2.I .- pt1.I
    # Length of u
    lu = norm(u)
    # Exclusive criterion - fast return
    nooverlap = coarse_r + coarse_r2 <= lu  
    nooverlap && return false
    # Vector u's angle to horizontal in x-y-up directions
    αu = atan(-u[1], u[2]) 
    # Intersection dist with pt1's boundary along u
    r1 = radial_distance_glyph(gs, v1, αu)
    # Intersection with glyph 2's boundary
    # A ray from pt2 to pt1 intersects glyph 2's extents at distance
    r2 = radial_distance_glyph(gs, v2, αu + π)
    # The glyphs overlap when
    overlap_direct = r1 + r2 >  lu
    overlap_direct && return true
    # All other checks works the same way: Find an angled path 
    # between anchor points which is entirely within the two glyphs.
    # The interesting paths depend on the GlyphSpec type, and this trait:
    
    overlap_indirect(gs, pt1, v1, pt2, v2, αu)
end
function overlap_indirect(gs::U, pt1, v1, pt2, v2, αu) where U<:GlyphSpec
    #
    # The order of checks 2 to 7 is optimized for a certain collection of vector glyphs.
    overlap3 = exit_of_first_is_in_second(gs, pt1, v1, pt2, v2, αu, π / 2)
    overlap3 && return true
    #
    # Main axis of glyph 2. v2 is given in (x,y) frame
    αm2 = atan(v2[2], v2[1])
    overlap7 = exit_of_first_is_in_second(gs, pt2, v2, pt1, v1,  αm2, 0)
    overlap7 && return true
    #
    overlap2 = exit_of_first_is_in_second(gs, pt1, v1, pt2, v2, αu, -π / 2)
    overlap2 && return true
    #
    overlap5 = exit_of_first_is_in_second(gs, pt2, v2, pt1, v1,  αu, π / 2)
    overlap5 && return true
    #
    # Main axis of glyph 1. v1 is given in (x,y) frame
    αm1 = atan(v1[2], v1[1])
    overlap6 = exit_of_first_is_in_second(gs, pt1, v1, pt2, v2, αm1, 0)
    overlap6 && return true
    #
    overlap4 = exit_of_first_is_in_second(gs, pt2, v2, pt1, v1, αu, -π / 2)
    overlap4 && return true
    #
    if length(gs.directions) == 1
        overlap10 = axis_intersect_within(gs, pt1, v1, pt2, v2, αm1, αm2)
        overlap10 && return true
    else
        overlap11 = dual_axes_intersect_within(gs, pt1, v1, pt2, v2, αm1, αm2)
        overlap11 && return true
    end
    # Additional 'glancing' checks for bidirectional vectors
    if U <: GSTensor
        # Main axis + π
        overlap8 = exit_of_first_is_in_second(gs, pt2, v2, pt1, v1,  αm2, π)
        overlap8 && return true
        overlap9 = exit_of_first_is_in_second(gs, pt1, v1, pt2, v2, αm1, π)
        overlap9 && return true
    end
    false
end

function axis_intersect_within(gs, pt1, v1, pt2, v2, αm1, αm2)
    # Angle (in x-y-up directions) of ray from pt1
    α1 = αm1
    α2 = αm2
    # Image coordinates
    i1 = pt1.I[1]
    j1 = pt1.I[2]
    i2 = pt2.I[1]
    j2 = pt2.I[2]
    # We're not checking vertically aligned anchor points
    j2 == j1 && return false 
    # Parallell glyphs don't intersect
    α1 ≈ α2 && return false
    # The intersection (floating) point is at 3
    #=
    i3f = i1 - (j3f - j1) * tan(α1)
    i3f = i2 - (j3f - j2) * tan(α2)
    =>
    i1 - (j3f - j1) * tan(α1) = i2 - (j3f - j2) * tan(α2)
    =>
    (-j3f + j1) * tan(α1) = i2 - i1 - (j3f - j2) * tan(α2)
    =>
    -j3f * tan(α1) + j3f * tan(α2) = i2 - i1 + j2 * tan(α2) - j1 * tan(α1)
    =>
    j3f * (tan(α1) - tan(α2)) =  i1 - i2 + j1 * tan(α1) - j2 * tan(α2)
    =>
    j3f  =  (i1 - i2 + j1 * tan(α1) - j2 * tan(α2)) / (tan(α1) - tan(α2))
    =#
    j3f  = (i1 - i2 + j1 * tan(α1) - j2 * tan(α2)) / (tan(α1) - tan(α2))
    i3f = i1 - (j3f - j1) * tan(α1)
    # Distances from anchor points to intersection
    l1 = hypot(j3f - j1, i3f - i1)
    l2 = hypot(j3f - j2, i3f - i2)
    # The intersection angle 3 may well be 180° to α
    # Signed angles from anchor to intersection (x-y up) 
    # Note we could probably find the angle to intersection without trigonometry.
    α1s = atan(-(i3f - i1), j3f - j1)
    α2s = atan(-(i3f - i2), j3f - j2)
    # Distance from pt1's anchor point to exit (the intersection with pt1's glyph boundary)
    r1 = radial_distance_glyph(gs, v1, α1s)
    r2 = radial_distance_glyph(gs, v2, α2s)
    # @show l1 l2 rad2deg(α1) rad2deg(α1s) rad2deg(α2) rad2deg(α2s) r1 r2 
    # Criterion
    # @show  (l1 <= r1 && l2 <= r2)
    l1 <= r1 && l2 <= r2
end

function dual_axes_intersect_within(gs, pt1, K1::TENSORMAP, pt2, K2::TENSORMAP, αp1, αp2)
    @assert length(gs.directions) == 2
    @assert gs.directions[1] == 1
    # Primary - primary. Both are bidirectional 2d vectors.
    axis_intersect_within(gs, pt1, K1, pt2, K2, αp1, αp2) && return true
    # Prepare for the other checks. `v` is 'primary' and `w` is secondary direction
    v1 = K1[:, 1] # Direction αp1
    v2 = K2[:, 1] # Direction αp2
    w1 = K1[:, 2] # Direction αs1
    w2 = K2[:, 2] # Direction αs2
    # Secondary axes of glyphs. w is given in (x,y) frame
    αs1 = atan(w1[2], w1[1])
    αs2 = atan(w1[2], w1[1])
    # 
    # Secondary - secondary.
    axis_intersect_within(gs, pt1, K1, pt2, K2, αs1, αs2) && return true
    # Primary - secondary. 
    #@show K1 v1 w1
    #@show K1[2] v1[2] w1[2]
    axis_intersect_within(gs, pt1, K1, pt2, K2, αp1, αs2) && return true
    #
    false
end




"""
     exit_of_first_is_in_second(gs::U, pt1, v1, pt2, v2, αu, Δαu) where U<:GlyphSpec

Less simple overlap calculation detection.
Checks if a path from `pt1` to `pt2` is within glyphs 1 or 2.
The path starts a straight line from `pt1` to `exit 1` along direction `αu + Δαu`.
It ends with a straight section from exit to `pt2`.
"""
function exit_of_first_is_in_second(gs::U, pt1, v1, pt2, v2, αu, Δαu) where U<:GlyphSpec
    #
    # Angle (in x-y-up directions) of ray from pt1
    α = αu + Δαu
    # Distance from pt1's anchor point to exit (the intersection with pt1's glyph boundary)
    r1 = radial_distance_glyph(gs, v1, α)
    # The exit of the ray from pt1's glyph 
    # Floating point, not CartesianIndex
    i_f = pt1.I[1] - r1 * sin(α)
    j_f = pt1.I[2] + r1 * cos(α)
    # q: Vector from the ray's exit of glyph 1,
    # pointing to pt2's anchor point.
    qi = pt2.I[1] - i_f
    qj = pt2.I[2] - j_f
    # Angle of vector q (in x-y up) c.s.
    β = atan(-qi, qj)
    # Length of vector q
    lq = hypot(qi, qj)
    # The length of a ray from glyph 2's anchor point along q 
    # which is inside of glyph 2.
    r2 = radial_distance_glyph(gs, v2, β + π)
    # Criterion for overlap 
    #@show pt1 rad2deg(α) qj qi rad2deg(β) lq r2 r1  (r2 > lq)
    #println()
    r2 >  lq
end

"""
    placements_and_values(b, gs, ppts)

Given a functor b and a glyph spec gs, 
select a subset of ppts. These are 
placements for non-overlapping glyphs.

Also returns 

```
passed_values == b.(passed_placements)
```
"""
function placements_and_values(b::T, gs::U, ppts)  where {T<: Union{BidirectionOnGrid, DirectionOnGrid}, U<:GlyphSpec}
    passed_placements = CartesianIndex{2}[]
    valtyp = typeof(b(first(ppts).I...))
    # Could be a vector of vectors, or a vector of TENSORMAP
    passed_values = Vector{valtyp}()
    passed_radii = Float64[]
    se = Set(ppts)
    while !isempty(se)
        # Set current point, drop it from the set.
        pt = pop!(se)
        value = b(pt.I...)
        r = coarse_radius_for_plotting(gs, value)
        if r > MAG_EPS
            if isempty(passed_placements)
                # Only this point to consider: There must be room.
                push!(passed_placements, pt)
                push!(passed_values, copy(value))
                push!(passed_radii, r)
            else
                any_too_close = any(zip(passed_placements, passed_radii, passed_values)) do (ppt, pr, pv)
                    too_close(gs, pt, r, value, ppt, pr, pv)
                end
                # None too close => add this glyph as passed
                if ! any_too_close
                    push!(passed_placements, pt)
                    push!(passed_values, copy(value))
                    push!(passed_radii, r)
                end
            end
        end
    end
    # @show length(passed_placements)
    passed_placements, passed_values
end


"""
pack_glyphs!(img, z::Matrix{<:AbstractFloat}, gs::GSTensor; scatterdist = gs.dashsize, seed = MersenneTwister(123))
pack_glyphs!(img, z::Matrix{<:AbstractFloat}, gs::GSVector; scatterdist = gs.dashsize, seed = MersenneTwister(123))
pack_glyphs!(img, b::T, gs::U; scatterdist = gs.dashsize, seed = MersenneTwister(123))  where
    {T<: Union{BidirectionOnGrid, DirectionOnGrid}, U<:GlyphSpec}
    
"""
function pack_glyphs!(img, z::Matrix{<:AbstractFloat}, gs::GSTensor; scatterdist = gs.dashsize, seed = MersenneTwister(123))
    pack_glyphs!(img, BidirectionOnGrid(𝐊!, z), gs; scatterdist, seed)
end
function pack_glyphs!(img, z::Matrix{<:AbstractFloat}, gs::GSVector; scatterdist = gs.dashsize, seed = MersenneTwister(123))
    pack_glyphs!(img, DirectionOnGrid(𝐧ₚ!, z), gs; scatterdist, seed)
end

function pack_glyphs!(img, b::T, gs::U; scatterdist = gs.dashsize, seed = MersenneTwister(123))  where
    {T<: Union{BidirectionOnGrid, DirectionOnGrid}, U<:GlyphSpec}
    ppts = potential_scattered_placements(b, gs; scatterdist, seed)
    filtered_placements, filtered_values = placements_and_values(b, gs, ppts)
    plot_glyphs_given_values!(img, filtered_placements, filtered_values, gs)
    img
end

