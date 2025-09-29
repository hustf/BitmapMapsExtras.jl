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
     @assert gs.directions == 1:2 "One direction not implemented. "
    # Extract the half-length of primary and secondary principal direction
    # glyphs
    if is_in_limits(gs, K)
        K1 = @view K[:, 1]
        K2 = @view K[:, 2]
        l1 = norm(K1) * gs.multip
        l2 = norm(K2) * gs.multip
    else
        # We return 0, because we do not want
        # this potential glyph to be one of the selected ones.
        l1 = 0.0
        l2 = 0.0
    end
    max(l1, l2) 
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
    a = hypot(x1, y1)             # major radius
    b = hypot(x2, y2)             # minor radius (approx)
    @assert a >= b
    θ = atan(y1, x1)              # orientation
    ϕ = α - θ
    denom = (cos(ϕ)^2)/a^2 + (sin(ϕ)^2)/b^2
    return 1 / sqrt(denom)
end

"""
    intersection_distance_glyph(gs::GSTensor, K, α)
    intersection_distance_glyph(gs::GSVector, v, α)

Approximate distance from origin along ray at angle `α` to eclipsing elliptical boundary .

Here, we assume that glyph values are within limits.
"""
function intersection_distance_glyph(gs::GSTensor, K, α)
     @assert gs.directions == 1:2 "One direction not implemented. "
    # Extract the half-length of primary and secondary principal direction
    # glyphs.
    # 
    # (x1, y1): Longest vector
    # (x2, y2): Shortest vector
    # The local 'y' is in a right handed, y up, x right coordinate system
    if norm(K[:, 1]) >= norm(K[:, 2])
        x1, x2 = K[1, :]
        y1, y2 = -K[2, :]
    else
        x2, x1 = K[1, :]
        y2, y1 = -K[2, :]
    end
    l1 = hypot(x1, y1)
    l2 = hypot(x2, y2)
    if l2 / l1 < 5 * VECTOR_REL_HALFWIDTH
        # The minor axis is so short compared to the other axis
        # that the determining 'width' is the half-cone axis of the vector glyph.
        # Since we don't want completely neighbouring parallel vectors,
        # we also include a factor of 3.
        # The direction is set to 90° off (x1, y1)
        #x2 = -y1 * VECTOR_REL_HALFWIDTH
        #y2 = x1 * VECTOR_REL_HALFWIDTH
        x2norm = x2 / l2
        y2norm = y2 / l2
        minlen = 5 * l1 * VECTOR_REL_HALFWIDTH
        x2 = x2norm * minlen
        y2 = y2norm * minlen
    end
    idist =  gs.multip * intersection_distance_ellipse(x1, y1, x2, y2, α)
    # 
end

"""
    intersection_distance_glyph(gs::GSVector, v, α)
    intersection_distance_glyph(::Type{GSVector}, multip, v, α)

v is the 2d vector.
α is the angle from screen horizontal (j or x-axis) and a ray to an external point.
"""
intersection_distance_glyph(gs::GSVector, v, α) = intersection_distance_glyph(GSVector, gs.multip, v, α)
function intersection_distance_glyph(::Type{GSVector}, multip, v, α)
    Δj = Int(round(multip * v[1]))
    Δi = -Int(round(multip * v[2]))
    l = multip * norm(v)
    # Max, min radius
    rA = Float32(max(0.5, VECTOR_REL_HALFWIDTH * l))
    rB = min(0.5f0, rA)
    # The glyph points in direction v, (x, y) coordinates
    # Glyph points in angle 
    β = atan(v[2], v[1]) 
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
    too_close(gs::U, pt, coarse_r, K, ptaccepted, coarse_raccepted, Kaccepted) where U<:GlyphSpec

..If the potential glyph at pt would intersect with the
already selected glyph at ptaccepted.

A coarse filtering should be done first. Here, we assume
that glyph values are within limits.
"""
function too_close(gs::U, pt, coarse_r, K, ptaccepted, coarse_raccepted, Kaccepted) where U<:GlyphSpec
    # Vector from pt to ptaccepted
    u = ptaccepted.I .- pt.I
    # Length of u
    lu = norm(u)
    # Exclusive criterion - fast return
    nooverlap = coarse_r + coarse_raccepted <= lu  
    nooverlap && return false
    # Vector u's angle to horizontal 
    αu = atan(-u[1], u[2]) 
    # Intersection dist with pt's boundary along u
    r = intersection_distance_glyph(gs, K, αu)
    # Intersection with already accepted glyph's boundary
    # A ray from ptaccepted to pt intesects the accepted glyph at distance
    raccepted = intersection_distance_glyph(gs, Kaccepted, αu + π)
    # The glyphs overlap when
    overlap1 = r + raccepted >  lu
    # (This criterion would be slightly different, since it does not compare with the
    #    shortest distance between anchor points, but between 'glancing points': 
    #    too_close_glancing(gs, pt, K, ptaccepted, Kaccepted, αu, 0)
    # )
    overlap1 && return true
    #
    # The order of checks 2 to 7 is optimized for a certain collection of vector glyphs.
    overlap3 = too_close_glancing(gs, pt, K, ptaccepted, Kaccepted, αu, π / 2)
    overlap3 && return true
    #
    # Main axis of accepted glyph. K is given in (x,y) frame
    αmaccepted = atan(Kaccepted[2], Kaccepted[1])
    overlap7 = too_close_glancing(gs, ptaccepted, Kaccepted, pt, K,  αmaccepted, 0)
    overlap7 && return true
    #
    overlap2 = too_close_glancing(gs, pt, K, ptaccepted, Kaccepted, αu, -π / 2)
    overlap2 && return true
    #
    overlap5 = too_close_glancing(gs, ptaccepted, Kaccepted, pt, K,  αu, π / 2)
    overlap5 && return true
    #
    # Main axis of glyph. K is given in (x,y) frame
    αm = atan(K[2], K[1])
    overlap6 = too_close_glancing(gs, pt, K, ptaccepted, Kaccepted, αm, 0)
    overlap6 && return true
    #
    overlap4 = too_close_glancing(gs, ptaccepted, Kaccepted, pt, K, αu, -π / 2)
    overlap4 && return true
    #
    false
end
function too_close_glancing(gs::U, pt, K, ptaccepted, Kaccepted, αu, Δαu) where U<:GlyphSpec
    # Less simple overlap calculation, in addition to checking 'straight line overlap'
    # This checks for 'glancing' overlaps:
    α = αu + Δαu
    # Distance from pt's anchor point to the intersection with pt's glyph boundary
    r = intersection_distance_glyph(gs, K, α)
    # The intersection with pt's glyph boundary
    # Floating point, not CartesianIndex
    y = pt.I[1] - r * sin(α)
    x = pt.I[2] + r * cos(α)
    # q: From the intersection with pt's glyph boundary to ptaccepted's anchor point.
    qy = ptaccepted.I[1] - y
    qx = ptaccepted.I[2] - x
    # Angle of vector q
    β = atan(-qy, qx)
    # Length of vector q
    lq = hypot(qy, qx)
    # The length of a ray from the accepted glyph's anchor point along q 
    # which is inside of the accepted glyph.
    raccepted = intersection_distance_glyph(gs, Kaccepted, β + π)
    #printstyled("pt $pt K $K ptaccepted $ptaccepted Kaccepted $Kaccepted rad2deg(αu) $(rad2deg(αu)) x $x y $y rad2deg(β) $(rad2deg(β)) r $r raccepted $raccepted lq $lq qx $qx qy $qy" ,"\n", color = 176)
    # Criterion for overlap (r would appear on both sides and is dropped)
    raccepted >  lq
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
    #@show length(passed_placements)
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
