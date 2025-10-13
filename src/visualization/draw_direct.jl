# Drawing functionality for partial coverage
# Additional to `BitmapMaps.jl/mark_utils.jl`.
#
# Contains 
# `draw_vector!`
# `draw_bidirectional_vector!`
# `draw_points!` 
# `spray!`
# `LogMapper`
# `apply_color_by_coverage!`


"""
    draw_vector!(cov::Matrix{Float32}, A::CartesianIndex{2}, Δi::Int, Δj::Int, strength::Float32)

Draws a line from A to B in the given coverage matrix. Edges are blurry. At A, maximum radius is 

    `VECTOR_REL_HALFWIDTH` * `l`

where `l` = √(Δi²+Δj²). 

See LogMapper for conversion of `cov`` to an image.
"""
function draw_vector!(cov::Matrix{Float32}, A::CartesianIndex{2}, Δi::Int, Δj::Int, strength::Float32)
    #
    i1, j1 = Tuple(A)
    li = abs(Δi)
    lj = abs(Δj)
    si = sign(Δi)
    sj = sign(Δj)
    # Get all valid indices for the image
    R = CartesianIndices(cov)
    # Length
    l = hypot(Δi, Δj)
    # Max radius
    rA = Float32(max(0.5, VECTOR_REL_HALFWIDTH * l))
    # Min radius
    rB = min(0.5f0, rA)
    # Zero vector
    if Δi == 0 && Δj == 0
        spray!(cov,
                CartesianIndex(i1, j1),
                rB,
                strength)
        return cov
    end
    # Non-zero vector
    if li > lj
        # Radius as a function of vertical position
        ρ = let rA = rA, rB = rB, nl = li
            i -> rA + (i - 1) * (rB - rA) / nl
        end
        err = li ÷ 2
        for i in 1:(li + 1)
            r = ρ(i)
            if isnan(r)
                @show  Δi  Δj rA rB li lj r
                throw("ah NaN")
            end 
            spray!(cov,
                CartesianIndex(i1, j1),
                r,
                strength)
            err -= lj
            if err < 0
                j1 += sj
                err += li
            end
            i1 += si
        end
    else
        # Radius as a function of horizontal position
        ρ = let rA = rA, rB = rB, nl = lj
            i -> rA + (i - 1) * (rB - rA) / nl
        end
        err = lj ÷ 2
        for i in 1:(lj + 1)
            r = ρ(i)
            if isnan(r)
                @show  Δi  Δj rA rB li lj r
                throw("ah, what happened?")
            end 
            spray!(cov,
                CartesianIndex(i1, j1),
                r,
                strength)
            err -= li
            if err < 0
                i1 += si
                err += lj
            end
            j1 += sj
        end
    end
    cov
end

"""
    draw_bidirectional_vector!(cov, p, v::AbstractVector, strength::Float32)

Opposing vectors. Outwards pointing means positive, inwards pointing means negative,
as is customary for e.g. tensile / compressive force.

- `v`[1:2] specifies both direction and sign and is given in an x-y up coordinate system. 

`v` points in direction `θ` from x around z in an right-handed coordinate system. Then:

    θ ∈ [0, π>:  Positive sign.   
    θ ∈ [π, 2π>: Negative sign.   
"""
function draw_bidirectional_vector!(cov, p, v::AbstractVector, strength::Float32)
    Δj = Int(round(v[1]))
    Δi = Int(round(v[2]))
    if is_bidirec_vect_positive(v)
        # First or second quadrant. Positive.
        draw_vector!(cov, p, -Δi, Δj, strength)
        draw_vector!(cov, p, Δi, -Δj, strength)
    else
        # Third or fourth quadrant. Negative.
        draw_vector!(cov, p + CartesianIndex(-Δi, Δj), Δi, -Δj, strength)
        draw_vector!(cov, p + CartesianIndex(Δi, -Δj), -Δi, Δj, strength)
    end
end


















"""
    draw_streamlines_points!(cov, spts, r::Float32, strength::Float32) 

`spts` is a nested array of CartesianIndex{2}. The outer is a collection of points in each streamline.    

The shape and ordering of streamlines do not currently matter,  but an extension might be
used to convey some property.
"""
function draw_streamlines_points!(cov, spts, r::Float32, strength::Float32) 
    @assert eltype(spts) <: Vector
    for streamline in spts
        draw_streamline_points!(cov, streamline, r, strength)
    end
    cov
end

"""
    draw_streamline_points!(cov, spts::Vector{CartesianIndex{2}}, r::Float32, strength::Float32)
"""
function draw_streamline_points!(cov, spts::Vector{CartesianIndex{2}}, r::Float32, strength::Float32)
    for pt in spts
        spray!(cov, pt,  r, strength)
    end
    cov
end


"""
    spray!(cov::Matrix{Float32},
                centre::CartesianIndex{2},
                r::Float32,
                strength::Float32)

Accumulate one spray hit on the coverage buffer `cov`
(a `Matrix{Float32}`).

* `centre`   – CartesianIndex of the hit  
* `r`        – radius ∈ (0, 7] (Float32)  
* `strength` – ink load ≥ 0 

Coverage added to each in‑circle pixel at offset (dx, dy):

    w = strength * max(0, 1 - 0.8*(dx^2+dy^2)/r^2)

This is modified for 0 < r < 1 to keep additional coverage roughly proportional to r^2.

Returns the mutated `cov` (for chaining).
"""
function spray!(cov::Matrix{Float32},
                centre::CartesianIndex{2},
                r::Float32,
                strength::Float32)
    @assert 0f0 < r ≤ 100f0        "radius r must be in (0,100], is $r"
    @assert strength ≥ 0f0         "strength must be non‑negative"
    m   = Int(floor(r))
    r2  = r^2
    # Note that the number of pixels affected
    # increases and decreases in rough steps with r.
    # Compensate strength to make 
    # applied coverage roughly linear with r.
    strength_eff = strength * if r2 < 1f0
        8 * r^3f0 
    elseif r2 < 2
        # at r = 1 , target  coverage 1.0 
        # at r = √2, target coverage 2.66
        0.645f0 * r^3f0 + 2.1f0 #1.0 * r^3f0 + 1.1f0
    else
        1f0 + r
    end
    r2min = (max(0f0, r - 2f0))
    i₀, j₀ = Tuple(centre)
    H,  W  = size(cov)
    @inbounds for dy in -m:m, dx in -m:m
        d2 = dx*dx + dy*dy
        d2 > r2 && continue     # outside
        d2 < r2min && continue  # skip centre
        w = strength_eff * max(0, 1 - 0.8 * d2 / r2)
        w ≤ 0f0 && continue
        i  = clamp(i₀ + dy, 1, H)
        j  = clamp(j₀ + dx, 1, W)
        cov[i, j] += w
    end
    cov
end

"""
struct LogMapper{T}
    scale::T
    offset::T
end

When spraying we add hits linearly, which soon accumulates to more than one.
Tweaking 'scale' and 'offset' affects if one will be able to distinguish between
one, two or more hits.

"""
struct LogMapper{T}
    scale::T
    offset::T
end
# Default constructor
LogMapper() = LogMapper(N0f8(1/log1p(10)), 0N0f8)
# Callable
@inline (lm::LogMapper{T})(c) where T = T(clamp((log1p(c) * lm.scale) + lm.offset, 0, 1))


"""
    apply_color_by_coverage!(img, cov::Matrix{Float32}, color::RGB{N0f8})

Apply 'cov' on top of 'img' in `color`.

`img` is an RGB or RGBA image on which we apply a one-color overlay.
`cov` represents the coverage at each pixel, think of each element value as the 
     number of times a spray brush passes over this pixel. Also see `spray!`

The overlay is a temporary RGBA image where RGB is set by `color`. 
The overlay's opaqueness (the 'A' channel) is 0 where coverage is 0, 
meaning that 'img' is unchanged by the overlay. 

Where the coverage value is 10 or above, `img` is completely covered by the overlay:

| Coverage | Overlay opaqueness |
|----------|---------------|
| 0.0      | 0.0           |
| 0.5      | 0.168627      |
| 1.0      | 0.286275      |
| 1.5      | 0.380392      |
| 2.0      | 0.454902      |
| 2.5      | 0.521569      |
| 3.0      | 0.576471      |
| 3.5      | 0.623529      |
| 4.0      | 0.670588      |
| 4.5      | 0.709804      |
| 5.0      | 0.745098      |
| 5.5      | 0.776471      |
| 6.0      | 0.807843      |
| 6.5      | 0.839216      |
| 7.0      | 0.862745      |
| 7.5      | 0.890196      |
| 8.0      | 0.913725      |
| 8.5      | 0.937255      |
| 9.0      | 0.956863      |
| 9.5      | 0.976471      |
| 10.0     | 0.996078      |
 
"""
function apply_color_by_coverage!(img, cov::Matrix{Float32}, color::RGB{N0f8})
    mapper = LogMapper()
    f = x -> RGBA{N0f8}(color.r, color.g, color.b, mapper(x))
    # Composite over img
    img .= blend.(img, f.(cov))
end

"""
    apply_color_by_any_coverage!(img, cov::Matrix{Float32}, color::RGB{N0f8})

See `apply_color_by_coverage!`, but here, coverage is on-off.


| Coverage | Overlay opaqueness |
|----------|---------------|
| 0.0      | 0.0           |
| 0.01     | 1.0           |
| 123.45   | 1.0           |

"""
function apply_color_by_any_coverage!(img, cov::Matrix{Float32}, color::RGB{N0f8})
    mapper = x -> x > 0 ? 1 : 0
    f = x -> RGBA{N0f8}(color.r, color.g, color.b, mapper(x))
    # Composite over img
    img .= blend.(img, f.(cov))
end