# We can simplify the quantified curvature K a bit by
# categorizing each of the principal directions to +, ~0 and -.
# See `convexity_rank`. 
#
# Relies on constant PALETTE_CONVEXITY
# TODO: Move some of these functions to a better place.



_order_rank_123(c1, c2) = (c1 - 1) * (8 - c1) ÷ 2 + (c2 - c1 +1)

"""
    convexity_rank(s::Float64, flatval)
    --> Int 1-3
    convexity_rank(K, flatval)
    --> Int 1-6

Lexicographic rank, not "matrix rank".

# Arguments

- `s`       : signed curvature component along principal direction
- `flatval` : Flat is when  `abs(curvature) < flatval`
- `K` is a 2 x 2 matrix, format defined by `principal_curvature_components!`


| Signed curvature component | Rank   | Description |
| -------------------------- | -------| ----------- |
|            s > flatval     | 1      |  Convex     |
| -flatval < s < flatval     | 2      |  Flat       |
|            s < - flatval   | 3      |  Concave    |


Rank of principal curvature components (`K`):


| Major rank | Minor rank | Rank   | Description       |
| ---------- | ---------- | ------ | ----------------  |
|  1         | 1          |  1     | Convex  - convex  |
|  1         | 2          |  2     | Convex  - flat    |
|  1         | 3          |  3     | Convex  - concave |
|  2         | 2          |  4     | Flat    - flat    |
|  2         | 3          |  5     | Flat    - Concave |
|  3         | 3          |  6     | Concave - concave |

"""
convexity_rank(s::Float64, flatval) = s > flatval ? 1 : s < -flatval ? 3 : 2
function convexity_rank(K, flatval)
    # K is ordered: Coloumn 1 is major, column 2 is minor.
    s1, s2 = signed_curvature_values(K)
    @assert s1 >= s2 "s1 $s1 >= s2 $s2"
    c1 = convexity_rank(s1, flatval)
    c2 = convexity_rank(s2, flatval)
    # Ordered combination no., 1-6
    _order_rank_123(c1, c2)
end


"""
    signed_curvature_values(K)
    --> {Float64, Float64}
"""
function signed_curvature_values(K)
    s1 = norm(K[:, 1]) * (is_bidirec_vect_positive(K[:, 1]) ? 1 : -1)
    s2 = norm(K[:, 2]) * (is_bidirec_vect_positive(K[:, 2]) ? 1 : -1)
    s1, s2
end
# Graphical output 

"""
paint_convexity_rank!(img::Matrix{<:RGB}, z, pts;  
    flatval = 0.0003, palette = PALETTE_CONVEXITY)

Overlays a color-map onto the image img.

See callees `color_point_by_convexity_rank!` and  `convexity_rank`!
"""
function paint_convexity_rank!(img::Matrix{<:RGB}, z, pts;  
    flatval = 0.0003, palette = PALETTE_CONVEXITY) 
    #
    buf = zeros(RGB{N0f8}, size(img)...)
    _paint_convexity_rank!(buf, z, pts, flatval, palette)
    # Take hue and chroma from buf
    chromaticity_over!(img, buf)
    img
end 

function _paint_convexity_rank!(buf, z, pts, flatval, palette)
    # Prepare
    bog = BidirectionOnGrid(𝐊!, z)
    R = CartesianIndices(z)
    Ω = CartesianIndices((-2:2, -2:2))
    minj = R[1][2] - Ω[1][2]
    maxi = R[end][1] - Ω[end][2]
    maxj = R[end][2] - Ω[end][2]
    mini = R[1][1] - Ω[1][2]
    Ri = CartesianIndices((mini:maxi, minj:maxj))
    # Color points one at a time
    for pt in filter(pt -> pt ∈ Ri, sort(vec(pts)))
        color_point_by_convexity_rank!(buf, pt, bog(pt.I...), flatval, palette)
    end
    buf
end


"""
    color_point_by_convexity_rank!(buf, pt,  K, flatval, palette)

K is a 2x2 matrix with two bidirectional vectors, returned from `components_matrix!`.
"""
function color_point_by_convexity_rank!(buf, pt,  K, flatval, palette)
    i = convexity_rank(K, flatval)
    buf[pt] = palette[i]
end


