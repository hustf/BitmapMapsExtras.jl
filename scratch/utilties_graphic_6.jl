# These, if proven useful, could
# possibly belong in BitmapMaps.jl or in a related package.

"""
    color_neighbors!(img::Matrix{GrayA{N0f8}}, R, C::CartesianIndex{2}, r)

Spray the area around a pixel, alpha channel variation is parabolic with distance.

All neighbors of index `C` within distance `r` and bounds `R` get 
value 1 and alpha channel varying from 1 at centre to 0 at r2. Parabolic variation.
"""
function color_neighbors!(img::Matrix{GrayA{N0f8}}, R, C::CartesianIndex{2}, r)
    # Pixel alpha function
    f = let r = r, α_cen = r < 1 ? r : 1
        (di, dj) -> α_cen - 0.8 * (di^2 + dj^2 ) / r^2
    end
    m = Int(floor(r))
    for di in -m:m
        for dj in -m:m
            neighbor = CartesianIndex(C[1] + di, C[2] + dj)
            if neighbor in R 
                # Not outside image.
                α1 = f(di, dj)
                if α1 >= 0.2
                    # We have a contribution for this pixel
                    # What's here already?
                    α0 = alpha(img[neighbor])
                    if α0 <= oneunit(N0f8)
                        img[neighbor] = GrayA{N0f8}(oneunit(N0f8), N0f8(min(1.0, α1 + float(α0))))
                    end
                end
            end
        end
    end
    img
end


"""
    vector!(img, A, B) = vector!(img, A, B[1] - A[1], B[2] - A[2])
    vector!(img, A::CartesianIndex{2}, Δi, Δj)

Draws a line from A to B in the given image, with a max thickness defined by `tol_dist`.
"""
function vector!(img, A::CartesianIndex{2}, Δi::Int, Δj::Int)
    #
    i1, j1 = Tuple(A)
    li = abs(Δi)
    lj = abs(Δj)
    si = sign(Δi)
    sj = sign(Δj)
    # Get all valid indices for the image
    R = CartesianIndices(img)
    # Length
    l = hypot(Δi, Δj)
    # Max radius
    rA = max(0.5, 0.075 * l)
    rB = min(0.5, rA)
    # Zero vector
    if Δi == 0 && Δj == 0
        color_neighbors!(img, R, CartesianIndex(i1, j1), rB)
        return img
    end
    # Non-zero vector
    if li > lj
        # Radius as a function of vertical position
        ρ = let rA = rA, rB = rB, nl = li
            i -> rA + (i - 1) * (rB - rA) / nl
        end
        err = li ÷ 2
        for i in 1:(li + 1)
            if isnan(ρ(i))
                @show  Δi  Δj rA rB li lj
                throw("ah")
            end 
            color_neighbors!(img, R, CartesianIndex(i1, j1), ρ(i))
            err -= lj
            if err < 0
                j1 += sj
                err += li
            end
            i1 += si
            if ρ(i) > 98
                @show   Δi  Δj rA rB li lj
                throw("Don't draw this thick")
            end
        end
    else
        # Radius as a function of horizontal position
        ρ = let rA = rA, rB = rB, nl = lj
            i -> rA + (i - 1) * (rB - rA) / nl
        end
        err = lj ÷ 2
        for i in 1:(lj + 1)
            if isnan(ρ(i))
                @show  Δi  Δj rA rB li lj
                throw("ah")
            end 
    
            spray_neighbors!(img, R, CartesianIndex(i1, j1), ρ(i) )
            err -= li
            if err < 0
                i1 += si
                err += lj
            end
            j1 += sj
            if ρ(i) > 98
                @show   Δi  Δj rA rB li lj
                throw("Don't draw this thick")
            end
        end
    end
    img
end


"""
    draw_two_arrow_glyph!(img, p, v::AbstractVector)

img is the image matrix
p is the location
F is the length of both arrows. If negative, arrows end at p.
v[1:2] specifies both direction and sign. 
    0-π:  Positive sign, arrows outward.   
    π-2π: Negative sign: arrows inward.   
"""
function draw_two_arrow_glyph!(img, p, v::AbstractVector)
    Δj = Int(round(v[1]))
    Δi = Int(round(v[2]))
    if is_quantity_positive(v)
        # First or second quadrant. Positive.
        vector!(img, p, -Δi, Δj)
        vector!(img, p, Δi, -Δj)
    else
        # Third or fourth quadrant. Negative.
        vector!(img, p + CartesianIndex(-Δi, Δj), Δi, -Δj)
        vector!(img, p + CartesianIndex(Δi, -Δj), -Δi, Δj)
    end
end

"""
    is_quantity_positive(v::AbsractVector) --> Bool

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
julia> is_quantity_positive([1.0, 0.0])
true

julia> is_quantity_positive([1.0, -0.0])
false

julia> is_quantity_positive([1.0, -1.0])
false

julia> is_quantity_positive([-1.0, -2.0])
false
```
"""
is_quantity_positive(v::AbstractVector) = !signbit(v[2])


"""
    draw_bidirectional_quantity_glyph!(img, p, K)

Draws two-arrow glyphs on an image at position `p` for principal bi-directional quantities
specified by the matrix `K`.

# Arguments
- `img`: The image matrix to render the glyphs onto.
- `p`: The location of the glyphs in 2D space.
- `K`: A principal bi-directional 4x4 matrix where `K[:,1]` specifies the first principal direction,
  magnitude, and sign, and `K[:,2]` specifies the second. See [`draw_two_arrow_glyph!`](@ref).

# Details
Due to symmetry, only 180° (π) is needed for the direction of each quantity. The additional
information is used to determine if the quantity is negative (arrows toward `p`) or positive
(arrows from `p`).

# Notes
- As an example, 2D principal stress directions are always orthonormal, whereas in `K`, the
  directions need not be orthonormal. If an orthonormal glyph is projected onto a plane that is
  rotated with respect to the stress plane, the projected axes would not be orthonormal. This
  function could be used to draw such a projection.
- This function modifies `img` in-place.

# Examples
```julia
img = zeros(100, 100)  # Example image matrix
p = [50.0, 50.0]       # Center position
K = rand(4, 4)         # Example principal bi-directional matrix
draw_bidirectional_quantity_glyph!(img, p, K)  # Draws two glyphs
```
"""
function draw_bidirectional_quantity_glyph!(img, p, K)
    @assert size(K) == (2, 2)
    draw_two_arrow_glyph!(img, p, view(K, :, 1))
    draw_two_arrow_glyph!(img, p, view(K, :, 2))
end
function draw_bidirectional_quantity_glyph!(img, p, v::AbstractVector)
    @assert size(v, 1) == 2
    draw_two_arrow_glyph!(img, p, v)
end
############################
# [-1, 1] diverging colormap
############################
"""
    mapto_diverging_color(M; palette = RGBA{N0f8}.(colormap("RdBu")))

For M elements in the range -1 to 1. `tanh` is a good limiter.

[-1:1] -> palette[1:length(palette)]

Non-linear at zero (for even-numbered palette). Not optimized.
"""
function mapto_diverging_color(M; palette = RGBA{N0f8}.(colormap("RdBu")))
    @assert minimum(M) >= -1
    @assert maximum(M) <= 1
    n0 = 1
    n1 = length(palette) ÷ 2
    n2 = length(palette)
    fcol = x -> palette[Int(round(x < 0 ? n0 + (x + 1) * (n1 - n0) : n1 + x * (n2 - n1)))]
    map(fcol, M)
end



###########################
# Grayscale
###########################
function to_img(z)
    fsc = scaleminmax00(minimum(z), maximum(z))
    map(z) do ζ
        pix = fsc(ζ)
        GrayA{N0f8}(pix, pix > N0f8(2/255) ? 1N0f8 : 0N0f8)
    end
end

# This relies on global r
function bluesc(z; mi = 0.001, ma = float(r))
    f = scaleminmax00(mi, ma)
    mapto_diverging_color(f.(z))
end
to_green(img) = map(x-> RGBA{N0f8}(0, gray(x), 0, alpha(x)), img)
to_blue(img) = map(x-> RGBA{N0f8}(gray(x), 0, 0, alpha(x)), img)