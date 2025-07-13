# Basis drawing functionality,
# builds on `BitmapMaps.jl/mark_utils.jl`.
# Contains `draw_vector!` and `spray_neighbors!`


"""
    draw_vector!(img, A::CartesianIndex{2}, Δi, Δj)

Draws a line from A to B in the given image, with a max thickness defined by `tol_dist`.
"""
function draw_vector!(img, A::CartesianIndex{2}, Δi::Int, Δj::Int)
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
        spray_neighbors!(img, R, CartesianIndex(i1, j1), rB)
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
            spray_neighbors!(img, R, CartesianIndex(i1, j1), ρ(i))
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
    spray_neighbors!(img::Matrix{GrayA{N0f8}}, R, C::CartesianIndex{2}, r)

Spray the area around a pixel, alpha channel variation is parabolic with distance.

All neighbors of index `C` within distance `r` and bounds `R` get 
value 1 and alpha channel varying from 1 at centre to 0 at r2. Parabolic variation.
"""
function spray_neighbors!(img::Matrix{GrayA{N0f8}}, R, C::CartesianIndex{2}, r)
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