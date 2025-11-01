# Contains 
# - indices_scattered
# - indices_on_grid


"""
    indices_on_grid((h, w)::Tuple{Int64, Int64}; Δ::Int = 100,
        offset::Tuple{Int64, Int64} = (0, 0))
    indices_on_grid(fxy; Δ = 100, offset = (0, 0))
    ---> CartesianIndices{2, Tuple{StepRange{Int64, Int64}, StepRange{Int64, Int64}}}

The first argument can be a tuple, an image, an AbstractXYFunctor or any x for which size(x) works.
"""
function indices_on_grid((h, w)::Tuple{Int64, Int64}; Δ::Int = 100,
    offset::Tuple{Int64, Int64} = (0, 0))
    # Center
    ci = h ÷ 2 + 1 + offset[1]
    cj = w ÷ 2 + 1 + offset[2]
    # Let's center the grid as far a possible.
    loi = ci + Δ * ceil(Int, -(h ÷ 2) / Δ)
    hii = ci + Δ * floor(Int,  (h ÷ 2) / Δ)
    rngi = loi:Δ:hii
    loj = cj + Δ * ceil(Int, -(w ÷ 2) / Δ)
    hij = cj + Δ * floor(Int, (w ÷ 2) / Δ)
    rngj = loj:Δ:hij
    grpts = CartesianIndices((rngi, rngj))
end
indices_on_grid(fxy; Δ = 100, offset = (0, 0)) = indices_on_grid(size(fxy); Δ, offset )

"""
    indices_scattered(R::CartesianIndices; scatterdist = 3.0, seed = MersenneTwister(123))
    --> Vector{CartesianIndex}
    indices_scattered(fij::AbstractIJFunctor; scatterdist = 3.0, seed = MersenneTwister(123))
    indices_scattered(z::Matrix{<:AbstractFloat}, gs::GSTangentBasis;
         scatterdist = 3.0, seed = MersenneTwister(123))

Return a vector of random CartesianIndex scattered within `R` such that the
mean nearest-neighbour distance is approximately `scatterdist`. Refer Poisson disk sampling.

# Arguments

- `R` defines the possible values of each returned index. 
- `scatterdist` (default keyword value 10.0) is the mean nearest-neighbour distance between points. The 
    resulting set of points will vary slightly from this target.
- `fz`: Restricts indices to the domain of the functor `fz`
"""
function indices_scattered(R::CartesianIndices; scatterdist = 3.0, seed = MersenneTwister(123))
    area = length(R)
    # Estimate number of points, within bounds.
    n = max(1, min(round(Int, area / (4 * scatterdist^2)), area))
    # Draw random unique indices
    chosen = randperm(seed, area)[1:n]
    [R[i] for i in chosen]
end
function indices_scattered(fij::AbstractIJFunctor;
     scatterdist = 3.0, seed = MersenneTwister(123))
    #
    R = CartesianIndices(fij.z)
    Ω = fij.Ω
    minj = R[1][2] - Ω[1][2]
    maxi = R[end][1] - Ω[end][2]
    maxj = R[end][2] - Ω[end][2]
    mini = R[1][1] - Ω[1][2]
    Ri = CartesianIndices((mini:maxi, minj:maxj))
    indices_scattered(Ri; scatterdist, seed)
end
function indices_scattered(z::Matrix{<:AbstractFloat}, gs::GSTangentBasis;
     scatterdist = 3.0, seed = MersenneTwister(123))
    #
    R = CartesianIndices(z)
    Δ = gs.halfsize
    minj = R[1][2] + Δ
    maxi = R[end][1] - Δ
    maxj = R[end][2] - Δ
    mini = R[1][1] + Δ
    Ri = CartesianIndices((mini:maxi, minj:maxj))
    indices_scattered(Ri; scatterdist, seed)
end
