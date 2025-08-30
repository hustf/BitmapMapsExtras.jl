# For testing graphically
# 
import BitmapMapsExtras.BitmapMaps
using BitmapMapsExtras.BitmapMaps: scaleminmax, RGBA, N0f8
using BitmapMapsExtras:  PALETTE_BACKGROUND, Lab, RGB
import BitmapMapsExtras.TestMatrices
using SHA # Test dependency, temporarily added as a direct dependency.

function hashstr(img)
    iob = IOBuffer()
    show(iob, img)
    bytes2hex(sha1(take!(iob)))
end


function background(z; α = 1.0)
    foo = scaleminmax(extrema(z)...)
    scaled_z = foo.(z)
    img = map(scaled_z) do z
        RGBA{N0f8}(get(PALETTE_BACKGROUND, z), α)
    end
    # Add simple contour lines, too
    Δc = -(-(extrema(z)...)) / 10 # elevation spacing
    wc = Δc / 10         # 'width' of contour lines, in height....
    map!(img, z, img) do zz, pix
        mod(zz, Δc) < wc ? RGBA{N0f8}(0.1, 0.1, 0.1, 1.0) : pix 
    end
end


"""
    grid_indices((h, w)::Tuple{Int64, Int64}; Δ::Int = 100,
        offset::Tuple{Int64, Int64} = (0, 0))
        h, w = size(img)
    ---> CartesianIndices{2, Tuple{StepRange{Int64, Int64}, StepRange{Int64, Int64}}}
"""
function grid_indices((h, w)::Tuple{Int64, Int64}; Δ::Int = 100,
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

"""
    grid_fcall!(f, vargs; 
        z = z_paraboloid(),
        Δ::Int = 100,
        centreoffset::Int = 0) where T <: AbstractVector{<:Tuple{Tuple{Vararg}, NamedTuple}}

`vargs` is mapped to f as arguments following 'img', 'z' and 'grpts':

    f(img, z, grpts, positional arguments; keyword arguments)

See `convert_to_list_of_arguments` for 'vargs'.
"""
function grid_fcall!(img, f, vargs; 
    z = z_paraboloid(),
    Δ::Int = 100,
    offset::Tuple{Int64, Int64} = (0, 0))
    #
    grpts = grid_indices(size(img); Δ, offset)
    fcall(f, vargs, img, z, grpts)
    img
end


"""
    grid_fcall_with_background(f, vargs; 
        z = z_paraboloid(),
        Δ::Int = 100,
        offset::Tuple{Int64, Int64} = (0, 0))

`vargs` is mapped to f as arguments following the ones produced by grid: 'img', 'z', 'grpts':

    f(img, z, grpts, positional arguments; keyword arguments)

See `convert_to_list_of_arguments` for 'vargs'.
"""
function grid_fcall_with_background(f, vargs; 
    z = z_paraboloid(),
    Δ::Int = 100,
    offset::Tuple{Int64, Int64} = (0, 0))
    #
    grpts = grid_indices(size(z); Δ, offset)
    fcall_with_background(f, vargs, z, grpts)
end

"""
    fcall_with_background(f, vargs, z, grpts)

`vargs` is mapped to f as arguments following: img, z, grpts.

Each call to f: f(img, z, grpts, positional arguments; keyword arguments)

See `convert_to_list_of_arguments` for 'vargs'.
"""
function fcall_with_background(f, vargs, z, grpts)
    img = background(z)
    fcall(f, vargs, img, z, grpts)
end


"""
    fcall(f!, vargs::T, img, z, grpts)
    fcall(f!, vargs, img, z, grpts)

Map an argument list (img, z, grpts, positional arguments; keyword arguments) to f!

See `convert_to_list_of_arguments` for 'vargs'.
"""
function fcall(f!, vargs::T, img, z, grpts)  where T <: AbstractVector{<:Tuple{Tuple{Vararg}, NamedTuple}}
    if !isempty(vargs)
        for (pos, kwargs) in vargs
            if ! (pos isa Tuple{Nothing})
                #println("\nCalling $f!(img, z, grpts, $(pos...); kwargs = ($kwargs...))")
                f!(img, z, grpts, pos...; kwargs...)
            else
                # Keyword arguments only
                #println("\nCalling $f!(img, z, grpts; kwargs = ($kwargs...))")
                f!(img, z, grpts; kwargs...)
            end
        end
    else
        #println("\nCalling `$f!(img, z, grpts)`")
        kwargs = (;)
        f!(img, z, grpts; kwargs...)
    end
    img
end
# Catchall for vargs
fcall(f!, vargs, img, z, grpts) =  fcall(f!, convert_to_list_of_arguments(vargs), img, z, grpts)

"""
    convert_to_list_of_arguments(args)

1) Function with no arguments, no keyword arguments. This is just the minimum case 2).
    [((f,), (;))]

2) Functions with no arguments, no keyword arguments
    [((f1,), (;)), 
    ((f2,), (;)) ]

3) Functions with arguments (and keyword arguments)
    [((f, 1), (; kwd1 = value1, kwd2 = value2)), 
    ((f, 2), (; kwd1 = value1, kwd2 = value2))]

Types for args that are converted to the above as default:    

4) Single function without arguments
    f

5) Single function with positional arguments
    (f, 1, 2)

6) Single function with keyword arguments
    ((f,), (;kwd = arg))

7) Single function with positional and keyword arguments
    ((f, 1,), (;kwd = arg))
"""
function convert_to_list_of_arguments(args)
    # Case 1-3: Already in canonical form (Vector of Tuples)
    if args isa AbstractVector{<:Tuple{Tuple{Vararg}, NamedTuple}}
        return args
    end
    # Case 4: Single function without arguments
    if args isa Function
        return [((args,), (;))]
    end
    # Case 5: Single function with positional arguments
    if args isa Tuple && !isempty(args) && args[1] isa Function
        return [((args...), (;))]
    end
    # Case 6 & 7: Single function with/without positional args and keyword args
    if args isa Tuple{<:Tuple, <:NamedTuple}
        pos_args, kw_args = args
        if !isempty(pos_args) && pos_args[1] isa Function
            return [(pos_args, kw_args)]
        end
    end
    # Case 8: Nothing
    if isnothing(args)
        return Vector{Tuple{Tuple{Function, Vararg}, NamedTuple}}()
    end
    # Case 9: Keywords
    if args isa NamedTuple
        return Vector{Tuple{Tuple{Nothing}, NamedTuple}}([((nothing,), args)])
    end
    # Case 10: Vector of Nothing
    if args isa Vector{Nothing}
        return [((nothing,), (;)) for i = 1:length(args)]
    end
    # Case 11: Vector of Keywords
    if args isa Vector{<:NamedTuple}
        return [((nothing,), kws) for kws in args]
    end
    # Unrecognized format
    @warn("Canonical types for 'args': 

    1) Function with no arguments, no keyword arguments. This is just the minimum case 2).
        [((f,), (;))]

    2) Functions with no arguments, no keyword arguments
       [((f1,), (;)), 
        ((f2,), (;)) ]
    
    3) Functions with arguments (and keyword arguments)
       [((f, 1), (; kwd1 = value1, kwd2 = value2)), 
        ((f, 2), (; kwd1 = value1, kwd2 = value2))]

    Types for args that are converted to the above as default:    
    
    4) Single function without arguments
        f

    5) Single function with positional arguments
        (f, 1, 2)

    6) Single function with keyword arguments
        ((f,), (;kwd = arg))

    7) Single function with positional and keyword arguments
        ((f, 1,), (;kwd = arg))
    ")
    #@show typeof(args)
    #@show args
    throw(ErrorException("Unrecognized format for `args`"))
end
nothing