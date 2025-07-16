# For testing graphically
# 
import BitmapMapsExtras.BitmapMaps
using BitmapMapsExtras.BitmapMaps: scaleminmax, RGBA, N0f8
using SHA # Test dependency

if ! @isdefined hashstr
    function hashstr(img)
        iob = IOBuffer()
        show(iob, img)
        bytes2hex(sha1(take!(iob)))
    end
end

if ! @isdefined grid_fcall 
    function grid_fcall(; 
        f = (z, grpts) -> plot_tangent_basis_glyphs(z, grpts), 
        z = z_paraboloid(;a = -0.5r, b = r), Δ::Int = 100)
        # Let's center the grid as far a possible.
        # The centex index is r + 1 in both directions
        lo = (r + 1) + Δ * ceil(Int, -r / Δ)
        hi = (r + 1) + Δ * floor(Int,  r / Δ)
        rng = lo:Δ:hi
        grpts = CartesianIndices((rng, rng))
        # Call it
        f(z, grpts)
    end
end


if ! @isdefined background 
    function background(z; α = 0.6)
        foo = scaleminmax(extrema(z)...)
        fi = x -> RGBA{N0f8}(x, x, x, α)
        img = fi.(foo.(z))
        # Add simple contour lines, too
        Δc = -(-(extrema(z)...)) / 10 # elevation spacing
        wc = Δc / 10         # 'width' of contour lines, in height....
        map!(img, z, img) do zz, pix
            mod(zz, Δc) < wc ? RGBA{N0f8}(0.1, 0.1, 0.1, 1.0) : pix 
        end
    end
end

if ! @isdefined grid_fcall_with_background
    function grid_fcall_with_background(; 
        f = plot_tangent_basis_glyphs!,
        z = z_paraboloid(;a = -0.5r, b = 0.5r),
        Δ::Int = 100)
        # Let's center the grid as far a possible.
        # The centex index is r + 1 in both directions
        lo = (r + 1) + Δ * ceil(Int, -r / Δ)
        hi = (r + 1) + Δ * floor(Int,  r / Δ)
        rng = lo:Δ:hi
        grpts = CartesianIndices((rng, rng))
        img = background(z)
        # Call it
        f(img, z, grpts)
    end
end
