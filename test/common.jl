# For testing graphically
# 
import BitmapMapsExtras.BitmapMaps
using BitmapMapsExtras.BitmapMaps: scaleminmax, RGBA, N0f8
using BitmapMapsExtras:  PALETTE_BACKGROUND, RGB, AbstractGlyphSpec
import BitmapMapsExtras.TestMatrices
using SHA # Test dependency, temporarily added as a direct dependency.
import Random
using Random: MersenneTwister

function hashstr(img)
    iob = IOBuffer()
    show(iob, img)
    bytes2hex(sha1(take!(iob)))
end

# TODO move this into TestMatrices. And make that a separate package?
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


# TODO move this into TestMatrices. And make that a separate package?
"""
    pack_glyphs_with_background(z::Matrix{<:AbstractFloat}, gs::AbstractGlyphSpec; scatterdist = 3.0, seed = MersenneTwister(123))
"""
function pack_glyphs_with_background(z::Matrix{<:AbstractFloat}, gs::AbstractGlyphSpec; scatterdist = 3.0, seed = MersenneTwister(123))
    img = background(z)
    pack_glyphs!(img, z, gs; scatterdist, seed)
end
nothing
