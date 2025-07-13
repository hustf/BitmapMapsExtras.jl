# This package doesn't implement a color scheme,
# and using this in practice would involve an arbitrary level for flatness anyway,
# and adapted filtering to focus on the right size of categories.
# This is based off 'calculate_and_draw_glyphs.jl/plot_curvature_glyphs!
using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: paint_curvature_types!
using BitmapMapsExtras: N0f8
using SHA # Test dependency
import BitmapMapsExtras.BitmapMaps
using BitmapMapsExtras.BitmapMaps: scaleminmax, RGB, RGBA
##########################################
# COPYPASTA - `t_calculate_and_draw_glyps`
##########################################
if ! @isdefined hashstr
    function hashstr(img)
        iob = IOBuffer()
        show(iob, img)
        bytes2hex(sha1(take!(iob)))
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

#################################
# Curvature types, graphical test
#################################





f = (img, z, grpts) -> paint_curvature_types!(img, z, grpts;
    maxcurv_flat = 0.00002)

# All types
img = grid_fcall_with_background(; f, z = z_ridge_peak_valleys(), Δ = 1)
@test hashstr(img) == "98ed0c914b98ebda87efa0ba3cb58f7b53fed303" 

# All convex (red)
img = grid_fcall_with_background(; f, z = z_paraboloid(; a= 0.6r, b = 0.5r), Δ = 1)
@test hashstr(img) == "8f33c5a74801dbe4bbd9b1e8a3f909fa6507aa91"

# All concave  (green) or flat (gray). Imperfect at edges.
img = grid_fcall_with_background(; f, z = z_ellipsoid(), Δ = 1)
@test hashstr(img) == "a73a45c827d6b4f3036c4bc581e1905e60b93f02"

# All convex-concave, saddle (blue)
img = grid_fcall_with_background(; f, z = z_paraboloid(;a= 0.6r, b = -0.4r), Δ = 1)
@test hashstr(img) == "4f4ecc0d7f08ce57bdb29e814a41760c60a5ff6e"

# A cylinder. Concave.
img = grid_fcall_with_background(; f, z = z_cylinder(1), Δ = 1)
@test hashstr(img) == "9113d7dcda9f5845e3e924419b7e49cea7a2fc0b"

# Lower part of a cylinder. Convex
img = grid_fcall_with_background(; f, z = -z_cylinder(π / 6), Δ = 1)
@test hashstr(img) == "911d663925b9c0b156f82b432c4018c43b9f6212"

