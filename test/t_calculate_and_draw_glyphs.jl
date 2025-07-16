using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: plot_tangent_basis_glyphs, plot_tangent_basis_glyphs!
using BitmapMapsExtras: plot_curvature_glyphs, plot_curvature_glyphs!
using BitmapMapsExtras: plot_𝐧ₚ_glyphs, plot_𝐧ₚ_glyphs!

include("common.jl")

# Tangent basis
img = grid_fcall()
@test hashstr(img) == "866a354a2d14df7cbb971e05b41cdc2a1b63eb09"
img = grid_fcall(;z = z_cylinder_offset(π/6)) 
@test hashstr(img) == "ca900c961cd455d743d67282a5cbced92aa05eab"

# Principal curvatures
f = (z, grpts) -> plot_curvature_glyphs(z, grpts; multglyph = 12000, minglyph = -100)
img = grid_fcall(;f)
@test hashstr(img) == "d36a3d7b6d61bda0cce4b9b7971603708f9aaa26"
img = grid_fcall(;f, z = z_cylinder_offset(π/6)) 
@test hashstr(img) == "9a9c3b21060af3683af2c38f2fefc3bf630856fb"

# Projected normal vector
f = (z, grpts) -> plot_𝐧ₚ_glyphs(z, grpts; multglyph = 100, maxglyph = 100)
img = grid_fcall(; f)
@test hashstr(img) == "26e586df4e59e3860ae43bdcceb0d4b560cde52c" 
img = grid_fcall(; f, z = z_cylinder_offset(π/6)) 
@test hashstr(img) == "4a22a722924cf650c518eecf68a1dd88ef15fbb4"


img = background(z_paraboloid(;a = -0.5r, b = r))
@test hashstr(img) == "7d773c1777716d688b77544844c475ae9cd3c0ca"


# Tangent basis
img = grid_fcall_with_background()
@test hashstr(img) == "77134cde5cb5a067f5a493d4a0d53463d3623650"
img = grid_fcall_with_background(;z = z_cylinder(π/6)) 
@test hashstr(img) == "ed5c36af5a5abb927d0bde50a604841e05244448"
img = grid_fcall_with_background(;z = z_cylinder_offset(π/6)) 
@test hashstr(img) == "01d7f73a5fffbc46f3ba231df9adc61d6ebb85c9"


# Principal curvatures
f = (img, z, grpts) -> plot_curvature_glyphs!(img, z, grpts; directions = 1:2, multglyph = 10000, minglyph = -50, maxglyph = 50)
img = grid_fcall_with_background(; f)
@test let 
    h = hashstr(img) 
    h == "c4b2e2a1c2417376b87ad792c33350ed1baf92eb" || h == "08cd7d8db355aeea69d360018147b25a85a0f29b"
end

# We see acceptable inaccuracy, more close to the edge of the cylinder
img = grid_fcall_with_background(; f, z = z_cylinder_offset(π/6)) 
@test hashstr(img) == "eb11951901678323df607a988e2611250af67557"

# Projected normal vector
f = (img, z, grpts) -> plot_𝐧ₚ_glyphs!(img, z, grpts; multglyph = 100, maxglyph = 100)
img = grid_fcall_with_background(; f)
@test hashstr(img) == "0ccfff2719b6f554a5eced8b8d5671ed51b8b804"
img = grid_fcall_with_background(; f, z = z_cylinder_offset(π/6)) 
@test hashstr(img) == "8e20b76aa47c1579a2eb4830515ca7f9377c46a2"


#######################################################################
# Illustrating ridgeness with glyps (sharpness increases down in image)
#######################################################################
img = background(z_ridge_peak_valleys())
@test hashstr(img) == "3784de2112187f571dc642dfe50899befc54f60a"

# Principal curvatures (the variation is too large for 
# these glyphs to be very useful in a grid)
f = (img, z, grpts) -> plot_curvature_glyphs!(img, z, grpts; 
    multglyph = 50000, minglyph = -50, maxglyph = 50, dashsize = 0)
img = grid_fcall_with_background(; f, z = z_ridge_peak_valleys(), Δ = 50) 
@test hashstr(img) == "105ce8676076a128aa2636453b1bba315f4be2f0" 


# Projected normal vector (less varied, more useful)
f = (img, z, grpts) -> plot_𝐧ₚ_glyphs!(img, z, grpts;
    multglyph = 5000, maxglyph = 50, dashsize = 0)
img = grid_fcall_with_background(; f, z = z_ridge_peak_valleys(), Δ = 50) 
@test hashstr(img) == "69a31c992dee85ed7b8c2b7776797b6af2596e95"
