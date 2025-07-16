using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: plot_streamlines!, 𝐧ₚ!

include("common.jl")


f = (img, z, grpts) -> plot_streamlines!(img, 𝐧ₚ!, z, grpts; tstop = 100000.0, r = 1.2)
img = grid_fcall_with_background(;f, z = z_cylinder(π / 6), Δ = 150)
@test hashstr(img) == "b38eb23d90f753d281354008847eae19aad28231"

img = grid_fcall_with_background(;f, z = z_sphere(), Δ = 50)
@test hashstr(img) == "cfa9e9ad98edcade553e81af25c3d3c5c41e9bb9"

img = grid_fcall_with_background(;f, z = z_paraboloid(a=0.5r, b=-r), Δ = 50)
@test hashstr(img) == "49f37d0f9b1196d8939e61ef6b7c959b81f8bda9"

img = grid_fcall_with_background(;f, z = z_ridge_peak_valleys())
@test hashstr(img) == "70ce31a011f68a98acfd2792aac084af156380eb"

f = (img, z, grpts) -> plot_streamlines!(img, 𝐧ₚ!, z, grpts; tstop = 100000.0, r = 0.2)
img = grid_fcall_with_background(;f, z = z_ridge_peak_valleys(), Δ = 15)
@test hashstr(img) == "ccbe6b3f8843f9a43a6a112208be94406477a211"

#= 
TODO: Improve spray, add parameter 
for r in 0.2:1:5.2
    println(r)
    f = (img, z, grpts) -> plot_streamlines!(img, 𝐧ₚ!, z, grpts; tstop = 100000.0, r)
    img = grid_fcall_with_background(;f, z = z_ridge_peak_valleys(), Δ = 15)
    open_as_temp_in_imgedit(img)
    sleep(2)
end

=#