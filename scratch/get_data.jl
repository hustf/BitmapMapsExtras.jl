using BitmapMaps

# We rely on existing data from following the BitmapDoc walkthrough
pth_data = joinpath(homedir(), "BitmapMaps", "Stetind")
@assert ENV["BMM_CONFDIR"] == pth_data
pth = joinpath(homedir(), "BitmapMaps", "Stetind", "VectorCalculus")
@assert isdir(pth)
# Focus on Prestinden.  
southwest_corner = (565407, 7558496)
# Width and height in utm metre
w, h = ( 568822,  7561011) .- southwest_corner
# Width of sheet on paper in metre (we won't change that)
wp = define_builder().sheet_width_mm / 1000
# Density to fit on one sheet with cell_to_utm factor = 1
density_pt_m⁻¹ = Int(round(w / wp))
# Height of sheet on paper (accomodate aspect ratio)
sheet_height_mm = Int(round(1000 * h / density_pt_m⁻¹ ))
# SheetMatrixBuilder and SheetBuilder
smb = define_builder(; pth, density_pt_m⁻¹, southwest_corner, sheet_height_mm)
sb = smb[1]
# More feedback
ENV["JULIA_DEBUG"] = "BitmapMaps"
# Try to make this map
run_bitmapmap_pipeline(smb)
# That exited quickly since data was not found. Based on the 'info' message:
fofo = full_folder_path(sb)
@assert isdir(fofo)
copy_relevant_tifs_to_folder(pth_data, fofo)
run_bitmapmap_pipeline(smb)
# As in the walkthrough, we change the sun azimuth angle, drop snow in steep places,
#  and delete existing toporelief.png

using BitmapMaps: func_directional_pallette, func_reflection_coefficient,
    topo_relief, cell_to_utm_factor, full_folder_path, TOPORELIEF_FNAM,
    lambert_shade, linterp
function my_steep_topo_relief(sb::SheetBuilder; 
    sun_deg = 156, 
    snow_lim_deg = 45, 
    exponent_factor = 15.0)
    #
    ffna_source = joinpath(full_folder_path(sb), TOPORELIEF_FNAM)
    if isfile(ffna_source)
        rm(ffna_source; force = true)
    end
    f_hypso = func_directional_pallette()
    f_reflect = my_func_refl(; sun_deg, snow_lim_deg, exponent_factor)
    topo_relief(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb), f_hypso, f_reflect)
end

function my_func_refl(; sun_deg = 156, 
    light_elev_deg = [9, 30, 30, 30],
    snow_lim_deg = 45, 
    exponent_factor = 15.0)
    #
    # Direction no: From sun, opposite sun, one side 60°, other side 60 °
    Δazim_deg = [0, -180 - 60, - 180, - 180 + 60]
    light_azim_deg = sun_deg .+ Δazim_deg
    light_azim = light_azim_deg .* (π / 180)
    light_elev = light_elev_deg .* (π / 180)
    #
    # Convert light azimuth and light elevation to x-y-z light unit vectors.
    # Azimuth of 0     <=> Light from north, vector points north
    # Azimuth of π / 2 <=> Light from east, vector points east
    # Azimuth of π     <=> Light from south, vector points south
    # Elevation 0      <=> Light from horizon
    # Elevation π / 2  <=> Light from above
    #
    # Unit vector components of all the light sources
    l_ew = cos.(light_elev) .* sin.(light_azim)
    l_sn = cos.(light_elev) .* cos.(light_azim)
    l_up = sin.(light_elev)
    # Save some time by doing this part just once
    cos_lim = cos(snow_lim_deg * (π / 180))
    # Closure on the light source vectors and snow limit
    f = let  l_ew = l_ew, l_sn = l_sn, l_up = l_up, cos_lim = cos_lim
        (dno, z, n_ew, n_sn, n_up) -> begin
            @assert 1 <= dno <= 4
            # The dot product of light and surface normal is the fraction of light
            # reflected towards the observer. This is 'lambert shade'.
            #
            # We modify the lambert shading where there is snow.
            # Lambert_reflection^shade_exponent 
            # naively does not exceed 1.0. 
            # Thrust, but check.
            min(1.0f0, convert(Float32, 
                lambert_shade(n_ew, n_sn, n_up, l_ew[dno], l_sn[dno], l_up[dno]) ^
                my_shade_exponent(z, dno, n_up, cos_lim, exponent_factor)))
        end
    end
end

function my_shade_exponent(z, direction_no, n_up, cos_lim, exponent_factor)
    z1 = 400
    z2 = 500
    y1 = 1.2
    y2 = direction_no == 3 ? y1 : 0.4
    if z <= z1
        ex = y1
    elseif z < z2
        ex = linterp(y1, y2, z1, z2, z)
    else
        ex = y2
    end
    # Here comes the part about removing snow in steep places
    if z > 450
        # This is a faster way to say:
        # acos(n_up) > snow_lim_deg * (π / 180)
        if n_up < cos_lim
            # A very high exponent makes light reflect in one direction only,
            # in effect making the surface dark.    
            ex *= exponent_factor
        end
    end
    ex
end

my_steep_topo_relief(smb[1])
run_bitmapmap_pipeline(smb)
# Load images for background
using BitmapMaps: COMPOSITE_FNAM, TOPORELIEF_FNAM, CONTOUR_FNAM, RIDGE_FNAM, load
img_topo = load(joinpath(fofo, TOPORELIEF_FNAM))
img_cont = load(joinpath(fofo, CONTOUR_FNAM))
img_ridge = load(joinpath(fofo, RIDGE_FNAM))
img_comp = load(joinpath(fofo, COMPOSITE_FNAM))


