#=
using Rasters, ArchGDAL

function resample_raster(r1_path::String, r2_path::String, output_path_r1::String, output_path_r2::String)
    r1 = Raster(r1_path)
    r2  = Raster(r2_path)
    r1 = resample(r1; to=r2, method=:bilinear)
    Rasters.write(output_path_r1, r1, force = true)
    Rasters.write(output_path_r2, r2, force = true)
    println("$(size(r1))")
    println("$(size(r2))")  
end


surface = "I:/doss/x8891/01_projets_agence/paprog/99_Production/JMeugnier/Surf_MdG_2025.tif"
bed = "I:/doss/x8891/01_projets_agence/paprog/99_Production/JMeugnier/Bed_MdG_Farinotti.tif"
new_surface = "T:/RTM/06_ROGP/01_SUIVI_PAPROG/02-ACTIONS_EN_COURS_ONF/2026/4_stage-poches-eau/06_notes/WWFS/Input/Surf_MdG_2025_25m.tif"
new_bed = "T:/RTM/06_ROGP/01_SUIVI_PAPROG/02-ACTIONS_EN_COURS_ONF/2026/4_stage-poches-eau/06_notes/WWFS/Input/Bed_MdG_Farinotti.tif"
 

resample_raster(surface, bed, new_bed, new_surface)
=#

using Rasters, ArchGDAL

function resample_raster(surface_path::String, bed_path::String, output_surface::String, output_bed::String; res=25.0)

    surface = Raster(surface_path)
    bed = Raster(bed_path)

    # Surface sur une grille carrée de 25 m
    surface_25 = resample(surface; res=res, method=:bilinear)

    # Bedrock aligné exactement sur cette grille
    bed_25 = resample(bed; to=surface_25, method=:bilinear)

    Rasters.write(output_surface, surface_25; force=true)
    Rasters.write(output_bed, bed_25; force=true)

    println("Surface : ", size(surface_25))
    println("Bedrock : ", size(bed_25))

    x = Vector(dims(surface_25, X).val)
    y = Vector(dims(surface_25, Y).val)

    println("dx = ", x[2] - x[1])
    println("dy = ", y[2] - y[1])
end


surface = "I:/doss/x8891/01_projets_agence/paprog/99_Production/JMeugnier/Surf_MdG_2025.tif"
bed     = "I:/doss/x8891/01_projets_agence/paprog/99_Production/JMeugnier/Bed_MdG_Farinotti.tif"
new_surface = "T:/RTM/06_ROGP/01_SUIVI_PAPROG/02-ACTIONS_EN_COURS_ONF/2026/4_stage-poches-eau/06_notes/WWFS/Input/Surf_MdG_2025_25m.tif"
new_bed = "T:/RTM/06_ROGP/01_SUIVI_PAPROG/02-ACTIONS_EN_COURS_ONF/2026/4_stage-poches-eau/06_notes/WWFS/Input/Bed_MdG_Farinotti_25m.tif"

resample_raster(surface, bed, new_surface, new_bed; res=25.0)