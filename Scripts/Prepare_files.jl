using Rasters, ArchGDAL

function resample_raster(surface::String, bed::String, output_surface::String, output_bed::String; res=25.0)

    surface = Raster(surface)
    bed = Raster(bed)

    surface = resample(surface; res=res, method=:bilinear)
    bed = resample(bed; to=surface, method=:bilinear)

    Rasters.write(output_surface, surface; force=true)
    Rasters.write(output_bed, bed; force=true)

end

#### Main ####
surface = "I:/doss/x8891/01_projets_agence/paprog/99_Production/JMeugnier/Surf_Bionnassay_2023.tif"
bed = "I:/doss/x8891/01_projets_agence/paprog/99_Production/JMeugnier/Bed_Bionnassay_Farinotti.tif"
new_surface = "T:/RTM/06_ROGP/01_SUIVI_PAPROG/02-ACTIONS_EN_COURS_ONF/2026/4_stage-poches-eau/06_notes/WWFS/Input/Surf_Bionnassay_2023_25m.tif"
new_bed = "T:/RTM/06_ROGP/01_SUIVI_PAPROG/02-ACTIONS_EN_COURS_ONF/2026/4_stage-poches-eau/06_notes/WWFS/Input/Bed_Bionnassay_Farinotti_25m.tif"

resample_raster(surface, bed, new_surface, new_bed; res=25.0)