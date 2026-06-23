using WhereTheWaterFlows, GLMakie, Rasters, ArchGDAL
using ImageMorphology, GeoInterface
using Statistics

const WWF = WhereTheWaterFlows
const WWFS = WhereTheWaterFlows.Subglacially

# ── Chemins ──────────────────────────────────────────────────────────────────

surface_path = "T:/RTM/06_ROGP/01_SUIVI_PAPROG/02-ACTIONS_EN_COURS_ONF/2026/4_stage-poches-eau/06_notes/WWFS/Input/Surf_MdG_2025_25m.tif"
bed_path = "T:/RTM/06_ROGP/01_SUIVI_PAPROG/02-ACTIONS_EN_COURS_ONF/2026/4_stage-poches-eau/06_notes/WWFS/Input/Bed_MdG_Gilbert_25m.tif"
outline_path = "I:/doss/x8891/01_projets_agence/paprog/99_Production/JMeugnier/Outlines_MdG.gpkg"
depression_path = "I:/doss/x8891/01_projets_agence/paprog/99_Production/JMeugnier/Depressions_MdG_2025.gpkg"
outdir = "T:/RTM/06_ROGP/01_SUIVI_PAPROG/02-ACTIONS_EN_COURS_ONF/2026/4_stage-poches-eau/06_notes/WWFS/Output"

# ── Lecture des rasters ──────────────────────────────────────────────────────

surface_r = Raster(surface_path)
bed_r     = Raster(bed_path)

x = Vector{Float32}(dims(surface_r, X).val)
y = Vector{Float32}(dims(surface_r, Y).val)

print(x[2]-x[1])
print(y[2]-y[1])

surface_r = replace_missing(surface_r, NaN32)
bed_r     = replace_missing(bed_r, NaN32)

surface = Float32.(Matrix(surface_r))
bed     = Float32.(Matrix(bed_r))

dx = abs(Float32(x[2] - x[1]))

print(size(surface))
print(size(bed))

@assert size(surface) == size(bed) "Les deux rasters doivent avoir la même taille !"

println("Taille : ", size(surface))
println("dx = ", dx, " m")


# ── Lissage de la surface ────────────────────────────────────────────────────
function smooth_surface_v2(x, y, surfdem, beddem, icethicknesses, mask=surfdem.>=beddem; minwindow=0)
    dx = x[2]-x[1]
    @assert abs(y[2]-y[1])==dx
    thickness = surfdem .- beddem
    thickness[isnan.(thickness)] .= 0
    window = round.(Int, icethicknesses .* thickness ./ dx)
    window[window.<minwindow] .= minwindow
    surface_smooth = WWFS.boxcar(surfdem, window, mask)
    return surface_smooth
end

smooth_frac = 0.1  # fraction de l'épaisseur locale de glace
mask_ice = isfinite.(surface) .& isfinite.(bed) .& (surface .>= bed)
surface = smooth_surface_v2(x, y, surface, bed, smooth_frac; minwindow = 1)
surface[.!mask_ice] .= NaN


# ── Calcul WWF ───────────────────────────────────────────────────────────────

@time out = WWFS.waterflows_subglacial(surface, bed, dx; gamma = WWFS.GAMMA, drain_pits = true, avoid_sc   = false)

# ── Variables WWF ────────────────────────────────────────────────────────────

phi    = out.routing.phi
flow   = out.routing.area.total
hwater = out.lakes.depth_free_surface

# ── Charge hydraulique ───────────────────────────────────────────────────────

head = copy(phi)

valid_head = head[isfinite.(head)]

head_step =10
head_min = floor(minimum(valid_head) / head_step) * head_step
head_max = ceil(maximum(valid_head) / head_step) * head_step
head_levels = head_min:head_step:head_max

# ── Chenaux principaux ───────────────────────────────────────────────────────

valid_flow = flow[isfinite.(flow)]

flow_quantile = 0.98
thr = quantile(valid_flow, flow_quantile)

mainflow = flow .> thr
mainflow_skeleton = thinning(mainflow)

# ── Poches d'eau ─────────────────────────────────────────────────────────────

water_threshold = 0.05
volume_threshold = 1000.0

lake_mask = hwater .> water_threshold
labels = label_components(lake_mask)

bigmask = falses(size(labels))

for lab in unique(labels)
    lab == 0 && continue

    idx = labels .== lab
    volume = sum(hwater[idx]) * dx^2

    if volume > volume_threshold
        bigmask .|= idx
    end
end

smallmask = lake_mask .& .!bigmask

# ── Rasters pour affichage / QGIS ────────────────────────────────────────────

hwater_bigonly = copy(hwater)
hwater_bigonly[.!bigmask] .= NaN
hwater_bigonly[.!isfinite.(surface)] .= NaN

hwater_smallonly = copy(hwater)
hwater_smallonly[.!smallmask] .= NaN
hwater_smallonly[.!isfinite.(surface)] .= NaN

hwater_all = copy(hwater)
hwater_all[.!lake_mask] .= NaN
hwater_all[.!isfinite.(surface)] .= NaN

# ── Fonction export GeoTIFF ──────────────────────────────────────────────────

function write_tif(path, data, ref_raster)
    r = Raster(data; dims=dims(ref_raster))
    Rasters.write(path, r; force=true)
    println("Écrit : ", path)
end

# Afficher vecteurs

function plot_polygon_layer!(ax, vector_path; color=(:limegreen, 0.35), strokecolor=:limegreen, strokewidth=1.5)

    ArchGDAL.read(vector_path) do ds
        layer = ArchGDAL.getlayer(ds, 0)

        for feature in layer
            geom = ArchGDAL.getgeom(feature)
            geomname = ArchGDAL.geomname(geom)

            if geomname == "POLYGON"
                ring = ArchGDAL.getgeom(geom, 0)

                pts = Point2f[]
                for i in 0:(ArchGDAL.ngeom(ring)-1)
                    p = ArchGDAL.getpoint(ring, i)
                    push!(pts, Point2f(p[1], p[2]))
                end

                p = poly!(ax, pts; color=color, strokecolor=strokecolor, strokewidth=strokewidth)
                
            elseif geomname == "MULTIPOLYGON"
                for j in 0:(ArchGDAL.ngeom(geom)-1)
                    polygeom = ArchGDAL.getgeom(geom, j)
                    ring = ArchGDAL.getgeom(polygeom, 0)

                    pts = Point2f[]
                    for i in 0:(ArchGDAL.ngeom(ring)-1)
                        p = ArchGDAL.getpoint(ring, i)
                        push!(pts, Point2f(p[1], p[2]))
                    end

                    poly!(ax, pts; color=color, strokecolor=strokecolor, strokewidth=strokewidth)
                end
            end
        end
    end
end

# ── Export TIFF ─────────────────────────────────────────────────

write_tif(joinpath(outdir, "WWFS_head_hydraulic.tif"),
          Float32.(head),
          surface_r)

write_tif(joinpath(outdir, "WWFS_hwater_pockets.tif"),
          Float32.(hwater_all),
          surface_r)

write_tif(joinpath(outdir, "WWFS_hwater_big_pockets_gt1000m3.tif"),
          Float32.(hwater_bigonly),
          surface_r)

write_tif(joinpath(outdir, "WWFS_hwater_small_pockets_lt1000m3.tif"),
          Float32.(hwater_smallonly),
          surface_r)

write_tif(joinpath(outdir, "WWFS_flow_accumulation.tif"),
          Float32.(flow),
          surface_r)

# ── Figure ─────────────────────────────────────────────────────────────

valid_water = hwater[lake_mask .& isfinite.(hwater)]
water_max = isempty(valid_water) ? 9.0 : ceil(maximum(valid_water))

fig = Figure(size = (1600, 1000))

ax = Axis(fig[1, 1], xlabel = "X (m)", ylabel = "Y (m)", title = "Hydrologie sous-glaciaire")

# ── Contours glaciers ────────────────────────────────────────────────────────

ArchGDAL.read(outline_path) do dataset
    layer = ArchGDAL.getlayer(dataset, 0)

    for feature in layer
        geom = ArchGDAL.getgeom(feature)

        # Cas polygone simple
        if ArchGDAL.geomname(geom) == "POLYGON"
            ring = ArchGDAL.getgeom(geom, 0)
            coords = ArchGDAL.getpoint.(Ref(ring), 0:(ArchGDAL.ngeom(ring)-1))

            xs = Float64[p[1] for p in coords]
            ys = Float64[p[2] for p in coords]

            lines!(ax, xs, ys; color=:black, linewidth=2)
        end

        # Cas multipolygone
        if ArchGDAL.geomname(geom) == "MULTIPOLYGON"
            for i in 0:(ArchGDAL.ngeom(geom)-1)
                poly = ArchGDAL.getgeom(geom, i)
                ring = ArchGDAL.getgeom(poly, 0)
                coords = ArchGDAL.getpoint.(Ref(ring), 0:(ArchGDAL.ngeom(ring)-1))

                xs = Float64[p[1] for p in coords]
                ys = Float64[p[2] for p in coords]

                lines!(ax, xs, ys; color=:skyblue, linewidth=2)
            end
        end
    end
end

plot_polygon_layer!(ax,depression_path;color=RGBAf(0,0,0,0),strokecolor=:green,strokewidth=0.75)

hm = heatmap!(ax,x,y,hwater_all;colormap = Reverse(:viridis),colorrange = (0, water_max),nan_color = RGBAf(0, 0, 0, 0),overdraw = true)

contour!(ax,x,y,head;levels = head_levels,color = (:grey, 0.9),linewidth = 1.2,labels = true,labelsize = 14,labelcolor = :black)
contour!(ax, x, y, Float64.(bigmask); levels=[0.5], color=:red, linewidth=2, overdraw = true)
contour!(ax,x,y,Float64.(smallmask);levels = [0.5],color = :purple,linewidth = 1.5,overdraw = true)
heatmap!(ax,x,y,Float64.(mainflow_skeleton);colormap = [:transparent, :blue],colorrange = (0, 1),interpolate = false,nan_color = RGBAf(0, 0, 0, 0),overdraw = true)

# ── Éléments pour la légende ────────────────────────────────────────────────

leg_glacier = LineElement(color=:skyblue, linewidth=2)
leg_depression = LineElement(color=:green, linewidth=1.5)
leg_channel = LineElement(color=:blue, linewidth=2)
leg_pocket = LineElement(color=:red, linewidth=3)
leg_head = LineElement(color=:grey, linewidth=1.2)

Legend(
    fig[1, 3],
    [leg_glacier, leg_depression, leg_channel, leg_pocket, leg_head],
    ["Contours du glacier", "Dépressions de surface", "Chenaux principaux", "Poches d'eau > 1000 m³", "Potentiel hydraulique"],
    orientation = :vertical,
    framevisible = true
)

Colorbar(fig[1, 2], hm, label = "Hauteur de poche d'eau (m)")

DataAspect()
display(fig)

save(joinpath(outdir, "WWFS_carte_hydrologie_sous_glaciaire_MdG.png"), fig)

println("Exports terminés.")