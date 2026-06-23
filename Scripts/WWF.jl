
using WhereTheWaterFlows, Makie, GLMakie, ArchGDAL
const AG = ArchGDAL

file_path = "I:/doss/x8891/01_projets_agence/paprog/99_Production/JMeugnier/MNT_Nord_MB_5m_2018.tif"

dataset = ArchGDAL.read(file_path)
dem_raw = ArchGDAL.read(dataset, 1)
dem = Float64.(dem_raw)
dem[dem .< -500] .= NaN

rows, cols = size(dem)
x = 1.0:Float64(rows)
y = 1.0:Float64(cols)

area, slen, dir, nout, nin, sinks, pits, c, bnds = waterflows(reverse(dem, dims = 2))

plt_area(x, y, area, pits)

#plt_it(x, y, dem)

#plt_catchments(x, y, c)

#plt_dir(x, y, dir)

#plt_sinks(x, y, sinks)

#plt_lakedepth(x, y, dem, dir, sinks)

#heatmap(x, y, slen)

#step = 20
#x_sub = x[1:step:end]
#y_sub = y[1:step:end]
#dir_sub = dir[1:step:end, 1:step:end]
#plt_dir(x_sub, y_sub, dir_sub)

#demf = fill_dem(dem, sinks, dir)
#fig = Figure()
#ax = Axis(fig[1,1])
#hm = heatmap!(ax, x, y, reverse(demf .- dem, dims = 2), colormap = :viridis, colorrange = (0,0.2))
#Colorbar(fig[1,2], hm)
#fig
