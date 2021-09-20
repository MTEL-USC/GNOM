#================================================
dust regions εNd values
================================================#
include("../plots_setup_Nd.jl")


# plot the εNd of sources
fig = Figure(; outer_padding, resolution = (1000, 600), backgroundcolor=:white)
ax = fig[1,1] = Axis(fig, backgroundcolor=:gray20)

# Polygons of dust regions
myrectangle(lon_tuple, lat_tuple) = [Point(lon_tuple[i], lat_tuple[j]) for (i,j) in zip([1,1,2,2,1], [1,2,2,1,1])]
POLYGON_DUST_REGIONS = Dict(
    # From email by Jasper Kok:
    #The coordinates of the nine source regions are:
    #(1) western North Africa: 20⁰ W – 7.5⁰ E; 18⁰ N – 37.5⁰ N
    :NWAf => myrectangle((-20, 7.5), (18, 37.5)),
    #(2) eastern North Africa: 7.5E – 35E; 18N – 37.5N
    :NEAf => myrectangle((7.5, 35), (18, 37.5)),
    #(3) the Sahel: 20W – 35E; 0N – 18N
    :Sahel => myrectangle((-20, 35), (0, 18)),
    #(4) the Middle East & Central Asia: 35⁰ E – 75⁰ E for 0⁰ N – 35⁰ N, and 35⁰ E – 70⁰ E for 35⁰ N – 50⁰ N
    :MECA => [Point(35, 0), Point(75, 0), Point(75, 35), Point(70, 35), Point(70, 50), Point(35, 50), Point(35, 0)],
    #(5) East Asia: 70⁰ E – 120⁰ E; 35⁰ N – 50⁰ N
    :EAsia => myrectangle((70, 120), (35, 50)),
    #(6) North America: 130⁰ W – 80⁰ W; 20⁰ N – 45⁰ N
    :NAm => myrectangle((-130, -80), (20, 45)),
    #(7) Australia: 110⁰ E – 160⁰ E; 10⁰ S – 40⁰ S
    :Aus => myrectangle((110, 160), (-40, -10)),
    #(8) South America: 80⁰ W – 20⁰ W; 0⁰ S – 60⁰ S
    :SAm => myrectangle((-80, -20), (-60, 0)),
    #(9) Southern Africa: 0⁰ E – 40⁰ E; 0⁰ S – 40⁰ S
    :SAf => myrectangle((0, 40), (-40, 0))
)
# eNd values
myεNd(r) = tp_opt[!,:Value][findfirst(tp_opt[!,:Symbol] .== Symbol("ε_", r, "_dust"))]
mypoly(r) = POLYGON_DUST_REGIONS[r]
mycenter(r) = mean(extrema(x -> x[1], mypoly(r))), mean(extrema(x -> x[2], mypoly(r)))


# Plot polygons of land
mapit!(ax, cl0, mypolys(cl0), color=:gray50)
# Plot polygons of dust regions
for r in keys(AEOL_Koketal)
    poly!(ax, mypoly(r), color=get(εcmap, myεNd(r), εclims))
    text!(ax, string(r, "\n", round(myεNd(r), digits=1)), position = mycenter(r), align = (:center, :center), textsize=15)
end
# polygons of land again (no fill)
mapit!(ax, cl0, mypolys(cl0), color=:transparent, strokecolor=RGBA(0,0,0,0.25), strokewidth=1)
# Better ticks and labels
mylatlons!(ax, latticks30, lonticks30)
#tightlimits!(ax)

cbar = fig[end+1, 1] = Colorbar(fig, colormap=εcmap, label="εNd (‱)", vertical=false, ticks=range(εclims..., step=5), flipaxis=false, ticklabelalign=(:center, :top))
cbar.limits=εclims
cbar.width = Relative(3/4)
cbar.height = 30
cbar.tellheight = true


#save(joinpath(output_path, "εNd_dust_regions.$EXT"), scene)
if use_GLMakie
    display(fig) # show the output wiht GLMakie
else
    save(joinpath(archive_path, "εNd_dust_region_$(lastcommit)_run$(run_num).pdf"), fig)
    nothing # just so that no output is spat out
end


