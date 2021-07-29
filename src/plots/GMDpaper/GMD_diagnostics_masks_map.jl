# Compute the diagnostics only once, unless rediagnose=true
(!isdefined(Main, :fDNd_wtags) || rediagnose) && include("GMD_diagnostics_setup.jl")

#=
function plot_diagnosis_mask_map!(fig)
    masks_map = permutedims(verticalmean(sum(i*m for (i,m) in enumerate(masks)), grd))
    ax = fig[1,1] = Axis(fig, backgroundcolor=seafloor_color)
    mapit!(ax, clon, mypolys(clon), color=continent_color)
    hm = heatmap!(ax, sclons, lats, view(masks_map, ilon, :), colormap=maskscmap, nan_color=nan_color)
    mapit!(ax, clon, mypolys(clon), color=:transparent, strokecolor=:black, strokewidth=1)
    mylatlons!(ax, latticks45, lonticks60)
    # Labels
    text!(ax, "IND",  position=( 70, -20), textsize=12)
    text!(ax, "NPAC", position=(180,  25), textsize=12)
    text!(ax, "SPAC", position=(200, -30), textsize=12)
    text!(ax, "NATL", position=(305,  20), textsize=12)
    text!(ax, "SATL", position=(330, -30), textsize=12)
    text!(ax, "SO",   position=(200, -65), textsize=12)
    text!(ax, "ARC",  position=(180,  72), textsize=12)
    #(topscene, bbox = axs[i,1].scene.px_area, "($(Δz[1])–$(Δz[2]))m Nd", textsize=20, halign=:left, valign=:top, padding=(30,0,0,30), font="Dejavu Sans", color=:white)
    #    Label(topscene, bbox = axs[i,2].scene.px_area, "($(Δz[1])–$(Δz[2]))m εNd", textsize=20, halign=:left, valign=:top, padding=(30,0,0,30), font="Dejavu Sans", color=:white)
    nothing
end
fig = Figure(resolution = (500, 250), backgroundcolor=:white)
plot_diagnosis_mask_map!(fig)
trim!(fig.layout)
if use_GLMakie
    fig # show the output wiht GLMakie
else
    #save(joinpath(output_path, "εNd_sources_v2.$EXT"), scene)
    #save(joinpath(archive_path, "εNd_sources_$(lastcommit)_v2.$EXT"), scene, px_per_unit=4)
    save(joinpath(archive_path, "masks_maps_$(lastcommit).pdf"), fig)
    nothing # just so that no output is spat out
end


function plot_diagnosis_mask2_map!(fig)
    masks_map = permutedims(verticalmean(sum(i*m for (i,m) in enumerate(masks2)), grd))
    ax = fig[1,1] = Axis(fig, backgroundcolor=seafloor_color)
    mapit!(ax, clon, mypolys(clon), color=continent_color)
    hm = heatmap!(ax, sclons, lats, view(masks_map, ilon, :), colormap=masks2cmap, nan_color=nan_color)
    mapit!(ax, clon, mypolys(clon), color=:transparent, strokecolor=:black, strokewidth=1)
    mylatlons!(ax, latticks45, lonticks60)
    # Labels
    text!(ax, "NPAC", position=(180,  25), textsize=12)
    text!(ax, "NATL", position=(305,  20), textsize=12)
    #(topscene, bbox = axs[i,1].scene.px_area, "($(Δz[1])–$(Δz[2]))m Nd", textsize=20, halign=:left, valign=:top, padding=(30,0,0,30), font="Dejavu Sans", color=:white)
    #    Label(topscene, bbox = axs[i,2].scene.px_area, "($(Δz[1])–$(Δz[2]))m εNd", textsize=20, halign=:left, valign=:top, padding=(30,0,0,30), font="Dejavu Sans", color=:white)
    nothing
end
fig = Figure(resolution = (500, 250), backgroundcolor=:white)
plot_diagnosis_mask2_map!(fig)
trim!(fig.layout)
if use_GLMakie
    fig # show the output wiht GLMakie
else
    #save(joinpath(output_path, "εNd_sources_v2.$EXT"), scene)
    #save(joinpath(archive_path, "εNd_sources_$(lastcommit)_v2.$EXT"), scene, px_per_unit=4)
    save(joinpath(archive_path, "masks2_maps_$(lastcommit).pdf"), fig)
    nothing # just so that no output is spat out
end
=#





function plot_diagnosis_Omega_map!(fig)
    masks_sum = sum(i*m for (i,m) in enumerate(masks))
    masks_map = permutedims(verticalmean(masks_sum, grd))
    ax = fig[1,1] = Axis(fig, backgroundcolor=seafloor_color, aspect=DataAspect())
    mapit!(ax, clon, mypolys(clon), color=continent_color)
    hm = heatmap!(ax, sclons, lats, view(masks_map, ilon, :), colormap=Ωcmap[[3, 2, 1]], nan_color=nan_color)
    mapit!(ax, clon, mypolys(clon), color=:transparent, strokecolor=:black, strokewidth=1)
    #xlims!(ax, (-110, 35))
    mylatlons!(ax, latticks45, lonticks60)
    # Labels
    text!(ax, "N tag", position=(360 - 30, latN), align = (:center, :bottom))
    text!(ax, "S tag", position=(360 - 25, latS), align = (:center, :top))
    #(topscene, bbox = axs[i,1].scene.px_area, "($(Δz[1])–$(Δz[2]))m Nd", textsize=20, halign=:left, valign=:top, padding=(30,0,0,30), font="Dejavu Sans", color=:white)
    #    Label(topscene, bbox = axs[i,2].scene.px_area, "($(Δz[1])–$(Δz[2]))m εNd", textsize=20, halign=:left, valign=:top, padding=(30,0,0,30), font="Dejavu Sans", color=:white)
    nothing
end
fig = Figure(resolution = (600, 325), backgroundcolor=:white)
plot_diagnosis_Omega_map!(fig)
trim!(fig.layout)
if use_GLMakie
    fig # show the output wiht GLMakie
else
    #save(joinpath(output_path, "εNd_sources_v2.$EXT"), scene)
    #save(joinpath(archive_path, "εNd_sources_$(lastcommit)_v2.$EXT"), scene, px_per_unit=4)
    save(joinpath(archive_path, "Omega_maps_$(lastcommit).pdf"), fig)
    nothing # just so that no output is spat out
end
