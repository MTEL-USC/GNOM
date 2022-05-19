include("../plots_setup_Nd.jl")





surfacemask = horizontalslice(ones(count(iswet(grd))), grd, depth=0)
function plot_horizontal_maps!(fig)
    depths_intervals = [(0,100), (100, 300), (300, 800), (800, 1500), (1500, 3000), (3000, 6000)]
    hms = Vector{Any}(undef, 2)
    axs = Array{Any,2}(undef, (length(depths_intervals), 2))
    # all maps (Nd source and εNd of source)
    for (i, Δz) in enumerate(depths_intervals)
        # s1 and s2 are the Nd source and the isotope source
        mask = Δz[1] .≤ depthvec(grd) .< Δz[2]
        ∫Nd = permutedims(verticalmean(DNd * upreferred(uDNd) .|> uDNd .|> ustrip, grd, mask))
        ∫εNd = permutedims(verticalmean(εNd * upreferred(uεNd) .|> uεNd .|> ustrip, grd, mask))
        # plot Nd source
        axs[i,1] = fig[i,1] = Axis(fig, backgroundcolor=seafloor_color)
        ax = axs[i,1]
        mapit!(ax, clon, mypolys(clon), color=continent_color)
        heatmap!(ax, sclons, lats, view(∫Nd, ilon, :), colormap=Ndcmap, nan_color=nan_color, colorrange=Ndclims)
        iobs = findall(Δz[1] .≤ ustrip.(DNdobs.depth) .< Δz[2])
        hms[1] = myscatter!(ax, centerlon.(DNdobs.lon[iobs]), DNdobs.lat[iobs]; color=DNdobs.value[iobs] * upreferred(uDNd) .|> uDNd .|> ustrip, colormap=Ndcmap, colorrange=Ndclims, markersize=markersize_maps)
        mapit!(ax, clon, mypolys(clon), color=:transparent, strokecolor=:black, strokewidth=1)
        mylatlons!(ax, latticks45, lonticks60)
        i≠length(depths_intervals) && hidexdecorations!(ax, ticks=false, grid=false)
        # plot εNd of source
        axs[i,2] = fig[i,2] = Axis(fig, backgroundcolor=seafloor_color)
        ax = axs[i,2]
        mapit!(ax, clon, mypolys(clon), color=continent_color)
        heatmap!(ax, sclons, lats, view(∫εNd, ilon, :), colormap=εcmap, nan_color=nan_color, colorrange=εclims)
        iobs = findall(Δz[1] .≤ ustrip.(εNdobs.depth) .< Δz[2])
        hms[2] = myscatter!(ax, centerlon.(εNdobs.lon[iobs]), εNdobs.lat[iobs]; color=εNdobs.value[iobs] * upreferred(uεNd) .|> uεNd .|> ustrip, colormap=εcmap, colorrange=εclims, markersize=markersize_maps)
        mapit!(ax, clon, mypolys(clon), color=:transparent, strokecolor=:black, strokewidth=1)
        mylatlons!(ax, latticks45, lonticks60)
        i≠length(depths_intervals) && hidexdecorations!(ax, ticks=false, grid=false)
        hideydecorations!(ax, ticks=false, grid=false)
    end
    # annotations (must come after?)
    for (i, Δz) in enumerate(depths_intervals)
        Label(fig, bbox = axs[i,1].scene.px_area, panellabels[i], textsize=20, halign=:left, valign=:bottom, padding=(10,0,5,0), font=labelfont, color=:white)
        Label(fig, bbox = axs[i,2].scene.px_area, panellabels[i+5], textsize=20, halign=:left, valign=:bottom, padding=(10,0,5,0), font=labelfont, color=:white)
        Label(fig, bbox = axs[i,1].scene.px_area, "Nd", textsize=20, halign=:left, valign=:top, padding=(60,0,0,35), font=labelfont, color=:white)
        Label(fig, bbox = axs[i,2].scene.px_area, "εNd", textsize=20, halign=:left, valign=:top, padding=(60,0,0,35), font=labelfont, color=:white)
        Label(fig, bbox = axs[i,1].scene.px_area, "($(Δz[1])–$(Δz[2]))m", textsize=10, halign=:left, valign=:top, padding=(60,0,0,60), font=labelfont, color=:white)
        Label(fig, bbox = axs[i,2].scene.px_area, "($(Δz[1])–$(Δz[2]))m", textsize=10, halign=:left, valign=:top, padding=(60,0,0,60), font=labelfont, color=:white)
    end
    # colorbars
    cbar1 = fig[end+1, 1] = Colorbar(fig, hms[1], label="Nd ($uDNd)", vertical=false, flipaxis=false, ticklabelalign = (:center, :top), ticks=Ndclims[1]:10:Ndclims[2])
    cbar1.width = Relative(3/4)
    cbar1.height = 30
    cbar1.tellheight = true
    cbar2 = fig[end, 2] = Colorbar(fig, hms[2], label="εNd ($uεNd)", vertical=false, flipaxis=false, ticklabelalign=(:center, :top), ticks=εclims[1]:5:εclims[2])
    cbar2.width = Relative(3/4)
    cbar2.height = 30
    cbar2.tellheight = true
    nothing
end
fig = Figure(resolution = (1200, 1800), backgroundcolor=:white)
plot_horizontal_maps!(fig)
trim!(fig.layout)
if use_GLMakie
    display(fig) # show the output wiht GLMakie
else
    #save(joinpath(output_path, "εNd_sources_v2.$EXT"), scene)
    #save(joinpath(archive_path, "εNd_sources_$(lastcommit)_run$(run_num)_v2.$EXT"), scene, px_per_unit=4)
    save(joinpath(archive_path, "horizontal_maps_$(lastcommit)_run$(run_num).pdf"), fig)
    nothing # just so that no output is spat out
end
