# Compute the diagnostics only once, unless rediagnose=true
(!isdefined(Main, :fDNd_wtags) || rediagnose) && include("GMD_diagnostics_setup.jl")


function plot_diagnosis_Omega_map!(fig)
    masks_sum = sum(i*m for (i,m) in enumerate(masks))
    masks_map = permutedims(verticalmean(masks_sum, grd))
    ax = fig[1,1] = Axis(fig, backgroundcolor=seafloor_color, aspect=DataAspect())
    mapit!(ax, clon, mypolys(clon), color=continent_color)
    hm = heatmap!(ax, sclons, lats, view(masks_map, ilon, :), colormap=[:lightgray; Ωcmap[[2, 1]]], nan_color=nan_color)
    mapit!(ax, clon, mypolys(clon), color=:transparent, strokecolor=:black, strokewidth=1)
    #xlims!(ax, (-110, 35))
    mylatlons!(ax, latticks45, lonticks60)
    # Labels
    text!(ax, L"\Omega_{\mathrm{N}}", position=(360 - 30, latN+4), align = (:center, :bottom), textsize=30)
    text!(ax, L"\Omega_{\mathrm{S}}", position=(360 - 25, latS-4), align = (:center, :top), textsize=30)
    xlims!(ax, (250, 360+wlon))
    #(fig, bbox = axs[i,1].scene.px_area, "($(Δz[1])–$(Δz[2]))m Nd", textsize=20, halign=:left, valign=:top, padding=(30,0,0,30), font="Dejavu Sans", color=:white)
    #    Label(fig, bbox = axs[i,2].scene.px_area, "($(Δz[1])–$(Δz[2]))m εNd", textsize=20, halign=:left, valign=:top, padding=(30,0,0,30), font="Dejavu Sans", color=:white)
    nothing
end
fig = Figure(resolution = (400, 440), backgroundcolor=:white)
plot_diagnosis_Omega_map!(fig)
trim!(fig.layout)
if use_GLMakie
    display(fig) # show the output wiht GLMakie
else
    save(joinpath(archive_path, "Omega_maps_$(lastcommit).pdf"), fig)
    nothing # just so that no output is spat out
end
