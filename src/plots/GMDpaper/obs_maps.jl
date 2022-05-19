#include("load.jl") # this is the model-independent load.jl script
#include("../../modelling/obs.jl")
include("../plots_setup_Nd.jl")


mcol = [RGBA(c.r, c.g, c.b, 1.0) for c in ColorSchemes.okabe_ito[1:3]]
mshp = [:dtriangle, :diamond, :hexagon]

surfacemask = horizontalslice(ones(count(iswet(grd))), grd, depth=0)
function plot_obs_maps!(fig)
    axs = Array{Any,2}(undef, (2, 1))
    axs[1,1] = fig[1,1] = Axis(fig, backgroundcolor=seafloor_color)
    ax = axs[1,1]
    scts = Vector{Any}(undef, 3) # to store scatter marker for legend
    mapit!(ax, clon, mypolys(clon), color=continent_color)
    sources = ["van de Flierdt", "GEOTRACES IDP17", "post IDP17"]
    for (i,source) in enumerate(sources)
        df = select(DNdobs, [:lat, :lon])[DNdobs.source .== source,:]
        scts[i] = myscatter!(ax, centerlon.(df.lon), df.lat; marker=mshp[i], color=mcol[i], strokewidth=1, markersize)
    end
    mapit!(ax, clon, mypolys(clon), color=:transparent, strokewidth=1, strokecolor=:black)
    mylatlons!(ax, latticks45, lonticks60)
    hidexdecorations!(ax, ticks=false, grid=false)
    # plot εNd of source
    axs[2,1] = fig[2,1] = Axis(fig, backgroundcolor=seafloor_color)
    ax = axs[2,1]
    mapit!(ax, clon, mypolys(clon), color=continent_color)
    for (i,source) in enumerate(sources)
        df = select(εNdobs, [:lat, :lon])[εNdobs.source .== source,:]
        myscatter!(ax, centerlon.(df.lon), df.lat; marker=mshp[i], color=mcol[i], strokewidth=1, markersize)
    end
    mapit!(ax, clon, mypolys(clon), color=:transparent, strokewidth=1, strokecolor=:black)
    mylatlons!(ax, latticks45, lonticks60)
    # annotations (must come after?)

    Label(fig, bbox = axs[1,1].scene.px_area, string("(a)"), textsize=20, halign=:left, valign=:bottom, padding=(10,0,5,0), font="Dejavu Sans", color=:white)
    Label(fig, bbox = axs[2,1].scene.px_area, string("(b)"), textsize=20, halign=:left, valign=:bottom, padding=(10,0,5,0), font="Dejavu Sans", color=:white)
    Label(fig, bbox = axs[1,1].scene.px_area, "Nd obs", textsize=20, halign=:left, valign=:top, padding=(50,0,0,50), font="Dejavu Sans", color=:white)
    Label(fig, bbox = axs[2,1].scene.px_area, "εNd obs", textsize=20, halign=:left, valign=:top, padding=(50,0,0,50), font="Dejavu Sans", color=:white)
    nothing

    # Legend
    markers = [MarkerElement(; marker, color, markersize=15, strokewidth=1) for (marker,color) in zip(mshp, mcol)]
    leg = Legend(fig, markers, sources)
    leg.orientation = :horizontal
    leg.tellheight = true
    fig[3,1] = leg
end
fig = Figure(resolution = (650, 700), backgroundcolor=:white)
plot_obs_maps!(fig)
trim!(fig.layout)
if use_GLMakie
    display(fig) # show the output wiht GLMakie
else
    save(joinpath(archive_path, "obs_maps.pdf"), fig)
    nothing # just so that no output is spat out
end
