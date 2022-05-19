include("../plots_setup_Nd.jl")


u∫dxdy_s = u"kmol/m/yr"
# plot the εNd of sources

surfacemask = horizontalslice(ones(count(iswet(grd))), grd, depth=0)
function plot_εNd_sinks!(fig, fun)
    islog = fun == log10
    u = islog ? u"mol/m^2/yr" : u"μmol/m^2/yr"
    hms = Vector{Any}(undef, 1)
    axs = Array{Any,2}(undef, (length(instances(ScavenginParticle)),2))
    # For profiles
    zbot = sort(unique(bottomdepthvec(grd)))
    ztop = sort(unique(topdepthvec(grd)))
    yticks = vcat(0, zbot)
    ypro = reduce(vcat, [zt, zb] for (zt, zb) in zip(ztop, zbot))
    for (i, t) in enumerate(instances(ScavenginParticle))
        Tx = T_D(t, p) * DNd * u"mol/m^3/s" # each sink is the vertical integral of Tx
        # Might freeze if DNd is a view of a SciML solution because it dispatches to dense matmul
        ∫dzsink = ∫dz(Tx, grd)
        ∫dxdysink = ∫dxdy(Tx, grd)
        # replace negative values from numerical noise when doing T * x
        ∫dzsink[findall(∫dzsink .≤ 0u"mol/m^2/s")] .= 1e-50 * u"mol/m^2/s"
        fun∫dzsink = fun.(ustrip.(u, permutedims(∫dzsink, (2,1))))
        ax = axs[i,1] = fig[i,1] = Axis(fig, backgroundcolor=seafloor_color)
        mapit!(ax, clon, mypolys(clon), color=continent_color)
        hms[1] = heatmap!(ax, sclons, lats, view(fun∫dzsink, ilon, :),
                          colormap = islog ? logσcmap : σcmap,
                          colorrange = islog ? logσclims : σclims,
                          nan_color=nan_color)
        mapit!(ax, clon, mypolys(clon), color=:transparent, strokecolor=:black, strokewidth=1)
        mylatlons!(ax, latticks45, lonticks60)
        i≠length(instances(ScavenginParticle)) && hidexdecorations!(ax, ticks=false, grid=false)
        # Profiles on the side
        ax = axs[i,2] = fig[i,2] = Axis(fig, width=300)
        x = vcat(0, repeat(ustrip.(u∫dxdy_s, ∫dxdysink), inner=2), 0)
        poly!(ax, Point2f0.(zip(x, vcat(0, ypro, maximum(zbot)))), color=ColorSchemes.colorschemes[:tableau_colorblind][1])
        ylims!(ax, (6000, -50))
        ax.yticks = 0:1000:6000
        ax.ylabel = "depth (m)"
        ax.xlabel = "$(u∫dxdy_s)"
        # Do not hide x ticks because different scales? Use log? (Should not use logscaled bars though)
        #i≠length(instances(ScavenginParticle)) && hidexdecorations!(ax, ticks=false, grid=false)
        # Cannot hide x ticks if axes are not linked...
        # # TODO maybe replace scavenging in ∫dxdy by the particle concentration instead? (good sanity check too)
    end
    # annotations (must come after?)

    for (i, t) in enumerate(instances(ScavenginParticle))
        Label(fig, bbox = axs[i,1].scene.px_area, panellabels[i], textsize=20, halign=:left, valign=:bottom, padding=(10,0,5,0), font=labelfont, color=:white)
        Label(fig, bbox = axs[i,2].scene.px_area, panellabels[i+length(instances(ScavenginParticle))], textsize=20, halign=:left, valign=:bottom, padding=(10,0,5,0), font=labelfont, color=:black)
        Label(fig, bbox = axs[i,1].scene.px_area, string(t)[2:end], textsize=20, halign=:left, valign=:top, padding=(70,0,0,50), font=labelfont, color=:white)
    end
    # colorbars
    label = islog ? "log₁₀(scavₖ / ($u))" : "sₖ ($u)"
    cbar = fig[end+1, 1] = Colorbar(fig, hms[1], label=label, vertical=false, flipaxis=false, ticklabelalign = (:center, :top))
    cbar.width = Relative(3/4)
    cbar.height = 30
    cbar.tellheight = true
    nothing
end

for fun in (log10, identity)
    local fig = Figure(resolution = (1000, 1500), backgroundcolor=:white)
    plot_εNd_sinks!(fig, fun)
    trim!(fig.layout)
    if use_GLMakie
        display(fig) # show the output wiht GLMakie
    else
        str = (fun == log10) ? "log_" : ""
        save(joinpath(archive_path, "sinks_$(str)maps_$(lastcommit)_run$(run_num).pdf"), fig)
        nothing # just so that no output is spat out
    end
end
