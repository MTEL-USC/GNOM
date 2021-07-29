
# Compute the diagnostics only once, unless rediagnose=true
(!isdefined(Main, :fDNd_wtags) || rediagnose) && include("GMD_diagnostics_setup.jl")

# Create the figure
fig = Figure(resolution=(600, 900))
use_GLMakie && display(fig)



function plot_conservative_ε!(fig)

    ΔεNd = ε_conservative - εNd
    δεlevels2 = δεlevels[δεlevels .≠ 0]

    axisopts = (xgridvisible=false, ygridvisible=false, backgroundcolor=land_color, ylabel="Depth (m)")
    opts = (mask=ATL, nan_color=nan_color, extendhigh=:auto, extendlow=:auto, linewidth=0.5)
    εopts = (colormap=εcmap, levels=εlevels, colorrange=εclims)
    δεopts = (colormap=δεcmap, levels=δεlevels2, colorrange=δεclims)
    
    labels = ("Modelled εNd", "Conservative εNd", "Difference")

    function commons!(ax)
        vlines!(ax, [latS, latN], linestyle=:dash, color=:black)
        myxlats!(ax, latticks30)
    end

    axs = Vector{Any}(undef, 3)

    # Model εNd
    ax = axs[1] = fig[1,1] = Axis(fig; axisopts...)
    _, εhm1 = generic_ZA!(ax, εNd .|> uεNd, grd; opts..., εopts...)
    myxlats!(ax, latticks30)
    commons!(ax)
    hidexdecorations!(ax, ticks=false, grid=false)

    # Conservative εNd
    ax = axs[2] = fig[2,1] = Axis(fig; axisopts...)
    _, εhm2 = generic_ZA!(ax, ε_conservative .|> uεNd, grd; opts..., εopts...)
    commons!(ax)
    hidexdecorations!(ax, ticks=false, grid=false)

    # Difference
    ax = axs[3] = fig[3,1] = Axis(fig; axisopts...)
    _, δεhm = generic_ZA!(ax, ΔεNd .|> uεNd, grd; opts..., δεopts...)
    commons!(ax)


    # colorbars
    cbar1 = fig[1:2, end+1] = Colorbar(fig, εhm1; label="εNd ($uεNd)", vertical=true, ticks=εlevels[1:5:end])
    cbar1.height = Relative(3/4)
    cbar1.width = 20
    cbar1.tellwidth = true
    cbar2 = fig[3, end] = Colorbar(fig, δεhm; label="Δ(εNd) ($uεNd)", vertical=true, ticks=δεlevels2)
    cbar2.height = Relative(1)
    cbar2.tickformat = x -> map(x -> x > 0 ? string("+", round(Int, x)) : string(round(Int, x)), x)
    cbar2.width = 20
    cbar2.tellwidth = true
    nothing

    # labels
    topscene = Scene(fig.scene)
    for i in 1:3
        Label(topscene, bbox = axs[i].scene.px_area, string(panellabels[i], "   ", labels[i]), textsize=20, halign=:left, valign=:bottom, padding=(10,0,5,0), font=labelfont, color=:white)
    end

    #zbot = sort(unique(bottomdepthvec(grd)))
    #ztop = sort(unique(topdepthvec(grd)))
    #yticks = vcat(0, zbot)
    #y = reduce(vcat, [zt, zb] for (zt, zb) in zip(ztop, zbot))
    #label_opts = (textsize=20, halign=:right, valign=:bottom, padding=(10,10,5,5), font=labelfont, color=:black)
    ## left panel is just plain source
    #ax = fig[1,1] = Axis(fig)
    #u = u"pmol/cm^2/yr"
    #upref = upreferred(u)

    #z = 0:1:6000
    #v = ustrip.(u, ϕ(p).(z) .* upref)

    ## Rearrange as a step function to match model source
    #lines!(ax, v, z)
    #ylims!(ax, (6000, -50))
    ##xlims!(ax, (0, 1.05maximum(v)))
    #xlims!(ax, (0, maximum(ax.finallimits[])[1]))
    #ax.yticks = 0:1000:6000
    #ax.xlabel = "local (per unit area)\nsedimentary source ($(u))"
    #ax.ylabel = "depth (m)"
    #Label(fig, bbox = ax.scene.px_area, panellabels[1]; label_opts...)
    ## panel for integrated source
    #ax = fig[1,2] = Axis(fig)
    #u∫dxdy = u"kmol/m/yr" 
    #∫dxdy_s_sed = ∫dxdy(s_sed(p) * upreferred(uDNd) / u"s", grd) .|> u∫dxdy
    #x = vcat(0, repeat(ustrip.(∫dxdy_s_sed), inner=2), 0)
    ##lines!(ax, x, vcat(0, y, maximum(zbot)))
    #poly!(ax, Point2f0.(zip(x, vcat(0, y, maximum(zbot)))), color=ColorSchemes.colorschemes[:tableau_colorblind][1])
    #ylims!(ax, (6000, -50))
    #xlims!(ax, (0, maximum(ax.finallimits[])[1]))
    #ax.yticks = 0:1000:6000
    #ax.xlabel = "horizontally integrated\nsedimentary source ($(u∫dxdy))"
    #hideydecorations!(ax, grid=false)
    #Label(fig, bbox = ax.scene.px_area, panellabels[2], ; label_opts...)

    nothing
end

# Create the plot
plot_conservative_ε!(fig)

# Label axes

# Add labels

#save(joinpath(output_path, "Nd_Makie_profiles.png"), scene)
#save(joinpath(archive_path, "Nd_profiles_$(lastcommit)_run$(run_num).png"), scene, px_per_unit=4)
save(joinpath(archive_path, "conservative_eNd_$(lastcommit)_run$(run_num).pdf"), fig)

nothing
