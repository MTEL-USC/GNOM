
# Compute the diagnostics only once, unless rediagnose=true
(!isdefined(Main, :fDNd_wtags) || rediagnose) && include("GMD_diagnostics_setup.jl")

# Create the figure
fig = Figure(resolution=(600, 900))
use_GLMakie && display(fig)



function plot_wtag_fNd!(fig)

    u = u"percent"
    ATL_ZA_fDNd_wtags = NamedTuple(j=>zonalaverage(ustrip.(v), grd, ATL) for (j,v) in pairs(fDNd_wtags))
	ATL_ZA_max = max.(collect(ATL_ZA_fDNd_wtags)...)
    
    #depth = [0; cumsum(ustrip.(grd.δdepth))]
    #lat = [ustrip(grd.lat[1]) - ustrip.(grd.δlat[1])/2; ustrip(grd.lat) + ustrip.(grd.δlat)/2]
    depth = ustrip.(grd.depth)
    lat = ustrip.(grd.lat)





    axisopts = (xgridvisible=false, ygridvisible=false, backgroundcolor=land_color, ylabel="Depth (m)")
    opts = (mask=ATL, nan_color=nan_color, linewidth=0.5, levels=10:10:100, filllevels=0:10:100, colorrange=(0,100))

    function commons!(ax)
    end

    axs = Vector{Any}(undef, 3)
    hms = Vector{Any}(undef, 3)
    hms2 = Vector{Any}(undef, 2)

    ax = axs[1] = fig[1,1] = Axis(fig; axisopts...)
    _, hms2[1] = generic_ZA!(ax, ATL_ZA_fDNd_wtags.N, lat, depth; opts..., colormap=cgrad([:white, Ωcmap[2]]))
    myxlats!(ax, latticks30)
    commons!(ax)
    vlines!(ax, [latS, latN], linestyle=:dash, color=:black)
    hidexdecorations!(ax)

    ax = axs[2] = fig[2,1] = Axis(fig; axisopts...)
    _, hms2[2] = generic_ZA!(ax, ATL_ZA_fDNd_wtags.S, lat, depth; opts..., colormap=cgrad([:white, Ωcmap[1]]))
    myxlats!(ax, latticks30)
    commons!(ax)
    vlines!(ax, [latS, latN], linestyle=:dash, color=:black)
    hidexdecorations!(ax)

    ax = axs[3] = fig[3,1] = Axis(fig; axisopts...)
    _, hm0 = generic_ZA!(ax, 0 * ATL_ZA_max, lat, depth; opts..., levels=0:10:100)
    # Model εNd
    for (i, (x, c)) in enumerate(zip(ATL_ZA_fDNd_wtags, Ωcmap[[2,1,3]]))
        x = copy(x)
		x[x .< ATL_ZA_max] .= NaN
        _, hms[i] = generic_ZA!(ax, x, lat, depth; opts..., colormap=cgrad([:white, c]))
    end
    myxlats!(ax, latticks30)
    commons!(ax)
    vlines!(ax, [latS, latN], linestyle=:dash, color=:black)



    # colorbars
    cbar1 = fig[1, end+1] = Colorbar(fig, hms2[1]; label="% Nd from N region", vertical=true, ticks=0:20:100)
    cbar1.height = Relative(1)
    cbar1.width = 20
    cbar1.tellwidth = true
    cbar2 = fig[2, end] = Colorbar(fig, hms2[2]; label="% Nd from S region", vertical=true, ticks=0:20:100)
    cbar2.height = Relative(1)
    cbar2.width = 20
    cbar2.tellwidth = true
    cbar2 = fig[3, end] = Colorbar(fig, hms[3]; label="% Nd from neither N or S", vertical=true, ticks=0:20:100)
    cbar2.height = Relative(1)
    cbar2.width = 20
    cbar2.tellwidth = true
    nothing

    # labels
    labels = ("% Nd from North", "% Nd from South", "Dominant fraction", )
    topscene = Scene(fig.scene)
    for i in 1:3
        Label(topscene, bbox = axs[i].scene.px_area, string(panellabels[i], "  ", labels[i]), textsize=20, halign=:left, valign=:bottom, padding=(10,0,5,0), font=labelfont, color=:white)
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
plot_wtag_fNd!(fig)

# Label axes

# Add labels

#save(joinpath(output_path, "Nd_Makie_profiles.png"), scene)
#save(joinpath(archive_path, "Nd_profiles_$(lastcommit)_run$(run_num).png"), scene, px_per_unit=4)
save(joinpath(archive_path, "water-tagged_Nd_$(lastcommit)_run$(run_num).pdf"), fig)

nothing
