
# Compute the diagnostics only once, unless rediagnose=true
(!isdefined(Main, :fDNd_wtags) || rediagnose) && include("GMD_diagnostics_setup.jl")

function plot_wtag_fNd!(fig)

    u = per100
    ATL_ZA_fDNd_wtags = NamedTuple(j=>zonalaverage(ustrip.(v), grd, ATL) for (j,v) in pairs(fDNd_wtags))
	ATL_ZA_max = max.(collect(ATL_ZA_fDNd_wtags)...)

    #depth = [0; cumsum(ustrip.(grd.δdepth))]
    #lat = [ustrip(grd.lat[1]) - ustrip.(grd.δlat[1])/2; ustrip(grd.lat) + ustrip.(grd.δlat)/2]
    depth = ustrip.(grd.depth)
    lat = ustrip.(grd.lat)

    axisopts = (xgridvisible=false, ygridvisible=false, backgroundcolor=land_color, ylabel="Depth (m)")
    opts = (mask=ATL, nan_color=nan_color, linewidth=0.5, levels=10:10:100, filllevels=0:10:100, colorrange=(0,100))

    axs = Vector{Any}(undef, 4)
    hms = Vector{Any}(undef, 3)
    hms2 = Vector{Any}(undef, 3)

    ax = axs[1] = fig[1,1] = Axis(fig; axisopts...)
    _, hms2[1] = generic_ZA!(ax, ATL_ZA_fDNd_wtags.N, lat, depth; opts..., colormap=cgrad([:white, Ωcmap[2]]))
    myxlats!(ax, latticks30)
    vlines!(ax, [latS, latN], linestyle=:dash, color=:black)
    hidexdecorations!(ax)

    ax = axs[2] = fig[2,1] = Axis(fig; axisopts...)
    _, hms2[2] = generic_ZA!(ax, ATL_ZA_fDNd_wtags.S, lat, depth; opts..., colormap=cgrad([:white, Ωcmap[1]]))
    myxlats!(ax, latticks30)
    vlines!(ax, [latS, latN], linestyle=:dash, color=:black)
    hidexdecorations!(ax)

    ax = axs[3] = fig[3,1] = Axis(fig; axisopts...)
    _, hms2[3] = generic_ZA!(ax, ATL_ZA_fDNd_wtags.U, lat, depth; opts..., colormap=cgrad([:white, Ωcmap[3]]))
    myxlats!(ax, latticks30)
    vlines!(ax, [latS, latN], linestyle=:dash, color=:black)
    hidexdecorations!(ax)

    ax = axs[4] = fig[4,1] = Axis(fig; axisopts...)
    _, hm0 = generic_ZA!(ax, 0 * ATL_ZA_max, lat, depth; opts..., levels=0:10:100)
    # Model εNd
    for (i, (x, c)) in enumerate(zip(ATL_ZA_fDNd_wtags, Ωcmap[[2,1,3]]))
        x = copy(x)
		x[x .< ATL_ZA_max] .= NaN
        _, hms[i] = generic_ZA!(ax, x, lat, depth; opts..., colormap=cgrad([:white, c]))
    end
    myxlats!(ax, latticks30)
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
    cbar3 = fig[3, end] = Colorbar(fig, hms2[3]; label="% Nd from neither N or S", vertical=true, ticks=0:20:100)
    cbar3.height = Relative(1)
    cbar3.width = 20
    cbar3.tellwidth = true
    nothing

    # labels
    labels = ("Nd from N", "Nd from S", "Nd not from N or S", "Dominant fraction")

    for i in 1:4
        Label(fig, bbox = axs[i].scene.px_area, string(panellabels[i], "  ", labels[i]), textsize=20, halign=:left, valign=:bottom, padding=(10,0,5,0), font=labelfont, color=:white)
    end



    nothing
end

# Create the plot
fig = Figure(resolution=(600, 1200))
plot_wtag_fNd!(fig)

if use_GLMakie
    display(fig) # show the output wiht GLMakie
else
    save(joinpath(archive_path, "water-tagged_Nd_$(lastcommit)_run$(run_num).pdf"), fig)
    nothing # just so that no output is spat out
end

nothing
