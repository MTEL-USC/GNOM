
# Compute the diagnostics only once, unless rediagnose=true
(!isdefined(Main, :fDNd_wtags) || rediagnose) && include("GMD_diagnostics_setup.jl")



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

    for i in 1:3
        Label(fig, bbox = axs[i].scene.px_area, string(panellabels[i], "   ", labels[i]), textsize=20, halign=:left, valign=:bottom, padding=(10,0,5,0), font=labelfont, color=:white)
    end

    nothing
end


# Create the figure
fig = Figure(resolution=(600, 900))

# Create the plot
plot_conservative_ε!(fig)

if use_GLMakie
    display(fig) # show the output wiht GLMakie
else
    save(joinpath(archive_path, "conservative_eNd_$(lastcommit)_run$(run_num).pdf"), fig)
    nothing # just so that no output is spat out
end