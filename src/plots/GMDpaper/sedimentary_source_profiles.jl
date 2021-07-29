#================================================
Profiles
================================================#
use_GLMakie = false
include("../plots_setup_Nd.jl")


# Create the figure
fig = Figure()
use_GLMakie && display(fig)

isbasins = [isatlantic2, ispacific2, isantarctic]
basins = ["ATL", "PAC", "SO"]
colors = ColorSchemes.colorschemes[:tableau_colorblind][[5,6,4]]


function sed_source_profiles!(fig)
    zbot = sort(unique(bottomdepthvec(grd)))
    ztop = sort(unique(topdepthvec(grd)))
    yticks = vcat(0, zbot)
    y = reduce(vcat, [zt, zb] for (zt, zb) in zip(ztop, zbot))
    label_opts = (textsize=20, halign=:right, valign=:bottom, padding=(10,10,5,5), font=labelfont, color=:black)
    # left panel is just plain source
    ax = fig[1,1] = Axis(fig)
    u = u"pmol/cm^2/yr"
    upref = upreferred(u)

    z = 0:1:6000
    v = ustrip.(u, ϕ(p).(z) .* upref)

    # Rearrange as a step function to match model source
    lines!(ax, v, z)
    ylims!(ax, (6000, -50))
    #xlims!(ax, (0, 1.05maximum(v)))
    xlims!(ax, (0, maximum(ax.finallimits[])[1]))
    ax.yticks = 0:1000:6000
    ax.xlabel = "local (per unit area)\nsedimentary source ($(u))"
    ax.ylabel = "depth (m)"
    Label(fig, bbox = ax.scene.px_area, panellabels[1]; label_opts...)
    # panel for integrated source
    ax = fig[1,2] = Axis(fig)
    u∫dxdy = u"kmol/m/yr" 
    ∫dxdy_s_sed = ∫dxdy(s_sed(p) * upreferred(uDNd) / u"s", grd) .|> u∫dxdy
    x = vcat(0, repeat(ustrip.(∫dxdy_s_sed), inner=2), 0)
    #lines!(ax, x, vcat(0, y, maximum(zbot)))
    poly!(ax, Point2f0.(zip(x, vcat(0, y, maximum(zbot)))), color=ColorSchemes.colorschemes[:tableau_colorblind][1])
    ylims!(ax, (6000, -50))
    xlims!(ax, (0, maximum(ax.finallimits[])[1]))
    ax.yticks = 0:1000:6000
    ax.xlabel = "horizontally integrated\nsedimentary source ($(u∫dxdy))"
    hideydecorations!(ax, grid=false)
    Label(fig, bbox = ax.scene.px_area, panellabels[2], ; label_opts...)

    nothing
end

# Create the plot
sed_source_profiles!(fig)

# Label axes

# Add labels

#save(joinpath(output_path, "Nd_Makie_profiles.png"), scene)
#save(joinpath(archive_path, "Nd_profiles_$(lastcommit)_run$(run_num).png"), scene, px_per_unit=4)
save(joinpath(archive_path, "sedimentary_source_profiles_$(lastcommit)_run$(run_num).pdf"), fig)

nothing
