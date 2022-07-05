#================================================
Profiles
================================================#
include("../plots_setup_Nd.jl")



isbasins = [isatlantic2, ispacific2, isantarctic]
basins = ["ATL", "PAC", "SO"]
colors = ColorSchemes.colorschemes[:tableau_colorblind][[5,6,4]]



function sed_source_profiles!(fig)
    zbot = sort(unique(bottomdepthvec(grd)))
    ztop = sort(unique(topdepthvec(grd)))
    yticks = vcat(0, zbot)
    y = reduce(vcat, [zt, zb] for (zt, zb) in zip(ztop, zbot))
    label_opts = (textsize=20, halign=:left, valign=:bottom, padding=(10,10,5,5), font=labelfont, color=:black)

    # left panel is just ϕ(z)
    ax = fig[1,2] = Axis(fig)
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
    ax.xlabel = "Base sed. flux, ϕ(z) ($(u))"
    ax.ylabel = "depth (m)"
    Label(fig, bbox = ax.scene.px_area, panellabels[3]; label_opts...)

    # panel for integrated source
    ax = fig[2,2] = Axis(fig)
    u∫dxdy = u"kmol/m/yr"
    ∫dxdy_s_sed = ∫dxdy(s_sed(p) * upreferred(uDNd) / s, grd) .|> u∫dxdy
    x = vcat(0, repeat(ustrip.(∫dxdy_s_sed), inner=2), 0)
    #lines!(ax, x, vcat(0, y, maximum(zbot)))
    poly!(ax, Point2f0.(zip(x, vcat(0, y, maximum(zbot)))), color=ColorSchemes.colorschemes[:tableau_colorblind][1])
    ylims!(ax, (6000, -50))
    xlims!(ax, (0, maximum(ax.finallimits[])[1]))
    ax.yticks = 0:1000:6000
    ax.ylabel = "Depth (m)"
    ax.xlabel = "∫dxdy sed. source ($(u∫dxdy))"
    #hideydecorations!(ax, grid=false)
    Label(fig, bbox = ax.scene.px_area, panellabels[4], ; label_opts...)

    # panel for alpha curve
    ax1 = fig[1,1] = Axis(fig)
    εs = upreferred.(collect(range(εclims..., length=1001) * per10000))
    αs = α_quad(εs, p)
    vlines!(ax1, [ustrip.(p.α_c)], linestyle=:dash, color=:gray)
    lines!(ax1, ustrip.(per10000, εs), αs)
    ax1.xlabel = "εNd (‱)"
    ax1.ylabel = "Scaling factor α(εNd)"
    ylims!(ax1, low=0.0)
    xlims!(ax1, εclims)
    Label(fig, bbox = ax1.scene.px_area, panellabels[1], ; label_opts...)


    # panel for shifted epsilon
    ax = fig[2,1] = Axis(fig)
    @unpack α_a, α_c, σ_ε = p
    ε_eff = shifted_ε.(εs, σ_ε, α_a, α_c, ε10)
    vlines!(ax, [ustrip.(p.α_c)], linestyle=:dash, color=:gray)
    hlines!(ax, [0.0], color=:gray)
    lines!(ax, ustrip.(per10000, εs), ustrip.(per10000, ε_eff - εs))
    ax.xlabel = "In situ εNd (‱)"
    ax.ylabel = "Released εNd − εNd (‱)"
    xlims!(ax, εclims)
    Label(fig, bbox = ax.scene.px_area, panellabels[2], ; label_opts...)


    # α map
    α2D = permutedims(rearrange_into_3Darray(α_quad(p), grd)[:,:,1])
    innan = findall(.!isnan.(α2D))
    colorrange = extrema(α2D[innan]) .* (0,1)
    ax = fig[3:4,:] = Axis(fig, backgroundcolor=:gray20)
    mapit!(ax, clon, mypolys(clon), color=:gray50)
    hm = heatmap!(ax, sclons, lats, view(α2D, ilon, :), colormap=αcmap; nan_color, colorrange)#, colorrange=αlims)
    mapit!(ax, clon, mypolys(clon), color=:transparent, strokecolor=:black, strokewidth=1)
    # Better lat/lon ticks
    mylatlons!(ax, latticks30, lonticks60)
    # colorbar
    cbar = fig[end+1,:] = Colorbar(fig, hm, label="Scaling factor α(εNd)", vertical=false, flipaxis=false, ticklabelalign = (:center, :top))#, ticks=αlevs)
    cbar.width = Relative(3/4)
    cbar.height = 30
    cbar.tellheight = true

    # label

    Label(fig, bbox = ax.scene.px_area, panellabels[5]; label_opts...)
    #textsize=20, halign=:left, valign=:bottom, padding=(10,0,5,0), font=labelfont, color=:black)

    nothing
end

# Create the figure
fig = Figure(resolution=(800,1000))
sed_source_profiles!(fig)

if use_GLMakie
    display(fig) # show the output wiht GLMakie
else
    save(joinpath(archive_path, "sedimentary_source_details_$(lastcommit)_run$(run_num).pdf"), fig)
    nothing # just so that no output is spat out
end

nothing
