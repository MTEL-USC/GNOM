include("../plots_setup_Nd.jl")

# levels and limits for colorbar
#αlevs = 0.99:0.002:1.01
#αlims = extrema(αlevs)

# α curve is defined by α_quad
α_curve(ε) = α_quad(ε * per10000 |> upreferred, p)
εs = range(εclims..., length=200)
αs = α_curve.(εs)

# α curve vs ε values
function αcurve(fig)
    ax = fig[1,1] = Axis(fig)
    lines!(ax, εs, αs)
    ylims!(ax, low=0)
    ax.xlabel = "εNd (‱)"
    ax.ylabel = "Reactivity scaling factor α"
    ax.xticks = range(εclims..., step=5)
    ax.yticks = 0:1:maximum(αs)

    Label(fig, bbox = ax.scene.px_area, panellabels[1], textsize=20, halign=:left, valign=:bottom, padding=(10,0,5,0), font=labelfont, color=:black)
end



#=====================================#
#          α map in logscale          #
#=====================================#

# α colormap
αcmap = cgrad(:solar)

surfacemask = horizontalslice(ones(count(iswet(grd))), grd, depth=0)
function plot_alpha_map!(fig)
    αcurve(fig)
    # α map
    α2D = permutedims(rearrange_into_3Darray(α_quad(p), grd)[:,:,1])
    innan = findall(.!isnan.(α2D))
    colorrange = extrema(log10.(α2D[innan]))
    ax = fig[2,1] = Axis(fig, backgroundcolor=:gray20)
    mapit!(ax, clon, mypolys(clon), color=:gray50)
    hm = heatmap!(ax, sclons, lats, view(log10.(α2D), ilon, :), colormap=αcmap; nan_color, colorrange)#, colorrange=αlims)
    mapit!(ax, clon, mypolys(clon), color=:transparent, strokecolor=:black, strokewidth=1)
    # Better lat/lon ticks
    mylatlons!(ax, latticks30, lonticks60)
    # colorbar
    cbar = fig[end+1, 1] = Colorbar(fig, hm, label="Log reactivity scaling factor log₁₀(α)", vertical=false, flipaxis=false, ticklabelalign = (:center, :top))#, ticks=αlevs)
    cbar.width = Relative(3/4)
    cbar.height = 30
    cbar.tellheight = true

    # label

    Label(fig, bbox = ax.scene.px_area, panellabels[2], textsize=20, halign=:left, valign=:bottom, padding=(10,0,5,0), font=labelfont, color=:black)
    nothing
end
fig = Figure(resolution = (700, 1000), backgroundcolor=:white)
plot_alpha_map!(fig)
trim!(fig.layout)
if use_GLMakie
    fig # show the output wiht GLMakie
else
    save(joinpath(archive_path, "logalpha_map_$(lastcommit)_run$(run_num).pdf"), fig)
    nothing # just so that no output is spat out
end



#=====================================#
#          α map in linear scale      #
#=====================================#

# α colormap
αcmap = cgrad(:lajolla, rev=true)

surfacemask = horizontalslice(ones(count(iswet(grd))), grd, depth=0)
function plot_alpha_map!(fig)
    αcurve(fig)

    α2D = permutedims(rearrange_into_3Darray(α_quad(p), grd)[:,:,1])
    innan = findall(.!isnan.(α2D))
    colorrange = extrema(α2D[innan]) .* (0,1)
    ax = fig[2,1] = Axis(fig, backgroundcolor=:gray20)
    mapit!(ax, clon, mypolys(clon), color=:gray50)
    hm = heatmap!(ax, sclons, lats, view(α2D, ilon, :), colormap=αcmap; nan_color, colorrange)#, colorrange=αlims)
    mapit!(ax, clon, mypolys(clon), color=:transparent, strokecolor=:black, strokewidth=1)
    # Better lat/lon ticks
    mylatlons!(ax, latticks30, lonticks60)
    # colorbar
    cbar = fig[end+1, 1] = Colorbar(fig, hm, label="Reactivity scaling factor α", vertical=false, flipaxis=false, ticklabelalign = (:center, :top))#, ticks=αlevs)
    cbar.width = Relative(3/4)
    cbar.height = 30
    cbar.tellheight = true

    # label

    Label(fig, bbox = ax.scene.px_area, panellabels[2], textsize=20, halign=:left, valign=:bottom, padding=(10,0,5,0), font=labelfont, color=:black)
    nothing
end
fig = Figure(resolution = (700, 1000), backgroundcolor=:white)
plot_alpha_map!(fig)
trim!(fig.layout)
if use_GLMakie
    display(fig) # show the output wiht GLMakie
else
    save(joinpath(archive_path, "alpha_map_$(lastcommit)_run$(run_num).pdf"), fig)
    nothing # just so that no output is spat out
end
