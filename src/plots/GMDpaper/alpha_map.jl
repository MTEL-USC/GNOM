use_GLMakie = false
include("../plots_setup_Nd.jl")


# Create the figure
fig = Figure()
use_GLMakie && display(fig)

# levels and limits for colorbar
#αlevs = 0.99:0.002:1.01
#αlims = extrema(αlevs)

# α colormap
αcmap = cgrad(:solar)

surfacemask = horizontalslice(ones(count(iswet(grd))), grd, depth=0)
function plot_alpha_map!(fig)
    α2D = permutedims(rearrange_into_3Darray(α_quad(p), grd)[:,:,1])
    innan = findall(.!isnan.(α2D))
    colorrange = extrema(log10.(α2D[innan]))
    ax = fig[1,1] = Axis(fig, backgroundcolor=:gray20)
    mapit!(ax, clon, mypolys(clon), color=:gray50)
    hm = heatmap!(ax, sclons, lats, view(log10.(α2D), ilon, :), colormap=αcmap; nan_color, colorrange)#, colorrange=αlims)
    mapit!(ax, clon, mypolys(clon), color=:transparent, strokecolor=:black, strokewidth=1)
    # Better lat/lon ticks
    mylatlons!(ax, latticks30, lonticks60)
    # colorbar
    cbar = fig[end+1, 1] = Colorbar(fig, hm, label="Reactivity scaling factor log₁₀(α)", vertical=false, flipaxis=false, ticklabelalign = (:center, :top))#, ticks=αlevs)
    cbar.width = Relative(3/4)
    cbar.height = 30
    cbar.tellheight = true
    nothing
end
fig = Figure(resolution = (700, 500), backgroundcolor=:white)
plot_alpha_map!(fig)
trim!(fig.layout)
if use_GLMakie
    fig # show the output wiht GLMakie
else
    save(joinpath(archive_path, "logalpha_map_$(lastcommit)_run$(run_num).pdf"), fig)
    nothing # just so that no output is spat out
end



# Create the figure
fig = Figure()
use_GLMakie && display(fig)
# α colormap
αcmap = cgrad(:lajolla, rev=true)

surfacemask = horizontalslice(ones(count(iswet(grd))), grd, depth=0)
function plot_alpha_map!(fig)
    α2D = permutedims(rearrange_into_3Darray(α_quad(p), grd)[:,:,1])
    innan = findall(.!isnan.(α2D))
    colorrange = extrema(α2D[innan]) .* (0,1)
    ax = fig[1,1] = Axis(fig, backgroundcolor=:gray20)
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
    nothing
end
fig = Figure(resolution = (700, 500), backgroundcolor=:white)
plot_alpha_map!(fig)
trim!(fig.layout)
if use_GLMakie
    fig # show the output wiht GLMakie
else
    save(joinpath(archive_path, "alpha_map_$(lastcommit)_run$(run_num).pdf"), fig)
    nothing # just so that no output is spat out
end
