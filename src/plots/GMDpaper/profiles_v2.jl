#================================================
Profiles
================================================#
include("../plots_setup_Nd.jl")

# Create the figure
fig = Figure(resolution=(1200, 900))


isbasins = [isatlantic2, ispacific2, isindian2, isantarctic]
basins = ["ATL", "PAC", "IND", "SO"]
basin_longnames = ["Atlantic", "Pacific", "Indian", "Southern Ocean"]
colors = ColorSchemes.colorschemes[:tableau_colorblind][[1,6,2,4]]
model_thickness = 2.5
obs_thickness = 1


# Function to average data
# that groups model and obs. by OCIM depths
function basin_profile_mean_and_std(isbasin, xmodel::Vector{T}, xmodelatobs, xdepthatobs, xobs, xwetobs, iobswet) where T
    ibasin = findall(isbasin(xobs.lat[iobswet], xobs.lon[iobswet], OCEANS))
    xmodelmean = T[]
    xmodelerror = T[]
    xobsmean = T[]
    xobserror = T[]
    depths = eltype(grd.depth)[]
    for (iz,z) in enumerate(grd.depth)
        idepth = findall(xdepthatobs[ibasin] .== ustrip(z))
        if !isempty(idepth)
            xmodel_atz = xmodelatobs[ibasin][idepth]
            xobs_atz = xwetobs[ibasin][idepth]
            push!(xmodelmean, mean(xmodel_atz))
            push!(xmodelerror, std(xmodel_atz, corrected=false))
            push!(xobsmean, mean(xobs_atz))
            push!(xobserror, std(xobs_atz, corrected=false))
            push!(depths, z)
        end
    end
    return ustrip.(xmodelmean), ustrip.(xmodelerror), ustrip.(xobsmean), ustrip.(xobserror), ustrip.(depths)
end

# Function for whole plot
function plot_profiles_v2!(fig)
    axs = Array{Any}(undef, (2, 4))
    for (itracer, (xmodel, xobs, ux)) in enumerate(zip((DNdmodel, εNdmodel), (DNdobs, εNdobs), (uDNd, uεNd)))
        xmodelatobs, xdepthatobs, iobswet, xwetobs = _locations(xmodel, xobs, ux)
        for (ibasin, isbasin) in enumerate(isbasins)
            ax = axs[itracer, ibasin] = fig[itracer, ibasin] = Axis(fig, xaxisposition = :top, xticklabelalign = (:center, :bottom))#, title = titles[itracer])
            xmodelmean, xmodelerror, xobsmean, xobserror, depths = basin_profile_mean_and_std(isbasin, xmodel, xmodelatobs, xdepthatobs, xobs, xwetobs, iobswet)
            xribbon!(ax, xmodelmean, xmodelerror, depths; color=colors[ibasin], αribbon, linewidth=model_thickness, ylims)
            xerrorbars!(ax, xobsmean, xobserror, depths, color=colors[ibasin], linewidth=obs_thickness)
            itracer == 1 && (ax.xlabel = "[Nd] ($(uDNd))")
            itracer == 2 && (ax.xlabel = "εNd ($(uεNd))")
            ibasin == 1 && (ax.ylabel = "Depth (m)")
            Label(fig, bbox = ax.scene.px_area, reshape(panellabels[1:8], (4,2))[ibasin, itracer], textsize=20, halign=:left, valign=:bottom, padding=(10,10,10,10), font=labelfont, color=:black)
            numobs =  count(isbasin(xobs.lat[iobswet], xobs.lon[iobswet], OCEANS))
            Label(fig, bbox = ax.scene.px_area, string(basins[ibasin], " ($numobs obs)"), halign=:right, valign=:bottom, padding=(10,10,10,10), font=labelfont, color=colors[ibasin])
        end
    end
    return axs
end

# Create the plot
axs = plot_profiles_v2!(fig)

# Label axes
#axs[1].xlabel = "[Nd] ($(uDNd))"
#axs[1].ylabel = "Depth (m)"
#axs[2].xlabel = "εNd ($(uεNd))"
#Nd_sublayout = GridLayout()
#ε_sublayout = GridLayout()
#fig[1,:] = Nd_sublayout
#fig[2,:] = ε_sublayout
#Ndlabel = Label(Nd_sublayout[0,:], "[Nd] ($(uDNd))")
#εlabel = Label(ε_sublayout[0,:], "εNd ($(uεNd))")
#leftlabel = Label(fig[1:2, 0], "Depth (m)", rotation=π/2)

# Add customized legend
obs_stl=(color=:black, linewidth=obs_thickness, linestyle=nothing)
group_modelobs = [
    [PolyElement(color=(:black, αribbon), strokecolor=:transparent),
     LineElement(color=:black, linewidth=model_thickness, linestyle=nothing, linepoints = Point2f0[(0.6, 0), (0.4, 1)])],
    [LineElement(;obs_stl...),
     LineElement(;obs_stl..., linepoints = Point2f0[(0.6, 0), (0.4, 1)]),
     LineElement(;obs_stl..., linepoints = Point2f0[(0, 0.35), (0, 0.65)]),
     LineElement(;obs_stl..., linepoints = Point2f0[(1, 0.35), (1, 0.65)])]
]

#basin_element_bot(ibasin) = PolyElement(color=colors[ibasin], strokecolor=:transparent, polypoints=Point2f0[(0, 0), (1, 0), (1, 1)])
#basin_element_top(ibasin) = PolyElement(color=(colors[ibasin], αribbon), strokecolor=:transparent, polypoints=Point2f0[(0, 0), (1, 1), (0, 1)])
#group_basin = [[basin_element_bot(ibasin), basin_element_top(ibasin)] for ibasin in 1:length(basins)]
#leg = Legend(fig, [group_modelobs, group_basin], [["Model", "Observations"], basin_longnames], [" ", "Basin"]) ;
leg = Legend(fig, group_modelobs, ["Model", "Observations"]) ;
leg.orientation = :horizontal
fig[end+1,:] = leg
leg.tellheight = true

# Link y-axes
linkyaxes!(axs...)
linkxaxes!(axs[1,:]...)
linkxaxes!(axs[2,:]...)
[hideydecorations!(ax, grid = false) for ax in axs[:,2:end]]

# Supertitle
#supertitle = layout[0, :] = LText(scene, "Vertical profiles of Nd concentrations & εNd values, basin averaged",
#    textsize = 30, font = "Noto Sans Bold", color = :black)



if use_GLMakie
    display(fig) # show the output wiht GLMakie
else
    save(joinpath(archive_path, "Nd_profiles_exploded_$(lastcommit)_run$(run_num).pdf"), fig)
    nothing # just so that no output is spat out
end

nothing
