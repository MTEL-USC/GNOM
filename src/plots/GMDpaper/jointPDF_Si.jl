#================================================
Joint PDF
================================================#
include("../plots_setup_Nd.jl")

DSi, tp_opt_Si = jldopen(joinpath(output_path, "optimized_Simodel_$circname.jld2")) do f
    f["DSi"], f["tp_opt"]
end
uDSi = mmol/m^3
DSimodel = uconvert.(uDSi, DSi * upreferred(uDSi))
using WorldOceanAtlasTools
DSiobs = let
    obs = WorldOceanAtlasTools.observations("silicate")
    obs.value = ustrip.(upreferred.(obs.silicate * ρSW))
    obs
end

cmap = cgrad(:lajolla, 10, categorical=true)
cmap2 = cgrad(cmap[2:end], categorical=true) # to skip the white color to leave the grid behind (prettier)

function myjointpdf_Si!(fig)
    xmodel = DSimodel
    xobs = DSiobs
    ux = uDSi
    boundary = (-5.0, 180.0)
    tracer_name = "[DSi]"

    xmodelatobs, xdepthatobs, iobswet, xwetobs = _locations(xmodel, xobs, ux)
    #@show size.([xwetobs, xmodelatobs])
    x, y = ustrip.(xwetobs), ustrip.(xmodelatobs)
    bw = (boundary[2]-boundary[1])/150
    D = kde((x, y); boundary=(boundary, boundary), bandwidth=(bw,bw))

    # calculate cumulative density from density
    δx = step(D.x)
    δy = step(D.y)
    Q = vec(D.density) * δx * δy
    idx = sortperm(Q)
    Q_sorted = Q[idx]
    Dcum = similar(D.density)
    Dcum[idx] .= 100cumsum(Q_sorted)

    ax = fig[1, 1] = Axis(fig, aspect = AxisAspect(1))
    co = contourf!(ax, D.x, D.y, Dcum, levels=10:10:100, colormap=cmap2)
    lines!(ax, collect(boundary), collect(boundary), linestyle=:dash, color=:black)
    #scatter!(ax, x, y, markersize=1)
    xlims!(ax, boundary)
    ylims!(ax, boundary)
    ax.xlabel = "Observed " * tracer_name * " ($ux)"
    ax.ylabel = "Modelled " * tracer_name * " ($ux)"

    # colorbar
    cbar = fig[2, 1] = Colorbar(fig, #cos[1],
                                  colormap=cmap, ticks = 0:10:100, limits=(0,100),
                                  label="Percentile", vertical=false, flipaxis=false, ticklabelalign = (:center, :top))
    cbar.width = Relative(3/4)
    cbar.height = 20
    cbar.tellheight = true

    # make top plots square
    sublayout = GridLayout()
    fig[1:1, 1] = sublayout
    colsize!(sublayout, 1, Aspect(1,1))

    # Root mean square error
    RMSE = sqrt(mean((x - y).^2))
    Label(fig, bbox = ax.scene.px_area, sprintf1("RMSE = %.1f $ux", RMSE), textsize=12, halign=:right, valign=:bottom, padding=(10,10,10,10), font=labelfont, color=:black)

    fig
end

# Create the plot
fig = Figure(resolution=(450, 450))
myjointpdf_Si!(fig)

#save(joinpath(archive_path, "jointPDF_$(lastcommit)_run$(run_num).png"), fig, px_per_unit=4)
if use_GLMakie
    display(fig) # show the output wiht GLMakie
else
    save(joinpath(archive_path, "jointPDF_Si_$(lastcommit)_run$(run_num).pdf"), fig)
    nothing # just so that no output is spat out
end


#======================================#
#          Print param tables          #
#======================================#
# Don't reload eveything everytime
tpSi = select(tp_opt_Si,
             :Symbol=>ByRow(symbol2latex)=>:Symbol,
             :Value,
             Symbol("Initial value"),
             :Prior=>ByRow(latexify)=>:Range,
             :Unit=>ByRow(latexify)=>:Unit,
             :Description=>ByRow(latexeNd)=>:Description)
             #:Optimizable=>ByRow(latexbool)=>:Optimized)
formatters = (v,i,j) -> j ∈ [2,3] ? string("\$", numformat(sprintf1("%.3g", v)), "\$") : (j == 4) ? string("\$", v, "\$") : v

println("Latex param table")

println(pretty_table(tpSi, tf=tf_latex_simple, formatters=formatters, nosubheader=true))
