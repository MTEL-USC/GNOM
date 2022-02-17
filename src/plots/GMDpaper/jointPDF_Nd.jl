#================================================
Joint PDF
================================================#
include("../plots_setup_Nd.jl")

# Create the figure
fig = Figure(resolution=(1200, 700))
use_GLMakie && display(fig)

Nd_bounds = (0.0, 40.0)
εNd_bounds = (-20.0, 5.0)

tracer_names = ("[Nd]", "εNd")
labels = ("(a)", "(b)")

cmap = cgrad(:oslo, 10, categorical=true, rev=true)
cmap2 = cgrad(cmap[2:end], categorical=true) # to skip the white color to leave the grid behind (prettier)

function myjointpdf!(fig)
    axs = Array{Any}(undef, 2)
    cos = []
    for (itracer, (xmodel, xobs, ux, boundary, tracer_name, label)) in enumerate(zip((DNdmodel, εNdmodel), (DNdobs, εNdobs), (uDNd, uεNd), (Nd_bounds, εNd_bounds), tracer_names, labels))
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

        ax = fig[1, itracer] = Axis(fig, aspect = AxisAspect(1))
        co = contourf!(ax, D.x, D.y, Dcum, levels=10:10:100, colormap=cmap2)
        lines!(ax, collect(boundary), collect(boundary), linestyle=:dash, color=:black)
        #scatter!(ax, x, y, markersize=1)
        xlims!(ax, boundary)
        ylims!(ax, boundary)
        ax.xlabel = "Observed " * tracer_name * " ($ux)"
        ax.ylabel = "Modelled " * tracer_name * " ($ux)"

        # Add label
        Label(fig, bbox = ax.scene.px_area, label, textsize=20, halign=:left, valign=:top, padding=(10,10,5,5), font=labelfont, color=:black)

        # Root mean square error
        RMSE = sqrt(mean((x - y).^2))
        Label(fig, bbox = ax.scene.px_area, sprintf1("RMSE = %.2f $ux", RMSE), textsize=15, halign=:right, valign=:bottom, padding=(10,10,10,10), font=labelfont, color=:black)

        push!(axs, ax)
        push!(cos, co)
    end
    axs

    # colorbar
    cbar = fig[2, 1:2] = Colorbar(fig, #cos[1],
                                  colormap=cmap, ticks = 0:10:100, limits=(0,100),
                                  label="Percentile", vertical=false, flipaxis=false, ticklabelalign = (:center, :top))
    cbar.width = Relative(2/4)
    cbar.height = 30
    cbar.tellheight = true

    # make top plots square
    sublayout = GridLayout()
    fig[1, 1:2] = sublayout
    rowsize!(sublayout, 1, Aspect(1,1))


    fig
end

# Create the plot
myjointpdf!(fig)

if use_GLMakie
    fig # show the output wiht GLMakie
else
    save(joinpath(archive_path, "jointPDF_$(lastcommit)_run$(run_num).pdf"), fig)
    nothing # just so that no output is spat out
end



function myjointpdf2!(fig)
    axs = Array{Any}(undef, 2)
    cos = []
    for (itracer, (xmodel, xobs, ux, boundary, tracer_name, label)) in enumerate(zip((DNdmodel, εNdmodel), (DNdobs, εNdobs), (uDNd, uεNd), (Nd_bounds, εNd_bounds), tracer_names, labels))
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

        ax = fig[itracer, 1] = Axis(fig, aspect = AxisAspect(1))
        co = contourf!(ax, D.x, D.y, Dcum, levels=10:10:100, colormap=cmap2)
        lines!(ax, collect(boundary), collect(boundary), linestyle=:dash, color=:black)
        #scatter!(ax, x, y, markersize=1)
        xlims!(ax, boundary)
        ylims!(ax, boundary)
        ax.xlabel = "Observed " * tracer_name * " ($ux)"
        ax.ylabel = "Modelled " * tracer_name * " ($ux)"

        # Add label
        Label(fig, bbox = ax.scene.px_area, label, textsize=20, halign=:left, valign=:top, padding=(10,10,5,5), font=labelfont, color=:black)

        # Root mean square error
        RMSE = sqrt(mean((x - y).^2))
        Label(fig, bbox = ax.scene.px_area, sprintf1("RMSE = %.2f $ux", RMSE), textsize=15, halign=:right, valign=:bottom, padding=(10,10,10,10), font=labelfont, color=:black)

        push!(axs, ax)
        push!(cos, co)
    end
    axs

    # colorbar
    cbar = fig[3, 1] = Colorbar(fig, #cos[1],
                                  colormap=cmap, ticks = 0:10:100, limits=(0,100),
                                  label="Percentile", vertical=false, flipaxis=false, ticklabelalign = (:center, :top))
    cbar.width = Relative(3/4)
    cbar.height = 30
    cbar.tellheight = true

    # make top plots square
    sublayout = GridLayout()
    fig[1:2, 1] = sublayout
    colsize!(sublayout, 1, Aspect(1,1))


    fig
end

# Create the plot
fig = Figure(resolution=(450, 1000))
use_GLMakie && display(fig)
myjointpdf2!(fig)

if use_GLMakie
    display(fig) # show the output wiht GLMakie
else
    save(joinpath(archive_path, "jointPDF_vertical_$(lastcommit)_run$(run_num).pdf"), fig)
    nothing # just so that no output is spat out
end
