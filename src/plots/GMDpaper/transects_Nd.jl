#================================================
Transects
================================================#
include("../plots_setup_Nd.jl")
# Basin profiles where obs only
# Only keep obs part of the cost function
# and interpolate model onto their location
# Function to extract indices of those model boxes (or obs.) that go in comparison




# Full figure with every transect
ts = Nd_transects.transects[Nd_t_sort]
nt = length(ts)
fig = Figure(resolution=(1600, 1600))
#display(s)
axs = Axis[]
scs = Any[]
# colors for cruises
tcol = ColorSchemes.tableau_10.colors
tcol = [:black; tcol]
# Minimap in top left
subl = GridLayout(1,2)
ax = subl[1,1] = Axis(fig, backgroundcolor=water_color)
fig.layout[1,1] = subl
mapit!(ax, clon, mypolys(clon), color=land_color, strokecolor=:black)
cts = Any[]
maptransects!(ax, cts, ts, tcol, wlon)
xlims!(ax, clon .+ (-180,180))
ylims!(ax, (-90,90))
Label(fig, bbox = ax.scene.px_area, panellabels[1], textsize=20, halign=:left, valign=:bottom, padding=(10,10,5,10), font=labelfont, color=labelcol)
hidedecorations!(ax)
hidespines!(ax)
# Legend
idx = permutedims(reshape(1:nt+1, (nt÷2+1, 2)), (2,1))[:]
labels = [" "; [t.cruise for t in ts]][idx]
invisible_element = LineElement(color = :transparent, linewidth = 0, linestyle = nothing)
legelems = [invisible_element; [PolyElement(color=tcol[i], strokecolor=:transparent) for i in eachindex(ts)]][idx]
leg = Legend(fig, legelems, labels, nbanks=2, labelsize=20, framevisible=false)
subl[1,2] = leg
# transect plots
for it in 1:nt
    i, j = Tuple(CartesianIndices((nt÷2+1, 2))[it+1])
    t = sort(OceanographyCruises.shiftlon(ts[it], baselon=wlon))
    ct = CruiseTrack(t)
    local ax = Axis(fig[i,j], xgridvisible=false, ygridvisible=false, backgroundcolor=land_color, halign=:left)
    push!(axs, ax)
    hm = makietransect!(ax, DNdmodel, ct, colormap=Ndcmap, nan_color=nan_color, levels=Ndlevels, colorrange=Ndclims, extendhigh = :auto) # note: added levels with contourf
    sc = makiescattertransect!(ax, t, colormap=Ndcmap, markersize=markersize)
    push!(scs, hm)
    sc.plots[end].colorrange = Ndclims # color of last overlaid scatter
    ylims!(ax, (6000,0))
    dis = ustrip(km, maximum(scattertransect(t)[1]))
    if dis > 1500
        Label(fig, bbox = axs[it].scene.px_area, t.cruise, textsize=20, halign=:right, valign=:bottom, padding=(0,35,5,0), font=labelfont, color=labelcol)
        Label(fig, bbox = axs[it].scene.px_area, panellabels[it+1], textsize=20, halign=:left, valign=:bottom, padding=(10,10,5,10), font=labelfont, color=labelcol)
        Label(fig, bbox = axs[it].scene.px_area, "■", textsize=20, halign=:right, valign=:bottom, padding=(0,10,5,0), font=labelfont, color=tcol[it])
    else # special treatment for tiny transects (labels don't fit)
        Label(fig, bbox = axs[it].scene.px_area, string(panellabels[it+1], " ", t.cruise), textsize=20, halign=:right, valign=:bottom, padding=(0,-115,5,0), font=labelfont, color=labelcol)
        Label(fig, bbox = axs[it].scene.px_area, "■", textsize=20, halign=:right, valign=:bottom, padding=(0,-140,5,0), font=labelfont, color=tcol[it])
    end
    ax.xlabel="Distance (km)"
    ax.ylabel="Depth (m)"
    !((i==nt÷2+1) || (it==nt)) && hidexdecorations!(ax, ticks=false)
    j==2 && hideydecorations!(ax, ticks=false)
    tightlimits!(ax, Bottom())
    ax.xticks = 0:1000:maximum(dis)
    ax.width = 700 * (dis+2t_ext) / 8000
    xlims!(ax, (-t_ext, dis+t_ext))
    ylims!(ax, (6001, -1))
    hidespines!(ax, :t, :r)
end
#linkaxes!(axs...)
# colorbar (need to use scs[1].plots[end] because doubly-overlaid-scatter plot)
cbar = fig[end+1,:] = Colorbar(fig, scs[1].plots[1], vertical=false, width=600, height=30,
                               flipaxis=false, ticklabelalign=(:center, :top), label="[Nd] (pM)",
                               labelsize=25, ticks=Ndlevels[1:5:end])
cbar.alignmode = Outside()
cbar.tellheight = true
# final touches
trim!(fig.layout)

if use_GLMakie
    display(fig) # show the output wiht GLMakie
else
    save(joinpath(archive_path, "Nd_transects_$(lastcommit)_run$(run_num).pdf"), fig)
    nothing # just so that no output is spat out
end
