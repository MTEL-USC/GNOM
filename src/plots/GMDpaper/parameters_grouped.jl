#===============================================#
#              Parameters plots                 #
#===============================================#
include("../plots_setup_Nd.jl")

# Create the figure
fig = Figure(resolution=(1200, 1500))


# layout of params:
# a) 10×2 ε and β for dust
# b) 1×2 ε and β for volc
# c) α curve
# d) river source ε and concentration
# e) hydro source ε and concentration
# f) sed. flux ϕ
# g) scavening K and f

#     +-------+-------+----+
#     |       |       |    |
#     |       |   c   |    |
#     |       |       |    |
#     |       +-------| f  |
#     |       |   d   |    |
#     |       +-------|    |
#     |       |   e   |    |
#     |   a   +------------+
#     |       |            |
#     |       |            |
#     |       |            |
#     |       |     g      |
#     +-------|            |
#     |   b   |            |
#     +-------+------------+

gab = fig[1, 1:2] = GridLayout()
ga = gab[1:9, 1] = GridLayout()
gb = gab[10, 1] = GridLayout()
gcg = fig[1, 3:5] = GridLayout()
gcf = gcg[1, 1] = GridLayout()
gg = gcg[2, 1] = GridLayout()
gce = gcf[1, 1:2] = GridLayout()
gf = gcf[1, 3] = GridLayout()
gc = gce[1:2, 1] = GridLayout()
gd = gce[3, 1] = GridLayout()
ge = gce[4, 1] = GridLayout()

# colors of densities

# a) 9×2 ε and β for dust and volc
ε_params = (:ε_EAsia_dust, :ε_NEAf_dust, :ε_NWAf_dust, :ε_NAm_dust, :ε_SAf_dust, :ε_SAm_dust, :ε_MECA_dust, :ε_Aus_dust, :ε_Sahel_dust)
β_params = (:β_EAsia_dust, :β_NEAf_dust, :β_NWAf_dust, :β_NAm_dust, :β_SAf_dust, :β_SAm_dust, :β_MECA_dust, :β_Aus_dust, :β_Sahel_dust)
for (i,s) in enumerate(ε_params)
    local ax = ga[i,1] = Axis(fig)
    plot_param!(ax, p, s; density_color=density_colors[1])
    xlims!(ax, (-35, 15))
    ax.ylabel = split(string(s), "_")[2]
    hideydecorations!(ax, label=false)
    i≠9 ? hidexdecorations!(ax, grid=false) : ax.xlabel="ε (‱)"
end
for (i,s) in enumerate(β_params)
    local ax = ga[i,2] = Axis(fig)
    plot_param!(ax, p, s; density_color=density_colors[1])
    xlims!(ax, (0,100))
    hideydecorations!(ax)
    ax.xticks=0:20:100
    i≠9 ? hidexdecorations!(ax, grid=false) : ax.xlabel="β (%)"
end

# b) volcanic ash params
let
    local ax = gb[1,1] = Axis(fig)
    plot_param!(ax, p, :ε_volc; density_color=density_colors[2])
    xlims!(ax, (-35, 15))
    ax.xlabel="ε (‱)"
    ax.ylabel = "volc"
    hideydecorations!(ax, label=false)
end
let
    local ax = gb[1,2] = Axis(fig)
    plot_param!(ax, p, :β_volc; density_color=density_colors[2])
    xlims!(ax, (0,100))
    hideydecorations!(ax)
    ax.xticks=0:20:100
    ax.xlabel="β (%)"
end

# c) α curve
for (i,s) in enumerate((:α_a, :α_c, :σ_ε, :α_GRL))
    I, J = Tuple(CartesianIndices((2,2))[i])
    local ax = gc[I,J] = Axis(fig)
    plot_param!(ax, p, s; density_color=density_colors[3])
    hideydecorations!(ax)
end

# d) river source ε and concentration
for (i,s) in enumerate((:c_river, :c_gw))
    local ax = gd[1,i] = Axis(fig)
    plot_param!(ax, p, s; density_color=density_colors[4])
    hideydecorations!(ax)
end

# e) hydro source ε and concentration
for (i,s) in enumerate((:σ_hydro, :ε_hydro))
    local ax = ge[1,i] = Axis(fig)
    plot_param!(ax, p, s; density_color=density_colors[5])
    hideydecorations!(ax)
end

# f) sed. flux ϕ
for (i,s) in enumerate((:ϕ_0, :ϕ_∞, :z_0))
    local ax = gf[i,1] = Axis(fig)
    plot_param!(ax, p, s; density_color=density_colors[6])
    hideydecorations!(ax)
end

# g) scavening K and f
scav_params = (:K_prec, :K_POC, :K_bSi, :K_dust,
               :f_prec, :f_POC, :f_bSi, :f_dust)
for (i,s) in enumerate(scav_params)
    I, J = Tuple(CartesianIndices((4,2))[i])
    local ax = gg[I,J] = Axis(fig)
    powers_of_ten = round(Int, log10(initial_value(p, s)))
    plot_param!(ax, p, s, powers_of_ten; density_color=density_colors[7])
    ax.ylabel = split(string(s), "_")[2]
    hideydecorations!(ax, label=J≠1)
end

# subplot labels
named_labels = ["Dust", "Volcanic ash", "Enhanced Nd release", "Riverine source", "Hyrdothermal source", "ϕ(z)", "Scavenging"]
for (label, layout, nlab) in zip(panellabels, [ga, gb, gc, gd, ge, gf, gg], named_labels)
    Label(layout[1, :, Top()], nlab * " parameters", valign = :bottom, padding = (0, 0, 5, 0))
    Label(layout[1, 1, TopLeft()], label,
        textsize = 18,
        font = labelfont,
        padding = (0, 5, 5, 0),
        halign = :right)
end

if use_GLMakie
    display(fig) # show the output wiht GLMakie
else
    save(joinpath(archive_path, "parameters_grouped_$(lastcommit)_run$(run_num).pdf"), fig)
    nothing # just so that no output is spat out
end


#======================================#
#          Print param tables          #
#======================================#
# Don't reload eveything everytime
tp2 = select(tp_opt,
             :Symbol=>ByRow(symbol2latex)=>:Symbol,
             :Value,
             Symbol("Initial value"),
             :Prior=>ByRow(latexify)=>:Range,
             :Unit=>ByRow(latexify)=>:Unit,
             :Description=>ByRow(latexeNd)=>:Description)
             #:Optimizable=>ByRow(latexbool)=>:Optimized)
formatters = (v,i,j) -> j ∈ [2,3] ? string("\$", numformat(sprintf1("%.3g", v)), "\$") : (j == 4) ? string("\$", v, "\$") : v

println("Latex param table")

println(pretty_table(tp2, tf=tf_latex_simple, formatters=formatters, nosubheader=true))
#open(joinpath(archive_path, "optimized_parameters.tex"), "w") do f
#    pretty_table(f, tp2, tf=tf_latex_simple, formatters=formatters, nosubheader=true)
#end

