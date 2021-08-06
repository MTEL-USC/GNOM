#================================================
Profiles
================================================#
use_GLMakie = false
include("../plots_setup_Nd.jl")

# Create the figure
fig = Figure(resolution=(1200, 1500))
use_GLMakie && display(fig)

# α parameters
for (i,s) in enumerate((:α_a, :α_c, :σ_ε, :α_GIC))
    local ax = fig[1,i] = Axis(fig)
    plot_params!(ax, p, s)
end
# river, gw, hydro params
for (i,s) in enumerate((:sol_river, :sol_gw, :σ_hydro, :ε_hydro))
    local ax = fig[2,i] = Axis(fig)
    plot_params!(ax, p, s)
end
# sed flux params
for (i,s) in enumerate((:ϕ_0, :ϕ_∞, :z_0))
    local ax = fig[3,i] = Axis(fig)
    plot_params!(ax, p, s)
end
# dust params 2×9 to 3×6
dust_params = (:ε_EAsia_dust, :ε_NEAf_dust, :ε_NWAf_dust, :ε_NAm_dust, :ε_SAf_dust, :ε_SAm_dust, :ε_MECA_dust, :ε_Aus_dust, :ε_Sahel_dust, :sol_EAsia_dust, :sol_NEAf_dust, :sol_NWAf_dust, :sol_NAm_dust, :sol_SAf_dust, :sol_SAm_dust, :sol_MECA_dust, :sol_Aus_dust, :sol_Sahel_dust)
for (i,s) in enumerate(dust_params)
    I, J = Tuple(CartesianIndices((3,6))[i])
    local ax = fig[3+I,J] = Axis(fig)
    plot_params!(ax, p, s)
end
# volc params
for (i,s) in enumerate((:sol_volc, :ε_volc))
    local ax = fig[7,i] = Axis(fig)
    plot_params!(ax, p, s)
end
# scavenging params 4×2 so 2×4
scav_params = (:K_prec, :f_prec, :K_POC, :f_POC, :K_bSi, :f_bSi, :K_dust, :f_dust)
for (i,s) in enumerate(scav_params)
    I, J = Tuple(CartesianIndices((2,4))[i])
    local ax = fig[7+I,J] = Axis(fig)
    plot_params!(ax, p, s)
end
ax = fig[10,1] = Axis(fig)
plot_params!(ax, p, :τ_ns)

save(joinpath(archive_path, "parameters_$(lastcommit)_run$(run_num).pdf"), fig)

