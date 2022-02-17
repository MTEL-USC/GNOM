# You are running a single run
# that will create DNd, εNd, and p.

# This lines sets up the model (and `F`) only once,
# so that you rerunning this file will not repeat the entire model setup
!isdefined(Main, :F) && include("model_setup.jl")

# This should create a new run name every time you run the file
# And add an empty runXXX file to the single_runs folder

allsingleruns_path = joinpath(output_path, "single_runs")
mkpath(output_path)
mkpath(allsingleruns_path)
# Check previous runs and get new run number
run_num = let
    previous_run_nums = [parse(Int, match(r"run(\d+)", f).captures[1]) for f in readdir(allsingleruns_path) if (contains(f, "run") && isdir(joinpath(allsingleruns_path, f)))]
    run_num = 1
    while run_num ∈ previous_run_nums
        run_num += 1
    end
    run_num
end
@info "This is run $run_num"
lastcommit = "single" # misnomer to call this lastcommit but simpler
archive_path = joinpath(allsingleruns_path, "run$run_num")
mkpath(archive_path)
reload = false # prevents loading other runs
use_GLMakie = true # Set to true for interactive mode if plotting with Makie later


# Chose your parameter values here. Optimized parameters
# — as published in Pasquier, Hines, et al. (2021) —
# are shown in comment (leave them there
# if you want to refer back to them)
p = Params(
    α_a = 6.79, # 6.79
    α_c = -12.7per10000, # -12.7per10000
    α_GRL = 1.57, # 1.57
    σ_ε = 0.379per10000, # 0.379per10000
    c_river = 376.0pM, # 376.0pM
    c_gw = 109.0pM, # 109.0pM
    σ_hydro = 0.792Mmol/yr, # 0.792Mmol/yr
    ε_hydro = 10.9per10000, # 10.9per10000
    ϕ_0 = 83.7pmol/cm^2/yr, # 83.7pmol/cm^2/yr
    ϕ_∞ = 1.11pmol/cm^2/yr, # 1.11pmol/cm^2/yr
    z_0 = 170.0m, # 170.0m
    ε_EAsia_dust = -7.6per10000, # -7.6per10000
    ε_NEAf_dust = -13.7per10000, # -13.7per10000
    ε_NWAf_dust = -12.3per10000, # -12.3per10000
    ε_NAm_dust = -4.25per10000, # -4.25per10000
    ε_SAf_dust = -21.6per10000, # -21.6per10000
    ε_SAm_dust = -3.15per10000, # -3.15per10000
    ε_MECA_dust = 0.119per10000, # 0.119per10000
    ε_Aus_dust = -4.03per10000, # -4.03per10000
    ε_Sahel_dust = -11.9per10000, # -11.9per10000
    β_EAsia_dust = 23.0per100, # 23.0per100
    β_NEAf_dust = 23.3per100, # 23.3per100
    β_NWAf_dust = 3.17per100, # 3.17per100
    β_NAm_dust = 82.8per100, # 82.8per100
    β_SAf_dust = 38.5per100, # 38.5per100
    β_SAm_dust = 2.52per100, # 2.52per100
    β_MECA_dust = 14.7per100, # 14.7per100
    β_Aus_dust = 11.6per100, # 11.6per100
    β_Sahel_dust = 2.95per100, # 2.95per100
    ε_volc = 13.1per10000, # 13.1per10000
    β_volc = 76.0per100, # 76.0per100
    K_prec = 0.00576/(mol/m^3), # 0.00576/(mol/m^3)
    f_prec = 0.124, # 0.124
    w₀_prec = 0.7km/yr, # 0.7km/yr
    K_POC = 0.524/(mol/m^3), # 0.524/(mol/m^3)
    f_POC = 0.312, # 0.312
    w₀_POC = 40.0m/d, # 40.0m/d
    K_bSi = 2.56/(mol/m^3), # 2.56/(mol/m^3)
    f_bSi = 0.784, # 0.784
    w₀_bSi = 714.0m/d, # 714.0m/d
    K_dust = 1.7/(g/m^3), # 1.7/(g/m^3)
    f_dust = 0.0861, # 0.0861
    w₀_dust = 1.0km/yr, # 1.0km/yr
)

tp_opt = AIBECS.table(p)# table of parameters
# "opt" is a misnomer but it is simpler for plotting scripts

# Save model parameters table and headcommit for safekeeping
jldsave(joinpath(archive_path, "model$(headcommit)_single_run$(run_num)_$(circname).jld2"); headcommit, tp_opt)

# Set the problem with the parameters above
prob = SteadyStateProblem(fun, x, p)

# solve the system
sol = solve(prob, CTKAlg(), preprint="Nd & εNd solve ", τstop=ustrip(s, 1e3Myr)).u

# unpack nominal isotopes
DNd, DRNd = unpack_tracers(sol, grd)

# compute εNd
εNd = ε.(DRNd ./ DNd)

# For plotting, you can either
# follow the plotting scripts from the GNOM repository and use Makie
# or use Plots.jl (not a dependency of GNOM)
# I would recommend installing Plots.jl in your default environment anyway,
# so that it can be called even from inside the GNOM environment.
# You can then use the Plots.jl recipes exported by AIBECS, e.g.,
#
# julia> plotzonalaverage(εNd .|> per10000, grd, mask=ATL)