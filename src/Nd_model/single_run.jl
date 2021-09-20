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
# are shown in comment
p = Params(
    α_a = 4.25896,                   # α_a = 2.52
    α_c = -13.4003εunit,             # α_c = -11.6εunit
    α_GRL = 2.45502,                 # α_GRL = 1.03
    σ_ε = 0.255468εunit,               # σ_ε = 3.0εunit
    c_river = 1022.22pM,            # c_river = 100.0pM
    c_gw = 60.7159pM,               # c_gw = 100.0pM
    σ_hydro = 0.452881Mmol/yr,         # σ_hydro = 1.0Mmol/yr
    ε_hydro = 5.3727εunit,          # ε_hydro = 10.0εunit
    ϕ_0 = 100.613pmol/cm^2/yr,       # ϕ_0 = 20.0pmol/cm^2/yr
    ϕ_∞ = 3.34434pmol/cm^2/yr,       # ϕ_∞ = 10.0pmol/cm^2/yr
    z_0 = 159.749m,                 # z_0 = 200.0m
    ε_EAsia_dust = -7.29218εunit,     # ε_EAsia_dust = -8.0εunit
    ε_NEAf_dust = -12.1979εunit,     # ε_NEAf_dust = -12.0εunit
    ε_NWAf_dust = -12.2858εunit,     # ε_NWAf_dust = -12.0εunit
    ε_NAm_dust = -9.17992εunit,       # ε_NAm_dust = -8.0εunit
    ε_SAf_dust = -9.36719εunit,      # ε_SAf_dust = -10.0εunit
    ε_SAm_dust = -3.42779εunit,       # ε_SAm_dust = -3.0εunit
    ε_MECA_dust = -1.74237εunit,      # ε_MECA_dust = -2.0εunit
    ε_Aus_dust = -3.73625εunit,       # ε_Aus_dust = -4.0εunit
    ε_Sahel_dust = -12.1482εunit,    # ε_Sahel_dust = -12.0εunit
    β_EAsia_dust = 2*0.590052u"percent", # β_EAsia_dust = 5.0u"percent"
    β_NEAf_dust = 2*1.04977u"percent",  # β_NEAf_dust = 5.0u"percent"
    β_NWAf_dust = 2*2.8137u"percent",  # β_NWAf_dust = 5.0u"percent"
    β_NAm_dust = 2*1.11589u"percent",   # β_NAm_dust = 5.0u"percent"
    β_SAf_dust = 2*0.566935u"percent",   # β_SAf_dust = 5.0u"percent"
    β_SAm_dust = 2*1.42008u"percent",   # β_SAm_dust = 5.0u"percent"
    β_MECA_dust = 2*1.80498u"percent",  # β_MECA_dust = 5.0u"percent"
    β_Aus_dust = 2*1.11362u"percent",   # β_Aus_dust = 5.0u"percent"
    β_Sahel_dust = 2*4.22232u"percent", # β_Sahel_dust = 5.0u"percent"
    ε_volc = 10.1533εunit,           # ε_volc = 10.0εunit
    β_volc = 2*4.03051u"percent",      # β_volc = 10.0u"percent"
    K_prec = 0.00848695,         # K_prec = 0.01
    f_prec = 0.0695125,          # f_prec = 0.4
    w₀_prec = 0.7km/yr,           # w₀_prec = 0.7km/yr
    K_POC = 1.07351e14,          # K_POC = 3e13
    f_POC = 0.205732,          # f_POC = 0.78
    w₀_POC = 40.0m/d,             # w₀_POC = 40.0m/d
    K_bSi = 1.17098e12,          # K_bSi = 3e13
    f_bSi = 0.711062,           # f_bSi = 0.5
    w₀_bSi = 714.069m/d,          # w₀_bSi = 714.069m/d
    K_dust = 9.49175e14,         # K_dust = 2e15
    f_dust = 0.0607855,        # f_dust = 0.073
    w₀_dust = 1.0km/yr            # w₀_dust = 1.0km/yr
)

tp_opt = AIBECS.table(p)# table of parameters
# "opt" is a misnomer but it is simpler for plotting scripts

# Save model parameters table and headcommit for safekeeping
jldsave(joinpath(archive_path, "model$(headcommit)_single_run$(run_num)_$(circname).jld2"); headcommit, tp_opt)

# Set the problem with the parameters above
prob = SteadyStateProblem(fun, x, p)

# solve the system
sol = solve(prob, CTKAlg(), preprint="Nd & εNd solve ", τstop=ustrip(u"s", 1e3u"Myr")).u

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
# julia> plotzonalaverage(εNd .|> εunit, grd, mask=ATL)