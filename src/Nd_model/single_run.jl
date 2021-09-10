# This lines sets up the model (and `F`) only once,
# so that you rerunning this file will not repeat the entire model setup
!isdefined(Main, :F) && include("model_setup.jl")

# Chose your parameter values here
# Optimized parameters as published in Pasquier, Hines, et al. (2021)
# are shown as a comment on the right
p = Params(
    α_a = 1.0NoUnits,             # α_a = 1.0NoUnits
    α_c = -10.0εunit,             # α_c = -10.0εunit
    α_GRL = 2.0NoUnits,           # α_GRL = 2.0NoUnits
    σ_ε = 3.0εunit,               # σ_ε = 3.0εunit
    c_river = 100.0pM,            # c_river = 100.0pM
    c_gw = 100.0pM,               # c_gw = 100.0pM
    σ_hydro = 1.0Mmol/yr,         # σ_hydro = 1.0Mmol/yr
    ε_hydro = 10.0εunit,          # ε_hydro = 10.0εunit
    ϕ_0 = 20.0pmol/cm^2/yr,       # ϕ_0 = 20.0pmol/cm^2/yr
    ϕ_∞ = 10.0pmol/cm^2/yr,       # ϕ_∞ = 10.0pmol/cm^2/yr
    z_0 = 200.0m,                 # z_0 = 200.0m
    ε_EAsia_dust = -8.0εunit,     # ε_EAsia_dust = -8.0εunit
    ε_NEAf_dust = -12.0εunit,     # ε_NEAf_dust = -12.0εunit
    ε_NWAf_dust = -12.0εunit,     # ε_NWAf_dust = -12.0εunit
    ε_NAm_dust = -8.0εunit,       # ε_NAm_dust = -8.0εunit
    ε_SAf_dust = -10.0εunit,      # ε_SAf_dust = -10.0εunit
    ε_SAm_dust = -3.0εunit,       # ε_SAm_dust = -3.0εunit
    ε_MECA_dust = -2.0εunit,      # ε_MECA_dust = -2.0εunit
    ε_Aus_dust = -4.0εunit,       # ε_Aus_dust = -4.0εunit
    ε_Sahel_dust = -12.0εunit,    # ε_Sahel_dust = -12.0εunit
    β_EAsia_dust = 5.0u"percent", # β_EAsia_dust = 5.0u"percent"
    β_NEAf_dust = 5.0u"percent",  # β_NEAf_dust = 5.0u"percent"
    β_NWAf_dust = 5.0u"percent",  # β_NWAf_dust = 5.0u"percent"
    β_NAm_dust = 5.0u"percent",   # β_NAm_dust = 5.0u"percent"
    β_SAf_dust = 5.0u"percent",   # β_SAf_dust = 5.0u"percent"
    β_SAm_dust = 5.0u"percent",   # β_SAm_dust = 5.0u"percent"
    β_MECA_dust = 5.0u"percent",  # β_MECA_dust = 5.0u"percent"
    β_Aus_dust = 5.0u"percent",   # β_Aus_dust = 5.0u"percent"
    β_Sahel_dust = 5.0u"percent", # β_Sahel_dust = 5.0u"percent"
    ε_volc = 10.0εunit,           # ε_volc = 10.0εunit
    β_volc = 10.0u"percent",      # β_volc = 10.0u"percent"
    K_prec = 0.01NoUnits,         # K_prec = 0.01NoUnits
    f_prec = 0.4NoUnits,          # f_prec = 0.4NoUnits
    w₀_prec = 0.7km/yr,           # w₀_prec = 0.7km/yr
    K_POC = 3e13NoUnits,          # K_POC = 3e13NoUnits
    f_POC = 0.78NoUnits,          # f_POC = 0.78NoUnits
    w₀_POC = 40.0m/d,             # w₀_POC = 40.0m/d
    K_bSi = 3e13NoUnits,          # K_bSi = 3e13NoUnits
    f_bSi = 0.5NoUnits,           # f_bSi = 0.5NoUnits
    w₀_bSi = 714.069m/d,          # w₀_bSi = 714.069m/d
    K_dust = 2e15NoUnits,         # K_dust = 2e15NoUnits
    f_dust = 0.073NoUnits,        # f_dust = 0.073NoUnits
    w₀_dust = 1.0km/yr            # w₀_dust = 1.0km/yr
)

# Set the problem with the parameters above
prob = SteadyStateProblem(fun, x, p)

# solve the system
sol = solve(prob, CTKAlg(), preprint="Nd & εNd solve ", τstop=ustrip(u"s", 1e3u"Myr"))

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