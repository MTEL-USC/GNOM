# This is a simple nutrient restoring model for the global Si cycle.
# The goal here is to generate a reasonable 3D field for particulate Si (PSi)
# concentrations (AKA opal, or biogenic silica), in order to be used as
# scavenging sites for the main Nd-cycling model.
#
# The Si model has a few parameters that are optimized to ensure a decent fit
# to World Ocean Atlas (WOA) observations.
#
# This optimization must be run prior to running the main Nd model.
# This script generates the PSi field and saves it in the `output/`
# directory, which the Nd model expects.

include("environment_setup.jl")

using Unitful: g, mg, kg, L, mol, mmol, μmol, pmol, Mmol, yr, kyr, Myr, d, m, km, cm

#============================
Circulation grid and operator
============================#
const Circulation = OCIM2
const grd, T = Circulation.load()
const nb = count(iswet(grd))

#===============
Model parameters
===============#
import AIBECS: @units, units
import AIBECS: @limits, limits
import AIBECS: @initial_value, initial_value
import AIBECS: @flattenable, flattenable, flatten
import AIBECS: @description, description
const ∞ = Inf
@initial_value @units @flattenable @limits @description struct SiParams{Tp} <: AbstractParameters{Tp}
    f::Tp       | 0.99  | NoUnits  | true  | (0,1)  | "Fraction of non-buried bSi"
    w::Tp       | 200.0 | m/d      | true  | (0,∞)  | "Settling velocity of bSi" 
    τ_up::Tp    | 30.0  | d        | true  | (0,∞)  | "Si(OH)₄ restoring timescale" 
    τ_rem::Tp   |  1.0  | d        | true  | (0,∞)  | "bSi remineralization timescale" 
    DSi_geo::Tp | 80.0  | mmol/m^3 | true  | (0,∞)  | "Si(OH)₄ geological restoring target" 
    τ_geo::Tp   | 1.0   | Myr      | false | (0,∞)  | "Si(OH)₄ geological restoring timescale" 
    z₀::Tp      | 80.0  | m        | false | (0,∞)  | "Depth of euphotic zone" 
end
# Assign a prior to each parameter based on initial value and limits
import AIBECS: @prior, prior
function prior(::Type{T}, s::Symbol) where {T<:AbstractParameters}
    if flattenable(T, s)
        lb, ub = limits(T, s)
        if (lb, ub) == (0,∞)
            μ = log(initial_value(T, s))
            LogNormal(μ, 1.0)
        elseif (lb, ub) == (-∞,∞)
            μ = initial_value(T, s)
            σ = 10.0 # Assumes that a sensible unit is chosen (i.e., that within 10.0 * U)
            Distributions.Normal(μ, σ)
        else
            LocationScale(lb, ub - lb, LogitNormal())
        end
    else
        nothing
    end
end
prior(::T, s::Symbol) where {T<:AbstractParameters} = prior(T,s)
prior(::Type{T}) where {T<:AbstractParameters} = Tuple(prior(T,s) for s in AIBECS.symbols(T))
prior(::T) where {T<:AbstractParameters} = prior(T)

# BGC functions
relu(x) = (x ≥ 0) * x
const z = depthvec(grd)
function uptake(DSi, p)
    @unpack τ_up, z₀ = p
    @. (z < z₀) * relu(DSi - Si_obs) / τ_up
end
function remin(PSi, p)
    @unpack τ_rem = p
    @. PSi / τ_rem
end
function geores(DSi, p)
    @unpack DSi_geo, τ_geo = p
    @. (DSi_geo - DSi) / τ_geo
end

# Problem setup
T_DSi(p) = T
const z_bot = bottomdepthvec(grd)
const z_top = topdepthvec(grd)
const frac_seafloor = ETOPO.fractiontopo(grd)
function T_PSi(p)
    @unpack w, f = p
    transportoperator(grd, z -> w ; z_top, z_bot, frac_seafloor, fsedremin=f) + T
end


G_DSi(DSi, PSi, p) = geores(DSi, p) - uptake(DSi, p) + remin(PSi, p)
G_PSi(DSi, PSi, p) = uptake(DSi, p) - remin(PSi, p)

# Start from previous BSON-saved run close to optimal values
# to run locally quicker to get an optimized 3D field because
# BSON file can't be used anymore (DataFrames update caused this?)
# and the USC HPC is taking ages with Pkg (so can't run remotely).
p = SiParams(f=0.5, w=700.0, τ_rem=3.5, DSi_geo=3e4)

# initial guess
x0 = ustrip(upreferred(10mmol/m^3)) * ones(2nb)

# state function and its Jacobian
#F, ∇ₓF = state_function_and_Jacobian(T_D, Gs, nb)
fun = AIBECSFunction((T_DSi, T_PSi), (G_DSi, G_PSi), nb, SiParams)
F, ∇ₓF = F_and_∇ₓF(fun)

# problem
prob = SteadyStateProblem(fun, x0, p)








# Optimize model
s = solve(prob, CTKAlg(), preprint="Nd & εNd solve ", τstop=ustrip(u"s", 1e3u"Myr"))

#================================================
Optimization setup
================================================#
# Take WOA18 for Si obs
const ρSW = 1.035kg/L
Si_obs_3D, σ²Si_obs_3D = WorldOceanAtlasTools.fit_to_grid(grd, 2018, "silicate", 0, "1°", "an")
Si_obs_3D_painted = inpaint(Si_obs_3D, 0)
Si_obs = Si_obs_3D_painted[iswet(grd)] * upreferred(mol/g) * ρSW .|> upreferred .|> ustrip

# Check obs
#plothorizontalslice(Si_obs, grd, depth=0)

# Model steup
# Two tracers, Si(OH)₄ + PSi



# isotope function
modify(DSi, PSi) = (DSi,)

# Weights for objective function
const ωs = (1.0,) # the weight for the mismatches (note the new API reduced to observed parts)
const ωp = 1e-4
const Si_obs_tuple = let
    obs = WorldOceanAtlasTools.observations("silicate")
    obs.value = ustrip.(upreferred.(obs.silicate * ρSW))
    (obs,)
end

f, ∇ₓf = f_and_∇ₓf(ωs, ωp, grd, modify, Si_obs_tuple, SiParams)

# Use F1 for gradient and Hessian
λ = p2λ(p)
τstop = ustrip(u"s", 1e3u"Myr")
mem = F1Method.initialize_mem(F, ∇ₓf, ∇ₓF, x0, λ, CTKAlg(); preprint="mem ", τstop=τstop)

function objective(λ)
    p = λ2p(SiParams, λ) ; @show p
    F1Method.objective(f, F, ∇ₓF, mem, λ, CTKAlg(), preprint="obj ", τstop=τstop)
end
gradient(λ) = F1Method.gradient(f, F, ∇ₓf, ∇ₓF, mem, λ, CTKAlg(), preprint="grad", τstop=τstop)
hessian(λ) = F1Method.hessian(f, F, ∇ₓf, ∇ₓF, mem, λ, CTKAlg(), preprint="hess ", τstop=τstop)

# Reduced g_tol for optimization
opt = Optim.Options(store_trace=false, show_trace=true, extended_trace=false, g_tol=1e-3)

#================================================
Optimization run
================================================#
results = optimize(objective, gradient, hessian, λ, NewtonTrustRegion(), opt; inplace=false)

p_optimized = λ2p(SiParams, results.minimizer)
prob_optimized = SteadyStateProblem(fun, x0, p_optimized)
s_optimized = solve(prob_optimized, CTKAlg(), τstop=ustrip(u"s", 1e3u"Myr")).u

# TODO find a more generic approach to save this data... Maybe I can use datadeps?

DSi, PSi = unpack_tracers(s_optimized, grd)
DSi, = modify(DSi, PSi)
tp_opt = AIBECS.table(p_optimized)

# create directories to write in
jldsave(joinpath(output_path, "optimized_Simodel.jld2"); DSi, PSi, s_optimized, tp_opt)

println("Done!")
# print optimized parameters in md table
open(joinpath(output_path, "optimized_Si_parameters.md"), "w") do f
    pretty_table(f, AIBECS.table(p_optimized)[:,[[1,2,3,4];6:end]], nosubheader=true, tf=tf_markdown)
end

#=
# and print a LaTeX table for paper
function symbol2latex(s)
    s = latexstring(s)
    s = replace(s, "σ" => "\\sigma")
    s = replace(s, "α" => "\\alpha")
    s = replace(s, "ε" => "\\varepsilon")
    s = replace(s, "_dst" => "")
    s = replace(s, "′" => "'") # replace spaces with small spaces
    s = replace(s, r"_([A-Za-z]*)" => s"_\\mathrm{\1}")
    s = replace(s, r"₀" => "")
end
function latexify(U)
    str = string(U)
    str = replace(str, "⁻" => s"^-") # replace spaces with small spaces
    str = replace(str, "‱" => s"") # replace spaces with small spaces
    str = replace(str, "¹" => "1") # replace spaces with small spaces
    str = replace(str, "²" => "2") # replace spaces with small spaces
    str = replace(str, r"\^-(?<exp>\d+)" => s"$^{-\g<exp>}$") # add brackets around exponents
    str = replace(str, r"\s" => s"\\,") # replace spaces with small spaces
    return str
end
latexeNd(s) = replace(s, "εNd"=>"\\eNd{}")
latexbool(optimized) = optimized ? "\\checkmark" : ""
tp2 = select(tp_opt,
             :Symbol=>ByRow(symbol2latex)=>:Symbol,
             :Value,
             Symbol("Initial value"),
             :Unit=>ByRow(latexify)=>:Unit,
             :Description=>ByRow(latexeNd)=>:Description,
             :Optimizable=>ByRow(latexbool)=>:Optimized)
open(joinpath(archive_path, "optimized_parameters.tex"), "w") do f
    pretty_table(f, tp2, tf=tf_latex_simple, formatters = ft_printf("%.3g", [2, 3]), nosubheader=true)
end
=#


# Plot optimal model skill
function _locations(xmodel, obs, ux)
    M = interpolationmatrix(grd, obs)
    xmodelatobs = M * xmodel
    xdepthatobs = M * depthvec(grd)
    iobswet = findall(iswet(grd, obs))
    xwetobs = uconvert.(ux, obs.value[iobswet] * upreferred(ux))
    return xmodelatobs, xdepthatobs, iobswet, xwetobs
end
ux = mmol/m^3
DSi_model = uconvert.(ux, DSi * upreferred(ux))
xmodelatobs, xdepthatobs, iobswet, xwetobs = _locations(DSi_model, Si_obs_tuple[1], ux)

boundary = (-10,210)
x, y = ustrip.(xwetobs), ustrip.(xmodelatobs)
bw = (boundary[2]-boundary[1])/256
D = kde((x, y); boundary=(boundary, boundary), bandwidth=(bw,bw))
# calculate cumulative density from density
δx = step(D.x)
δy = step(D.y)
Q = vec(D.density) * δx * δy
idx = sortperm(Q)
Q_sorted = Q[idx]
Dcum = similar(D.density)
Dcum[idx] .= 100cumsum(Q_sorted)
ux, D, Dcum

cmap = cgrad(:oslo, categorical=true, rev=true)

contourf(D.x * ux .|> mmol/m^3, D.y * ux .|> mmol/m^3, Dcum,
         color=cmap, lc=:black, la=0, lw=0, levels=0:5:100)
plot!(collect(boundary)*ux, collect(boundary)*ux, c=:red, label="1:1")
plot!(collect(boundary)*ux, collect(boundary)*ux*0.9, c=:red, linestyle=:dot, label="±10%")
plot!(collect(boundary)*ux, collect(boundary)*ux*1.1, c=:red, linestyle=:dot, label="", 
      xlab="Observed Si(OH)₄", ylab="Modelled Si(OH)₄", legend=:bottomright)
xlims!(boundary)
ylims!(boundary)

# Add label






# Plot bSi distribution
WOA_blue   = RGB(0.231, 0.631, 0.922)
WOA_red    = RGB(0.957, 0.271, 0.188)
WOA_yellow = RGB(0.984 , 0.929, 0.298)
cmap = cgrad([WOA_blue, WOA_yellow, WOA_red])
plotverticalmean(PSi * mol/m^3 .|> μmol/m^3, grd, c=cmap, title="vertical mean opal", clim=(0,10))




# Save bSi distribution



