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
include("../common_setup.jl")

# Extra packages required for Si model only
using WorldOceanAtlasTools

#==================#
# Model parameters #
#==================#
@initial_value @units @flattenable @limits @description struct SiParams{Tp} <: AbstractParameters{Tp}
    f::Tp       |   1.0 | NoUnits  | false | (0,1)  | "Fraction of non-buried bSi"
    w::Tp       | 200.0 | m/d      | true  | (0,∞)  | "Settling velocity of bSi"
    τ_up::Tp    |  30.0 | d        | true  | (0,∞)  | "Si(OH)₄ restoring timescale"
    α_up::Tp    |   1.0 | NoUnits  | true  | (0,∞)  | "Si(OH)₄ scaling for restoring"
    τ_rem::Tp   |   1.0 | d        | true  | (0,∞)  | "bSi remineralization timescale"
    DSi_geo::Tp |  80.0 | mmol/m^3 | true  | (0,∞)  | "Si(OH)₄ geological restoring target"
    τ_geo::Tp   |   1.0 | Myr      | false | (0,∞)  | "Si(OH)₄ geological restoring timescale"
    z₀::Tp      |  80.0 | m        | false | (0,∞)  | "Depth of euphotic zone"
end


#=========================================#
# Transport operator and particle sinking #
#=========================================#
T_DSi(p) = T
function T_PSi(p)
    @unpack w, f = p
    transportoperator(grd, z -> w ; z_top, z_bot, frac_seafloor=f_topo, fsedremin=f) + T
end


#=========================#
# Local sources and sinks #
#=========================#
# Uptake
relu(x) = (x ≥ 0) * x
include("obs.jl") # loads Si_obs
function uptake(DSi, p)
    @unpack α_up, τ_up, z₀ = p
    @. (z < z₀) * relu(DSi - α_up * Si_obs) / τ_up
end
# Remineralization
function remin(PSi, p)
    @unpack τ_rem = p
    @. PSi / τ_rem
end
# Slow (geological) restoring of dissolved Si
function geores(DSi, p)
    @unpack DSi_geo, τ_geo = p
    @. (DSi_geo - DSi) / τ_geo
end


#====================================#
# Right-hand side of tracer equation #
#====================================#
RHS_DSi(DSi, PSi, p) = geores(DSi, p) - uptake(DSi, p) + remin(PSi, p)
RHS_PSi(DSi, PSi, p) = uptake(DSi, p) - remin(PSi, p)


#===============#
# Problem setup #
#===============#
# Initialize parameters with initial values.
p = SiParams()
# initial guess
x0 = ustrip(upreferred(10mmol/m^3)) * ones(2nb)
# state function and its Jacobian
F = AIBECSFunction((T_DSi, T_PSi), (RHS_DSi, RHS_PSi), nb, SiParams)
# problem
prob = SteadyStateProblem(F, x0, p)






