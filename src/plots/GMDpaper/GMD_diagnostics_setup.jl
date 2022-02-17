
use_GLMakie = false
include("../plots_setup_Nd.jl")


println("GMD Diagnostics computations (may take a minute)...")
# These diagnostics are shown in the GMD paper.

# List of items that are needed for diagnostics
# - H and sₖ such that H * x = ∑ₖ sₖ gives back original solution
# - sources
# - grd
# - x original (DNd concentrations and εNd here)
# - masks

# For Nd the model is affine as in the AO, so in the form A x = b
# SO I factorize A and compute all the contributions
# However, I use the notation from my papers with Mark, where the operator is called H

# Put the source vectors in a named tuple for convenience
source_functions = (s_dust, s_volc, s_sed, s_river, s_gw, s_hydro, s_tot)
source_iso_functions = (s_dust_iso, s_volc_iso, s_sed_iso, s_river_iso, s_gw_iso, s_hydro_iso, s_tot_iso)
source_names = Tuple(Symbol(string(sₖ)[3:end]) for sₖ in source_functions)
source_vectors = Tuple(sₖ(p) for sₖ in source_functions)
source_iso_vectors = Tuple(sₖ(p) for sₖ in source_iso_functions)
sources = (; zip(source_names, source_vectors)...)
sources_iso = (; zip(source_names, source_iso_vectors)...)



#===========================#
#          Masks            #
#===========================#
latN = 40
latS = -30
ATL = isatlantic(latvec(grd), lonvec(grd), OCEANS)
PAC = ispacific(latvec(grd), lonvec(grd), OCEANS)
IND = isindian(latvec(grd), lonvec(grd), OCEANS)
Indo = isindonesian(latvec(grd), lonvec(grd), OCEANS)
masks = let
    args = (latvec(grd), lonvec(grd), OCEANS)
    (N = (isatlantic(args...) .& (latN .< latvec(grd))) .| isarctic(args...) .| ismediterranean2(args...),
     S = (latvec(grd) .< latS) .| IND .| PAC .| Indo)
end



#===========================#
#    Operators H and M      #
#===========================#
# H for simple partition according to source
H = if !isdefined(Main, :neph_sink)
    sum(T_D(p).ops)
else
    sum(T_D(p).ops) + sparse(Diagonal(neph_sink(DNd, p) ./ DNd))
end
print("│ Factorizing H...")
Hf = factorize(H)
println(" Done!")
# M for partitions according to region of origin
Ms = NamedTuple(i => sparse(Diagonal(mᵢ)) for (i,mᵢ) in pairs(masks))
M = sum(Ms)
# H + M for Nd transport + relaxation in Ω
print("│ Factorizing H + M...")
HMf = factorize(H + M)
println(" Done!")
# T + M for water transport + relaxation in Ω
print("│ Factorizing T + M...")
TMf = factorize(T + M)
println(" Done!")

#=================================#
#  Partition according to source  #
#=================================#
DNdₖs = NamedTuple(k => Hf \ sₖ for (k,sₖ) in pairs(sources))
fDNdₖs = NamedTuple(k => DNdₖ ./ DNd .|> per100 for (k,DNdₖ) in pairs(DNdₖs))





#=================================#
#   Diagnosis 1: Conservative ε   #
#=================================#
ε_conservative = TMf \ (M * εNd)

#============================================#
#   Diagnosis 2: Nd water-tagged partition   #
#============================================#
# Here I infer the untagged component from
DNd_wtags = let
    tmp = NamedTuple(j => HMf \ (Mⱼ * DNd) for (j,Mⱼ) in pairs(Ms))
	(;tmp..., U = DNd .- sum(tmp))
end
fDNd_wtags = NamedTuple(j => DNdⱼ ./ DNd .|> per100 for (j,DNdⱼ) in pairs(DNd_wtags))

rediagnose = false

println("└─> Done!")

