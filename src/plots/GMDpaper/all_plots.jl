
# Set flags for re-setting up, reloading data, and reloading tools to false
reload = false
retools = false
resetup = false

include("jointPDF_Nd.jl")
include("transects_Nd.jl")
include("transects_eNd.jl")
#include("transects_extended.jl")
#include("transects_Nd_mismatches.jl")
#include("transects_eNd_mismatches.jl")
include("profiles_v2.jl")
#include("profiles_mismatch.jl")
include("horizontal_Nd_and_eNd_maps.jl")
#include("horizontal_Nd_and_eNd_mismatch_maps.jl")
include("dust_regions_eNd.jl")
include("source_maps.jl")
include("sedimentary_source_profiles.jl")
include("alpha_map.jl")
include("sink_maps.jl")

# Diagnostics
include("GMD_diagnostics_masks_map.jl")
include("GMD_diagnostics_conservative_eNd.jl")
include("GMD_diagnostics_wtag_Nd.jl")
#include("GMD_diagnostics_tables.jl")

# model-output-independent plots
include("obs_maps.jl")
#include("modelindependent/sedimentary_source_diagnostics.jl")

# parameters
include("parameters_grouped.jl")

# Parameters from our optimized runs on USC cluster (dirty)
# include("parameters_trajectories_dirty.jl")

# Si model (will only work if you have optimized it yourself)
# include("jointPDF_Si.jl")
