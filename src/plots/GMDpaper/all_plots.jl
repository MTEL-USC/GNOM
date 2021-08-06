

reload = false # setting this to true will reload model, data, plot functions and so on...

include("jointPDF_Nd.jl")
include("transects.jl")
#include("transects_extended.jl")
#include("transects_mismatch.jl")
include("profiles_v2.jl")
#include("profiles_mismatch.jl")
include("horizontal_Nd_and_eNd_maps.jl")
#include("horizontal_Nd_and_eNd_mismatch_maps.jl")
include("dust_regions_eNd.jl")
include("source_maps.jl")
include("sedimentary_source_profiles.jl")
include("alpha_map.jl")
include("sink_maps.jl")

include("parameters.jl")

# Diagnostics
include("GMD_diagnostics_masks_map.jl")
include("GMD_diagnostics_conservative_eNd.jl")
include("GMD_diagnostics_wtag_Nd.jl")
include("GMD_diagnostics_tables.jl")

# model-output-independent plots
include("obs_maps.jl")
#include("modelindependent/sedimentary_source_diagnostics.jl")
