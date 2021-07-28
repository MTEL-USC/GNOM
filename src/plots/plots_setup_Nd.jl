
# Load the model functions
(!isdefined(Main, :prob) || reload) && include("../Nd_model/model_setup.jl")

# Load the plots setup
(!isdefined(Main, :εNd_transects) || reload) && include("load.jl")
(!isdefined(Main, :numformat) || reload) && include("tools.jl")


