# Load model functions and all

#================================================
Loading input and output data for plots
================================================#
plot_path = joinpath(root_path, "src", "plots")
# Reloading some packages in case the optimization was already ru
using AIBECS
using JLD2
using GEOTRACES
using Statistics
using OceanBasins
using OceanographyCruises
using DataFrames
using Distributions
using DataDeps
using DataStructures
using Interpolations
using Dates
using Colors
using ColorSchemes
using ColorSchemeTools
using KernelDensity
using GeometryBasics
using Inpaintings
using Shapefile
using RemoteFiles
using Unitful
using NCDatasets
using XLSX
using Formatting
using NearestNeighbors
# Chose the plotting backend (CairoMakie for PDFs, GLMakie for "live" windows)
if use_GLMakie
    using GLMakie; GLMakie.activate!()
    GLMakie.WINDOW_CONFIG.vsync[] = false
else
    using CairoMakie; CairoMakie.activate!()
end
BACKEND = use_GLMakie ? GLMakie : CairoMakie

if isdefined(Main, :lastcommit) && lastcommit == "single"
    include("../Nd_model/obs.jl") # creates DNdobs and εNdobs
    # All the other variables should be the output of single_run.jl
else
    # Path for loading data and saving figures
    using LibGit2
    # Note that this archive path is built differently from the model setup:
    # When running the model, I use the head commit, but when plotting,
    # I use the last commit with a saved model run.
    archive_path, lastcommit = let
        allarchives_path = joinpath(output_path, "archive")
        # if ARGS is provided it should contain the commit's first 8 characters
        lastcommit = get(ARGS, 1, splitpath(first(sort(map(f -> (joinpath(allarchives_path, f), Dates.unix2datetime(mtime(f))), filter(isdir, readdir(allarchives_path, join=true))), by=last, rev=true))[1])[end])
        @show lastcommit
        archive_path = joinpath(output_path, "archive", lastcommit)
        archive_path, lastcommit
    end
    DNd, εNd, εNdobs, DNdobs, s_optimized, tp_opt = jldopen(joinpath(archive_path, "optimized_output_$(lastcommit)_run$(run_num)_$(circname).jld2")) do f
        f["DNd"], f["εNd"], f["εNdobs"], f["DNdobs"], f["s_optimized"], f["tp_opt"]
    end
    # Remake optimized parameters (note the parameters must exist in the current commit otherwise you're in trouble)
    p = Params(; zip(tp_opt.Symbol,tp_opt.Value)...)
end




iwet = findall(vec(iswet(grd)))
ρSW = 1.035u"kg/L" # approximate mean sea water density to convert mol/kg to mol/m^3
OCEANS = oceanpolygons()

#mkpath(output_path)
#mkpath(archive_path)
BG = :white
EXT = :png


fig_path = output_path
# convert model and obs to match units
uDNd = pM
DNdmodel = uconvert.(uDNd, DNd * upreferred(uDNd))
uεNd = u"pertenthousand"
εNdmodel = uconvert.(uεNd, εNd * upreferred(uεNd))
# TODO maybe use `Transects` to convert `obs` directly
# instead of reloading through GEOTRACES
# and having to reconvert the obs
using OceanographyCruises
Nd_transects = uconvert(uDNd, GEOTRACES.transects("Nd") * ρSW)
εNd_transects = uconvert(uεNd, GEOTRACES.transects("εNd"))

# Load shapefile for land
#using GeoDatasets
#landlon, landlat, landdata = GeoDatasets.landseamask(;resolution='l', grid=5)
#segments = GeoDatasets.gshhg('c', [1,6])

reload = false

