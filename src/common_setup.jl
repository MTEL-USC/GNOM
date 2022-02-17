#===================#
# Environment setup #
#===================#

# This is environment setup and packages/constants, and so on
# that are common to the Nd model and the Si model (required
# for scavenging by opal).

# Ensure that downloads are accepted through Datadeps
ENV["DATADEPS_ALWAYS_ACCEPT"] = true

# To ensure that the model always runs with the same environment,
# we activate the environment in the Project and Manifest files
# located at the root of the repository
# to avoid surprises!
root_path = let
    current_split_dir = splitpath(@__DIR__)
    iroot = findfirst((current_split_dir .== "GNOM") .| (current_split_dir .== "GNOM-main"))
    joinpath(current_split_dir[1:iroot]...)
end
# Now we activate the environment
using Pkg
Pkg.activate(root_path)
# We instantiate it
Pkg.instantiate()
# and print its status to standard output
Pkg.status()

# Packages for Si and Nd models
using AIBECS
using SuiteSparse
using SparseArrays
using LinearAlgebra
using Unitful
using Inpaintings
using Distributions
using JLD2
using F1Method
using Optim
using PrettyTables
using LaTeXStrings
using DataFrames
using DataDeps
using Downloads

# Load some unit variables from unitful (for more readable code)
using Unitful: pM, g, mg, kg, L, mol, mmol, μmol, pmol, Mmol, s, d, yr, kyr, Myr, m, km, cm
const per10000 = u"pertenthousand"
const per100 = u"percent"

# Other paths
output_path = joinpath(root_path, "output")
data_path = joinpath(root_path, "data")
# We use the head commit's 8 first characters to keep a record of which commit gave which output
using LibGit2
headcommit = if isdir(".git")
    string(LibGit2.GitHash(LibGit2.peel(LibGit2.GitCommit, LibGit2.head(GitRepo(root_path)))))[1:8]
else
    @warn """
        You probably just downloaded the GNOM without using git.
        This is OK, but be aware that the version you are using is not tracked!
        The GNOM is preferably used with git to track versions by naming runs after commit ids.
        Use `git clone https://github.com/MTEL-USC/GNOM.git` for version-controlled runs.
        Output from this run will be in folder named "nocommit".
    """
    "nocommit"
end
archive_path = joinpath(output_path, "archive", headcommit)

#===============================================#
# Circulation grid, matrix, and OCIM2 He fluxes #
#===============================================#
# otherwise, use Circulation
const Circulation = OCIM2
# if debug=true, use an unpublished coarse OCIM (faster)
# Note to external users: You will not be able to use this debug version
# because the matrix exists in my local data and is not public (and not mine to share!)
const debug = false # set to false to run with full circulation
const grd, T = let
    if debug
        circ_file = joinpath("/Users/benoitpasquier/Data/OceanGrids/OCIM0.1-lowres.jld2")
        JLD2.@load circ_file grid T
        grid, ustrip.(T)
    else
        Circulation.load()
    end
end
circname = debug ? "debug" : split(string(Circulation), '.')[2]
const nb = count(iswet(grd))
# He is tehcnically not required for the Si model
# but having it here reduces boilerplate.

const z = depthvec(grd)
const ρSW = 1.035kg/L

# Topography and constants for transport operators
const z_bot = bottomdepthvec(grd)
const z_top = topdepthvec(grd)
const f_topo = float.(isseafloorvec(grd)) # ETOPO.fractiontopo(grd)

#======================#
# Parameters interface #
#======================#
import AIBECS: @units, units
import AIBECS: @limits, limits
import AIBECS: @initial_value, initial_value
import AIBECS: @flattenable, flattenable, flatten
import AIBECS: @description, description
const ∞ = Inf
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
        else # LogitNormal with median as initial value and bounds
            m = initial_value(T, s)
            f = (m - lb) / (ub - lb)
            LocationScale(lb, ub - lb, LogitNormal(log(f/(1-f)), 1.0))
        end
    else
        nothing
    end
end
prior(::T, s::Symbol) where {T<:AbstractParameters} = prior(T,s)
prior(::Type{T}) where {T<:AbstractParameters} = Tuple(prior(T,s) for s in AIBECS.symbols(T))
prior(::T) where {T<:AbstractParameters} = prior(T)

# fallbaclk download function (for DataDeps)
function fallback_download(remotepath, localdir)
    @assert(isdir(localdir))
    filename = basename(remotepath)  # only works for URLs with filename as last part of name
    localpath = joinpath(localdir, filename)
    Downloads.download(remotepath, localpath)
    return localpath
end


