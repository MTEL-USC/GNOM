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
    iroot = findfirst(current_split_dir .== "GNOM")
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

# Load some unit variables from unitful (for more readable code)
using Unitful: pM, g, mg, kg, L, mol, mmol, μmol, pmol, Mmol, yr, kyr, Myr, d, m, km, cm

# Other paths
output_path = joinpath(root_path, "output")
data_path = joinpath(root_path, "data")
# We use the head commit's 8 first characters to keep a record of which commit gave which output
using LibGit2
headcommit = string(LibGit2.GitHash(LibGit2.peel(LibGit2.GitCommit, LibGit2.head(GitRepo(root_path)))))[1:8]
archive_path = joinpath(output_path, "archive", headcommit)

#===============================================#
# Circulation grid, matrix, and OCIM2 He fluxes #
#===============================================#
# otherwise, use Circulation
const Circulation = OCIM2
# if debug=true, use an unpublished coarse OCIM (faster)
const debug = true # set to false to run with full circulation
const grd, T = let
    if debug
        circ_file = joinpath("/Users/benoitpasquier/Data/OceanGrids/OCIM0.1-lowres.jld2")
        JLD2.@load circ_file grid T
        grid, ustrip.(T)
    else
        Circulation.load()
    end
end
circname = debug ? "debug" : string(Circulation)
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

