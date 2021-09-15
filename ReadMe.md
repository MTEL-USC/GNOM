# GNOM

***A Global Neodymium Ocean Model.***

This repository holds the code and data for the global steady-state model of the marine neodymium (Nd) cycle as described in *Pasquier, Hines, et al.* (in prep.)[^Pasquier_Hines_etal_GMD_2021].

This ReadMe serves as documentation for running the GNOM model in [Julia](https://julialang.org/).
Thanks to Julia's excellent built-in package manager and the [DataDeps.jl](https://github.com/oxinabox/DataDeps.jl) package, no setup should be required apart from installing Julia, downloading the GEOTRACES dataset, and activating the GNOM environment.

## Installation

1. Download and install the latest stable version of Julia from [julialang.org](https://julialang.org/) (v1.6.2 as of 1 Sep 2021).
1. Download or clone this repository on your local machine.
1. (Optional) If you are planning to optimize model parameters, you will need neodymium data from GEOTRACES, which is the only piece of data that cannot be programmatically downloaded in this project.
    The GEOTRACES Intermediate Data Product 2017 (IDP17) data must be downloaded manually from [www.geotraces.org](https://www.geotraces.org/) and, following the [GEOTRACES.jl](https://github.com/briochemc/GEOTRACES.jl) recommendation, you should preferably save the NetCDF file locally at

    ```bash
    $HOME/Data/GEOTRACES/GEOTRACES_IDP2017_v2/discrete_sample_data/netcdf/GEOTRACES_IDP2017_v2_Discrete_Sample_Data.nc
    ```

    where `$HOME` is your "home" directory.
    The remainder of the data required for running the model and its optimization (including ocean-circulation models, non-IDP17 Nd data, World Ocean Atlas data, and so on) are downloaded programmatically when the model runs.


## Single model run

We refer to GNOM v1 as the Nd-cycling (and isotope) model with "optimized" parameters as described in [*Pasquier, Hines, et al.* (2021)](), which you can modify and run in just a few seconds on your laptop.
For a single model run, all you need is to set the model up and solve the linear system.
(For instructions on how to run the optimization of the model parameters, head over to the next section.)
This is easily done in 3 steps:

1. Go to ([`cd`](https://en.wikipedia.org/wiki/Cd_%28command%29)) to the GNOM folder.
2. Start Julia
3. In Julia, type

    ```julia
    include("src/Nd_model/single_run.jl")
    ```

    to run a single simulation of the GNOM model with optimized parameters, as reported in [*Pasquier, Hines, et al.* (2021)]().

Once run, the [`src/Nd_model/single_run.jl`](src/Nd_model/single_run.jl) file should have created two variables, `DNd` and `εNd` which contain the vectors for the Nd concentration (in mol/m<sup>−3</sup>) and ε<sub>Nd</sub> (unitless). Because the `εNd` vector is the last computed variable, the output should end with something like

```julia
200160-element Vector{Float64}:
 -0.0006020528949600701
 -0.000595097968841718
  ⋮
 -0.0013770493809266426
 -0.001364787766450326
```

although the actual numerical values might be different.

These are the ε<sub>Nd</sub> values that you can convert to parts per ten thousand (‱) by typing, e.g.,

```julia
julia> εNd .|> εunit
200160-element Vector{Quantity{Float64, NoDims, Unitful.FreeUnits{(‱,), NoDims, nothing}}}:
  -6.020528949600701 ‱
 -5.9509796884171795 ‱
                     ⋮
 -13.770493809266426 ‱
  -13.64787766450326 ‱
```

where `εunit` is the ‱ unit.

Similarly, you can convert `DNd` to pM by doing:

```julia
julia> DNd * mol/m^3 .|> pM
200160-element Vector{Quantity{Float64, 𝐍 𝐋⁻³, Unitful.FreeUnits{(pM,), 𝐍 𝐋⁻³, nothing}}}:
  1.450921127576513 pM
 1.4497239997396394 pM
                     ⋮
 38.646836627631856 pM
 34.017114010235964 pM
```

where here we additionally had to apply the default unit (mol m<sup>−3</sup>) before converting.

You can edit the parameter values in [`src/Nd_model/single_run.jl`](src/Nd_model/single_run.jl) and run it again to simulate the [Nd] and ε<sub>Nd</sub> fields for different parameter values.

If you encounter any errors, please share what you did and the stacktrace in a GitHub issue and we will try to troubleshoot the issue with you as fast as possible.

## Plotting

There many packages for plotting in Julia.

### with Plots.jl

The most popular plotting package is [Plots.jl](https://github.com/JuliaPlots/Plots.jl) and AIBECS.jl provides recipes for it.
To use the Plots.jl package, simply type

```julia
julia> using Plots
```

The vectors of tracers, such as `DNd` and `εNd`, can be easily plotted with the AIBECS recipes.
For example, after running the single run above, you can use

```julia
julia> plotverticalmean(εNd .|> εunit, grd, c=:balance)
```

which will output something like

![plot_GNOM_demo](https://user-images.githubusercontent.com/4486578/133375278-623b1e7d-4b52-4806-beec-db5be82bac32.png)

For more plot types, see the [AIBECS documentation](https://juliaocean.github.io/AIBECS.jl/stable/), which contains [example tutorials](https://juliaocean.github.io/AIBECS.jl/stable/#.-Tutorials) and [how-to guides](https://juliaocean.github.io/AIBECS.jl/stable/#.-How-to-guides) with simple plots.

### With Makie.jl

The figures in  were created with [Makie.jl](https://github.com/JuliaPlots/Makie.jl) because it provides finer control to create detailed publication-quality PDFs.

AIBECS.jl does *not* provide recipes for Makie.jl, but there are a number of underlying functions to rearrange the 1D column vectors into 3D and take slices/averages over given regions/depths.
These functions are used by the plotting scripts in [`src/plots/GMDpaper/`](src/plots/GMDpaper/), which you can directly use to reproduce the plots in [*Pasquier, Hines, et al.* (2021)[^Pasquier_Hines_etal_GMD_2021]]().


## Optimization

To run the optimization, you can type

```julia
include("src/Nd_model/setup_and_optimization.jl")
```

This optimization script will randomize the initial value for the parameter values by taking a random sample from the prior distributions determined by the initial guess and the range defined in the parameter type (the `Params` struct in the `model_setup.jl` file).

## Running on a SLURM cluster

Clone your repository on your cluster, have Julia installed, and run

```bash
sbatch src/slurm/optimize_Si.sh
```

to optimize the Si model and run

```bash
sbatch src/slurm/optimize_Nd.sh
```

to optimize the Nd model. These are SLURM batch files that will request 1 node with 64GB for 20 hours to run each process.

## Si model

When running the Nd-cycling single run, the optimized Si-cycling fields required for opal scavenging are automatically downloaded from FigShare.
However, is you want, you can also edit the Si-cycle code and re-optimize it.

To run the Si-model optimization, call

```julia
include("src/Si_model/setup_and_optimization.jl")
```

Just like for the Nd-cycle model, you can modify the Si-cycle model in [`src/Si_model/model_setup.jl`](src/Si_model/model_setup.jl)

## Citation

To cite the GNOM v1 model, please cite [*Pasquier, Hines, et al.* (2021)](link_to_GMD_paper).
To cite the more general GNOM model and future versions, please cite [*Pasquier et al. (2021)*](zenodo_link?).

## Changelog

Currently a WIP, planned release v1.0 soon.



[^Pasquier_Hines_etal_GMD_2021]: Pasquier, B., Hines, S., Liang, H., Wu, Y., John, S., and Goldstein, S.: *GNOM v1.0: An optimized steady-state model of the modern marine neodymium cycle*, in preparation for submission to Geosci. Model Dev.