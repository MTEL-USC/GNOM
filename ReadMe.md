# GNOM

***A Global Neodymium Ocean Model.***

This repository holds the code and data for our global steady-state model of the marine Nd cycle as described in the following paper

Pasquier, B., Hines, S., Liang, H., Wu, Y., John, S., and Goldstein, S.: *GNOM v1.0: An optimized steady-state model of the modern marine neodymium cycle*, in preparation for submission to Geosci. Model Dev.

This ReadMe serves as documentation for running the GNOM model in Julia.
Thanks to Julia's excellent built-in package manager and the DataDeps.jl package, no setup should be required apart from installing Julia, downloading the GEOTRACES dataset, and activating the GNOM environment.

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

If you want to edit the parameter values, you can edit the [src/Nd_model/single_run.jl](src/Nd_model/single_run.jl), which contains a list of all the optimal parameters values and follow the same steps above.

(Note that the optimized bSi field required for opal scavenging is automatically downloaded from FigShare, but you can also edit the Si-cycle code and re-optimize it if you wish to.)

## Plotting

The vectors, matrices, and 3D arrays that you can extract from the GNOM model through the [AIBECS.jl](https://github.com/JuliaOcean/AIBECS.jl) interface can be visualized with your plotting package of choice.
Although AIBECS.jl provides recipes for [Plots.jl](https://github.com/JuliaPlots/Plots.jl), each figure in the GNOM v1.0 paper was created with [Makie.jl](https://github.com/JuliaPlots/Makie.jl) because it provides finer control on the layout.

You can reproduce the same plots as in the paper by using the code in the `src/plots/GMDpaper/`.

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

To run the Si-model optimization, call

```julia
include("src/Si_model/run.jl")
```

## Citation

To cite the GNOM v1 model, please cite [*Pasquier, Hines, et al.* (2021)](link_to_GMD_paper).
To cite the more general GNOM model and future versions, please cite [*Pasquier et al. (2021)*](zenodo_link?).

## Changelog

Currently a WIP, planned release v1.0 soon.