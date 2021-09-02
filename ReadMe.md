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

We refer to GNOM v1 as the Nd-cycling (and isotope) model with "optimized" parameters as described in [*Pasquier, Hines, et al.* (2021)](), which you can simply modify and run in just a few seconds on your laptop. 
For a single model run, all you need is to set the model up and solve the resulting equations. 
(For instructions on how to run the optimization of the model parameters, head over to the next section.)
This is easily done in 3 steps:

1. Go to ([`cd`](https://en.wikipedia.org/wiki/Cd_%28command%29)) to the GNOM folder.
2. Start Julia
3. In Julia, type

    ```julia
    include("src/Nd_model/model_setup.jl")
    ```

    to run the Julia code that will set the GNOM model up. At this stage, your model is ready to be run with the parameters in `p`.

4. Chose the parameters you want. For the optimal parameters, you can load them via XXX. Otherwise, choose parameter values, via

    ```julia
    p = Params(...)
    ```

5. Solve for the system by typing

    ```julia
    sol = solve(prob, CTKAlg(), preprint="Nd & εNd solve ", τstop=ustrip(u"s", 1e3u"Myr"))
    ```

    which will return `sol`, a large vector containing the vectors for the two isotopes nominally tracking <sup>143</sup>Nd and <sup>144</sup>Nd. To get the total Nd concentration and the ε<sub>Nd</sub> values, type

    ```julia
    DNd1, DNd2 = unpack_tracers(s_optimized, grd)
    DNd, εNd = modify(DNd1, DNd2)
    ```

    where `modify` is a function defined in the model setup to do exactly that conversion.


The optimized bSi field required for opal scavenging is automatically downloaded from XXX.

## Plotting

The vectors, matrices, and 3D arrays that you can extract from the GNOM model through the [AIBECS.jl](https://github.com/JuliaOcean/AIBECS.jl) interface can be visualized with your plotting package of choice.
Although AIBECS.jl provides recipes for [Plots.jl](https://github.com/JuliaPlots/Plots.jl), each figure in the GNOM v1.0 paper was created with [Makie.jl](https://github.com/JuliaPlots/Makie.jl) because it provides finer control on the layout.

You can reproduce the same plots as in the paper by using the code in the `src/plots/GMDpaper/`.



## Optimization

To run one optimization, you just need to call

```julia
include("src/Nd_model/run.jl")
```

## Running on a SLURM cluster

Clone your repository on your cluster, have Julia installed, and run

```bash
sbatch src/slurm/optimize_Nd.sh
```

## Si model

To run the Si-model optimization, call

```julia
include("src/Si_model/run.jl")
```

## Citation



## Changelog

Currently a WIP, planned release v1.0 soon.