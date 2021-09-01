# GNOM

***A Global Neodymium Ocean Model.***

This repository holds the code and data accompanying our paper

Pasquier, B., Hines, S., Liang, H., Wu, Y., John, S., and Goldstein, S.: *GNOM v1.0: An optimized steady-state model of the modern marine neodymium cycle*, in preparation for submission to Geosci. Model Dev.


## Documentation

This ReadMe serves as documentation for running the GNOM model in Julia.
Thanks to Julia's excellent built-in package manager and the DataDeps.jl package, no setup should be required apart from installing Julia, downloading the GEOTRACES dataset, and activating the GNOM environment.

### Installation

1. Download and install the latest stable version of Julia from [julialang.org](https://julialang.org/) (v1.6.2 as of 1 Sep 2021).
2. If you are planning to run some optimization, you will need neodymium data from the GEOTRACES Intermediate Data Product 2017 (IDP17), which is the only piece of data that cannot be programmatically downloaded in this project.
    IDP17 data must be downloaded manually (from [www.geotraces.org](https://www.geotraces.org/)) and, following the [GEOTRACES.jl](https://github.com/briochemc/GEOTRACES.jl) recommendation, that the NetCDF file be located at

    ```bash
    $HOME/Data/GEOTRACES/GEOTRACES_IDP2017_v2/discrete_sample_data/netcdf/GEOTRACES_IDP2017_v2_Discrete_Sample_Data.nc
    ```

    The remainder of the data (Ocean circulation models, pre- and post- GEOTRACES IDP17 data, World Ocean Atlas data, and so on) are downloaded programmatically when the model runs.

1. Download or clone this repository on your local machine.

### Single model run

We refer to GNOM v1.0 as the Nd-cycling (and isotope) model with "optimized" parameters. If you want to run the optimization yourself, head over to the next section.

For a single model run, all you need is to set the model up and solve the resulting equations. This is easily done in 3 steps:

1. Go to ([`cd`](https://en.wikipedia.org/wiki/Cd_%28command%29)) to the GNOM folder.
2. Start Julia
3. Type

    ```julia
    include("src/Nd_model/model_setup.jl")
    ```

    to set the model up. At this stage, your model is ready to be run with the parameters in `p`.

4. Chose the parameters you want. For the optimal parameters, you can load them via XXX. Otherwise, choose parameter values, via

    ```julia
    p = Params()
    ```

5. Solve for the system

The optimized bSi field required for opal scavenging is automatically downloaded from XXX.

### Plotting

The vectors, matrices, and 3D arrays that you can extract from the GNOM model through the [AIBECS.jl](https://github.com/JuliaOcean/AIBECS.jl) interface can be visualized with your plotting package of choice.
Although AIBECS.jl provides recipes for [Plots.jl](https://github.com/JuliaPlots/Plots.jl), each figure in the GNOM v1.0 paper was created with [Makie.jl](https://github.com/JuliaPlots/Makie.jl) because it provides finer control on the layout.

You can reproduce the same plots as in the paper by using the code in the `src/plots/GMDpaper/`.



### Optimization

To run one optimization, you just need to call

```julia
include("src/Nd_model/run.jl")
```

### Running on a SLURM cluster

Clone your repository on your cluster, have Julia installed, and run

```bash
sbatch src/slurm/optimize_Nd.sh
```

### Si model

To run the Si-model optimization, call

```julia
include("src/Si_model/run.jl")
```

## Changelog

Release v1.0