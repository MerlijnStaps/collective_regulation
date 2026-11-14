# Collective regulation

Authors: Merlijn Staps (Washington University in St. Louis), Deborah M. Gordon (Stanford University), Corina E. Tarnita (Princeton University)

Corresponding author: Merlijn Staps (merlijnstaps@gmail.com).

This repository includes the code to perform analyses and recreate the figures. 

## Activity distributions (Fig. 2abc)

To perform the analysis, run:

```console
$ julia activity_distributions.jl
```

This will create CSV files ``motif_...`` in the ``results`` folder. To create the figures, run:

```console
$ Rscript R/activity_distributions.R
```

This will create plots ``activity_distributions.pdf`` in the ``plots`` folder. 

## Optimal strategies (Fig. 2def)


To perform the analysis, run:

```console
$ julia optimal_strategies.jl
```

The results are stored in ``optimal_strategies.csv`` in the ``results`` folder. To create the figures, run:

```console
$ Rscript R/optimal_strategies.R
```

This will create plots ``optimal_strategies_motif_I.pdf`` and  ``optimal_strategies_motif_II.pdf`` in the ``plots`` folder. 

## Simulations

The simulation code was produced on Julia v1.11.1. The dependencies can be found in ``Project.toml`` and installed using 

```console
$ julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The code for the simulations can be found in `flexibility.jl`. 

The script `setup.jl` creates CSV files with parameter settings for the simulations, which are stored in the folder `sims`. 

The script `run.jl` can be used to run simulations as follows:

```console
$ julia --project=. run.jl evo1 100 123
```

This will run simulations with parameters from lines 100 through 123 from `sims/evo1.csv` (be sure to unzip the file first), and store the outputs in the folder `results/evo1`. It will also output the winning strategy found in the simulation and its fitness.

## Figures

The outputs of the simulations can be found in the folder `results`. 

The code to produce the figures was written in R (v4.3.2). The R scripts to generate the figures can be found in the folder `R`. The only dependency is RColorBrewer (version 1.1.3).
