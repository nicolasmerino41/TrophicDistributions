# TrophicDistributions 🗺️💻🐸

[![Paper](https://img.shields.io/badge/Paper-Open_Access-blue)](link_to_paper)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE.txt)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19053674.svg)](https://doi.org/10.5281/zenodo.19053674)

Official repository for the paper:
> **"Connectance and niche overlap determine when trophic interactions affect species distributions"**  
> *Nicolàs Merino, Núria Galiana, Miguel B. Araújo*  
> Journal, Year

## 📌 Overview

This repository contains the Julia/R code and data needed to reproduce the analyses and figures presented in the article. The revised simulation design treats focal-consumer degree as the species-level predictor and consumer-resource niche correlation as a community-level treatment. Communities contain a balanced mixture of consumer degrees, while random, modular, and cascade architectures determine how prey identities are allocated conditional on those degrees. The repository also provides the code used to analyse empirical consumer-resource thermal niche correlation.

Consumer niche optima are calibrated jointly across the food web so that both the pooled community correlation and every degree-specific correlation target the same treatment value. Outputs report separate overall and degree-level calibration checks and the maximum degree-specific error for each community.

## 📓 Instructions
We want this repository to be as easy-to-use and transparent as possible. The repo is built so you can run everything straight from cloning, as long as you set the right environment (see ⚙️ Installation). Additionally, we structured the code in a modular way, so each part of the model is defined in an individual script that makes it more digestable (for instance, if you'd like to know how we built the different types of networks, you can go to Networks.jl).

If you only want to run a simulation, start with two scripts. `Parameters.jl` contains the degree classes, community-correlation sweep, and other model settings. `MainScript.jl` executes the simulations and saves consumer-level, replicate-level, degree-summary, and community-level outputs. After that, `PlottingHeatmaps.R` reads the degree summary and produces the heatmaps.

The main degree-based outputs are written to `Outputs/degreeResults/`:

- `consumer_results.tsv`: one row per consumer and simulated community;
- `degree_replicate_summary.tsv`: one row per community and degree class;
- `degree_summary.tsv`: replicate-weighted degree-by-correlation estimates used for plotting;
- `degree_correlations.tsv`: achieved niche-correlation diagnostics by degree class.

Secondary whole-community metrics are written to `Outputs/communityMetrics/community_results.tsv`.

NOTE: The compiled thermal data is now available in this repository for easy access but it'll be removed after the reviewing process since we don't own the data. After the reviewing period, you can email me at nicolasmerino41@gmail.com and I'll provide the .csv's. They can also be found in the original sources (ThermoFresh, GlobTherm, Comte & Olden 2017, GloBI & TETRA‐EU 1.0) but that would require you some data mining you probably want to avoid :)

## 🗂️ Repository structure
```bash
TROPHICDISTRIBUTIONS/
├── Data/
│   ├── Comte_Olden_Data_Imputed.csv
│   ├── GlobalTherm_upload_02_11_17.csv
│   ├── outputs_imputed_globi_edges.csv
│   ├── TetraEU_pairwise_interactions.csv
│   ├── thermofresh_globi_metaweb_fish_predators.csv
│   ├── thermtol_comb_final.csv
│   └── SupplementaryTable1.xlsx
├── Outputs/
│   ├── heatmaps/
│   ├── degreeResults/
│   ├── communityMetrics/
│   └── thermal_metrics/
│       ├── CombinedFigures/
│       ├── ctmax/
│       ├── ctmin/
│       ├── lt50/
│       ├── ltmax/
│       └── summary_all_metrics.csv
├── SimulationsCode/
│   ├── Functions/
│   │   ├── Connectivity.jl
│   │   ├── Dynamics.jl
│   │   ├── Environment.jl
│   │   ├── Grid.jl
│   │   ├── IO.jl
│   │   ├── MechanisticCorrelation.jl
│   │   ├── Metrics.jl
│   │   ├── Networks.jl
│   │   ├── Niches.jl
│   │   ├── Simulation.jl
│   │   └── Sweep.jl
│   ├── Functions.jl
│   ├── MainScript.jl # THIS IS THE ONLY SCRIPT YOU HAVE TO RUN
│   ├── PackageLoading.jl
│   ├── Parameters.jl # TWEAK THIS SCRIPT TO CHANGE THE PARAMETRISATION
│   └── PlottingHeatmaps.R # Run this script after running the simulations
├── ThermalAnalysis/
│   ├── MainThermalAnalysis.jl
│   └── PlotThermalMetrics.R
├── .gitignore
├── Manifest.toml
├── Project.toml
└── README.md
```

## ⚙️ Installation (you must have Julia and R installed)
To set up the environment and install all dependencies:
```bash
# Clone the repository
git clone https://github.com/nicolasmerino41/TrophicDistributions.git
cd TrophicDistributions

# Start Julia with the project environment
julia --project=. -e "using Pkg; Pkg.instantiate()"

code . 
```
After running this, VSCode will open and you can run SimulationsCode/MainScript.jl. 

Before a full sweep, run the structural checks and one-community smoke simulation:

```bash
julia --project=. test/runtests.jl
```

In case needed:
 - See for Julia installation https://julialang.org/downloads/
 - See for git installation https://git-scm.com/install/
 - See for VSCode installation https://code.visualstudio.com/download

## 🔥 Computing power
If you plan on executing extensive parameter sweeps, you will have to parallelise the simulations. The easiest way to do so is by running the following code in the console after defining your parameter configuration in Parameters.jl. 

Windows (Powershell)
```
$env:JULIA_NUM_THREADS = 7 # Choose as many as cores your computer has minus 1 
julia --project=. SimulationsCode\MainScript.jl
```
MacOS🍎/Linux🐧
```
export JULIA_NUM_THREADS=$((7)) # Choose as many as cores your computer has minus 1 
julia --project=. SimulationsCode/MainScript.jl
```


