# TrophicDistributions 🗺️💻🐸

[![Paper](https://img.shields.io/badge/Paper-Open_Access-blue)](link_to_paper)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Official repository for the paper:
> **"Connectance and niche overlap modulate trophic effects on species distributions"**  
> *Nicolàs Merino, Núria Galiana, Miguel B. Araújo*  
> Journal, Year

## 📌 Overview

This repository contains all the Julia/R code and data needed to reproduce the analyses and figures presented in the article -. It includes the scripts to run the simulations and calculating distribution divergence metrics. It also provides the code to analyse the empirical data on consumer-resource thermal niche correlation. You can also find the MetawebCompilation.xlsx for table S1, which created Figure 2 of the main text. 

## 📓 Instructions
We want this repository to be as easy-to-use and transparent as possible. The repo is built so you can run everything straight from cloning, as long as you set the right environment (see ⚙️ Installation). Additionally, we structured the code in a modular way, so each part of the model is defined in an individual script that makes it more digestable (for instance, if you'd like to know how we built the different types of networks, you can go to Networks.jl).

If you only want to run a simulation (most likely), you only need to access two scripts. Parameters.jl gathers all possible parameters you could want to tweak (the default parameter configuration is a light version of the one used for the article). MainScript.jl is the script that executes the simulation and saves the outputs. After that, PlottingHeatmaps.R will read the outputs and produce all the heatmaps.  

NOTE: The compiled thermal data is now available in this repository for easy access but it'll be removed after the reviewing process since we don't own the data. After the reviewing period, you can email me at nicolasmerino41@gmail.com and I'll provide the .csv's. 

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
│   └── MetawebCompilation.xlsx # TABLE S1
├── Outputs/
│   ├── heatmaps/
│   ├── meanMetrics/
│   ├── tailMetrics/
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

## 🔥 Computing power
If you plan on executing extensive parameter sweeps, you will have to parallelise the simulations. The easiest way to do so is by running the following code in the console after defining your parameter configuration in Parameters.jl. 

Windows (Powershell)
```
$env:JULIA_NUM_THREADS = 7
julia --project=. SimulationsCode\MainScript.jl
```
MacOS🍎/Linux🐧
```
export JULIA_NUM_THREADS=$(7)
julia --project=. SimulationsCode/MainScript.jl
```


