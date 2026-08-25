# TrophicDistributions 🗺️💻🐸

[![Paper](https://img.shields.io/badge/Paper-Open_Access-blue)](link_to_paper)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE.txt)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19053674.svg)](https://doi.org/10.5281/zenodo.19053674)

Official repository for the paper:
> **"When trophic interactions affect species distributions"**  
> *Nicolàs Merino, Núria Galiana, Miguel B. Araújo*  
> Ecology Letters, 2026

## 📌 Overview

This repository contains all the Julia/R code and data needed to reproduce the analyses and figures presented in the article *When trophic interactions affect species distributions*. It includes the scripts used to run the simulations, calculate A-AB distribution mismatch, test their application in sdm's, and evaluate the sensitivity of the results to modelling choices. It also provides the code used to analyse empirical consumer-resource thermal niche correlation.

## 📓 Instructions

We want this repository to be as easy-to-use and transparent as possible. The repo is built so you can run everything straight from cloning, as long as you set the right environment (see ⚙️ Installation). Additionally, we structured the code in a modular way, so each part of the model is defined in an individual script that makes it more digestible. For instance, the construction of niches can be found in `Niches.jl` or the construction of the environments is in `Environment.jl`.

If you only want to run the main simulations, you only need to access three scripts. `Parameters.jl` gathers the parameters you may want to tweak. `MainScript.jl` executes the simulations and saves the results. After that, `PlottingHeatmaps.R` reads those results and produces the main heatmap figure. `PlotScenarioContrasts.R` uses the same simulation outputs to quantify and plot the secondary niche-breadth and environmental-structure contrasts.

The `sdm` folder contains the application linking the simulation framework to fitted species distribution models. Run `MainSDM.jl` followed by `PlotFigure3.R` to produce Figure 3 of the article.

The `sensitivity` folder contains the sensitivity analysis for grid dimension, minimum viable patch size, suitability threshold, and species richness. Run `MainSensitivity.jl` followed by `PlotSensitivity.R` to produce the combined sensitivity figure.

The `EmpiricalDegrees` folder estimates recorded consumer degree for the species represented in the thermal analysis using GloBI, and for all consumers in the TetraEU metaweb. Run `MainEmpiricalDegrees.R` followed by `PlotEmpiricalDegrees.R`.

NOTE: The compiled thermal data is now available in this repository for easy access but it'll be removed after the reviewing process since we don't own the data. After the reviewing period, you can email me at nicolasmerino41@gmail.com and I'll provide the .csv's. They can also be found in the original sources (ThermoFresh, GlobTherm, Comte & Olden 2017, GloBI & TETRA‐EU 1.0) but that would require you some data mining you probably want to avoid :)

## 🗂️ Repository structure

```text
TROPHICDISTRIBUTIONS/
├── Data/
│   └── Input data for the empirical analyses
├── Outputs/
│   ├── heatmaps/
│   │   └── Main simulation heatmap
│   ├── scenarioContrasts/
│   │   └── Supplementary treatment contrasts
│   └── thermal_metrics/
│       └── Empirical thermal-analysis figures
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
│   ├── MainScript.jl          # THIS IS THE ONLY SCRIPT YOU HAVE TO RUN
│   ├── Parameters.jl          # TWEAK THIS SCRIPT TO CHANGE THE PARAMETRISATION
│   ├── PlottingHeatmaps.R     # Run this script after running the simulations
│   └── PlotScenarioContrasts.R # Produce the supplementary treatment contrasts
├── sdm/
│   ├── MainSDM.jl             # Run the SDM application
│   ├── Parameters.jl
│   ├── PlotFigure3.R          # Produce Figure 3
│   └── README.md
├── sensitivity/
│   ├── MainSensitivity.jl     # Run all sensitivity scenarios
│   ├── Parameters.jl
│   ├── PlotSensitivity.R      # Produce the sensitivity figure
│   └── README.md
├── ThermalAnalysis/
│   ├── MainThermalAnalysis.jl
│   └── PlotThermalMetrics.R
├── EmpiricalDegrees/
│   ├── Data/
│   │   └── thermal_consumers.csv
│   ├── Functions.R
│   ├── MainEmpiricalDegrees.R
│   ├── PlotEmpiricalDegrees.R
│   └── README.md
├── test/
│   ├── runtests.jl
│   └── sdm_runtests.jl
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

After running this, VSCode will open and you can run `SimulationsCode/MainScript.jl`.

In case needed:

- See [Julia installation](https://julialang.org/downloads/)
- See [Git installation](https://git-scm.com/install/)
- See [VSCode installation](https://code.visualstudio.com/download)

## 🔥 Computing power

If you plan on executing extensive parameter sweeps, you will have to parallelise the simulations. The easiest way to do so is to run the following code in the console after defining your parameter configuration in `Parameters.jl`.

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
