# BioDist 🗺️🕸️🐸

[![Paper](https://img.shields.io/badge/Paper-Open_Access-blue)](link_to_paper)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Official code for the paper:
> **"Connectance and niche overlap modulate trophic effects on species distributions"**  
> *Nicolàs Merino, Núria Galiana, Miguel B. Araújo*  
> Journal, Year

## 📌 Overview

This repository contains all the Julia/R code and data needed to reproduce the analyses and figures presented in the article -. It includes the scripts to run the simulations and calculating distribution divergence metrics. It also provides the necessary code to analyse the empirical data on consumer-resource thermal niche correlation.

## 🗂️ Repository Structure
```bash
├── README.md               
├── Code/
  |── Functions.jl
  |── MainScript.jl # Only this script needs to be run
  |── PackageLoading.jl
  |── Plotting.jl                
├── Figures/
  |── Correlation_results_for_scenarios_ER_PL_MOD.png   # Figure 2 of the paper
  |── error_vs_structure.png                            # Figure 3 of the paper             
├── Outputs/                # .jls objects to be saved
├── paper.pdf               # PDF of the paper
├── Project.toml            # Package dependencies
├── Manifest.toml           # Pinned package versions for exact reproducibility. 
├── LICENSE                 # License information
```

## ⚙️ Installation
To set up the environment and install all dependencies:
```bash
# Clone the repository
git clone https://github.com/your-username/FromStructureToDynamics.git
cd FromStructureToDynamics

# Start Julia with the project environment
julia --project=.

# Inside Julia:
using Pkg
Pkg.instantiate()

# Finally, run:
Code/MainScript.jl
```
## 📊 Figures




