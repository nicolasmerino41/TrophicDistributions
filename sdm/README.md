# Figure 5: SDMs under imperfect information

This folder contains the complete SDM application of the virtual trophic worlds.
It uses Julia for simulation and model fitting and dependency-free R scripts for
the figures.

## Experimental design

The experiment mirrors the main simulation grid:

- 15 focal-consumer degrees;
- 15 community niche-correlation treatments from 0 to 0.95;
- random and autocorrelated environments;
- random, modular, and cascade food webs;
- all four niche-breadth regimes;
- the replicate count defined by `NREP` in `SimulationsCode/Parameters.jl`.

For each community, one eligible focal consumer is selected from every degree
class. Each focal species is fitted with the same spatially biased presence
sample and target-group background under six models:

1. abiotic-only baseline;
2. true prey availability and complete trophic information (`oracle`);
3. prey availability estimated from finite, biased prey occurrences with 100%
   of trophic links known;
4. the same estimated prey distributions with 75% of links known;
5. the same with 50% of links known;
6. the same with 25% of links known.

Unknown-link treatments use nested subsets of each consumer's true prey list.
Every prey SDM is fitted only once per community and reused by all focal
consumers. Predicted resource availability is the probability that at least one
known prey is present. All focal models are evaluated against the complete true
interaction-constrained (`AB`) distribution.

## Data handoff

`SimulationsCode/Functions/Simulation.jl` exposes `simulate_world!`, which
returns the environment, food web, degrees, niches, abiotic distributions (`A`),
interaction-constrained distributions (`AB`), and correlation diagnostics in
memory. The SDM workflow calls this function directly. It does not require the
main sweep to write thousands of full distribution rasters to disk.

## Run the analysis

From the repository root, first run the reduced end-to-end check:

```powershell
& 'C:\Users\nicol\.julia\juliaup\julia-1.12.5+0.x64.w64.mingw32\bin\julia.exe' --project=. sdm\SmokeTestSDM.jl
```

Then run the complete analysis. Using multiple Julia threads is strongly
recommended:

```powershell
$env:JULIA_NUM_THREADS = 'auto'
& 'C:\Users\nicol\.julia\juliaup\julia-1.12.5+0.x64.w64.mingw32\bin\julia.exe' --project=. sdm\MainSDM.jl
```

Every completed community is checkpointed in `sdm/Outputs/checkpoints`. Restart
the same command after an interruption and completed communities are loaded
instead of rerun. Delete that versioned checkpoint directory only when changing
the experimental design or model definitions.

## Outputs

- `community_diagnostics.tsv`: simulation and prey-model diagnostics;
- `focal_metadata.tsv`: focal identity, true mismatch, degree, and sample sizes;
- `model_results.tsv`: absolute metrics for every fitted model;
- `model_comparisons.tsv`: paired biotic-minus-abiotic comparisons;
- `final_summary.tsv`: replicate-level degree-by-correlation summaries;
- `relationship_summary.tsv`: binned relationship between true mismatch and SDM gain;
- `exclusions.tsv`: communities or degree classes that could not be analysed;
- `run_metadata.tsv`: complete parameter record for the run.

Positive `delta_auc`, `delta_brier`, `delta_logloss`, and `delta_jaccard` always
mean that the biotic model outperformed the abiotic-only baseline.

## Create figures

```powershell
& 'C:\Program Files\R\R-4.3.2\bin\Rscript.exe' sdm\PlotFigure5.R
& 'C:\Program Files\R\R-4.3.2\bin\Rscript.exe' sdm\PlotAllSurfaces.R
```

`PlotFigure5.R` combines the perfect-information mismatch surface, four SDM
improvement surfaces, and the mismatch-performance relationship. The reference
environment, architecture, and breadth regime are set at the top of that file.
`PlotAllSurfaces.R` creates the complete set of treatment-specific supplementary
heatmaps.
