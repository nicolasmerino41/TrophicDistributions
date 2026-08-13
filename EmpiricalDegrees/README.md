# Empirical consumer degrees

This folder estimates recorded diet breadth at the same species level used by the simulation framework.

The analysis has two independent empirical samples:

- **Thermal consumers (GloBI):** the fixed list in `Data/thermal_consumers.csv` contains every consumer retained in at least one of the four thermal metrics (`ctmax`, `ctmin`, `lt50`, and `ltmax`). `MainEmpiricalDegrees.R` queries GloBI's `eats` umbrella for every listed consumer and caches the responses.
- **TetraEU consumers:** all consumers in the regional TetraEU metaweb are processed from `Data/TetraEU_pairwise_interactions.csv` file.

## Run

From the repository root:

```r
source("EmpiricalDegrees/MainEmpiricalDegrees.R")
source("EmpiricalDegrees/PlotEmpiricalDegrees.R")
```

The first run queries GloBI and may take several minutes. Completed requests are cached under `Outputs/empirical_degrees/globi_cache`, so an interrupted run can safely be restarted. Set `REFRESH_GLOBI <- TRUE` in `MainEmpiricalDegrees.R` only when intentionally replacing with current GloBI data.

## Degree definitions

All returned GloBI resource names are retained in `globi_resource_links.csv`. Two degrees are reported:

- `degree_all_taxa`: unique recorded resource names at any taxonomic resolution;
- `degree_species`: unique species-level resources supported by GloBI's resolved taxonomic path, with subspecies collapsed to their binomial species name.

The figure and combined table use `degree_species`, excluding self-links, because this most closely matches the simulation's degree definition. Self-links remain flagged in the link tables. A value of zero means that no species-level resource was returned.

Panel B reports the proportion of consumers above every degree class used in the simulations.
