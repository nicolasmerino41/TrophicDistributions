# Empirical consumer degrees

This folder estimates recorded diet breadth at the same species level used by the simulation framework.

The analysis has two independent empirical samples:

- **GloBI consumers:** the fixed list in `Data/thermal_consumers.csv` contains the 482 verified animal consumers retained in at least one of the four thermal-alignment analyses (`ctmax`, `ctmin`, `lt50`, and `ltmax`). Every listed consumer has thermal information and at least one recorded resource with the same metric. `MainEmpiricalDegrees.R` queries GloBI's `eats` umbrella for the complete recorded resource set of each consumer, rather than restricting degree to resources with thermal data.
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

The figure and combined table use `degree_species`, excluding self-links, because this most closely matches the simulation's degree definition. Self-links remain flagged in the link tables.

Panel B reports the proportion of consumers above every degree class used in the simulations.
