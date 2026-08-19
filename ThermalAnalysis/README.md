# Empirical consumer-resource thermal alignment

This folder builds and analyses one auditable dataset linking documented trophic
interactions to four thermal metrics: `ctmax`, `lt50`, `ctmin`, and `ltmax`.

## Workflow

1. Run `BuildThermalDataset.R`. It harmonizes ThermoFresh, GlobTherm, and
   Comte–Olden; queries the GloBI `eats` interaction umbrella for every species
   with a usable thermal metric; incorporates the shared TetraEU metaweb; and
   writes the canonical files below.
2. Run `MainThermalAnalysis.jl`. It reads only the canonical interaction file,
   produces the four metric-specific analyses, and reports recorded-prey trait
   coverage.
3. Run `PlotThermalMetrics.R` to recreate the empirical thermal-alignment plots.
4. Run `Figure5.R` to create the manuscript figure from the updated
   consumer-level summaries. The script preserves the layout of the original
   Figure 5 and saves `Outputs/thermal_metrics/Figure5.png`.

From the repository root:

```r
source("ThermalAnalysis/BuildThermalDataset.R")
```

```julia
include("ThermalAnalysis/MainThermalAnalysis.jl")
```

```r
source("ThermalAnalysis/PlotThermalMetrics.R")
```

```r
source("ThermalAnalysis/Figure5.R")
```

`Figure5.R` requires the R packages `ggplot2`, `patchwork`, and `viridis`.

The first build can take several minutes because it queries GloBI. Each species
response is cached under `Outputs/thermal_metrics/globi_cache`; interrupted runs
can be restarted without repeating completed queries. Set the environment
variable `TD_REFRESH_GLOBI=true` only when intentionally replacing the GloBI
snapshot.

## Canonical data files

- `Data/thermal_traits_long.csv`: all usable source records, without discarding
  overlapping measurements.
- `Data/thermal_species_master.csv`: one row per harmonized species, containing
  taxonomy, habitats, source-specific summaries, selected values, interaction
  roles, and recorded degree.
- `Data/thermal_recorded_interactions.csv`: all species-level GloBI and TetraEU
  resource records for the assembled species set, including records without
  thermal information.
- `Data/thermal_interactions_master.csv`: the analysis-ready subset, with both
  partners' selected thermal values and metric-specific eligibility flags.
- `Data/thermal_consumer_coverage.csv`: recorded prey count and same-metric trait
  coverage for every consumer.
- `Data/thermal_data_manifest.csv`: source, extraction, and row-count metadata.

When the same species and metric occur in more than one trait source, every
record remains in the long table. The analysis selects a single source-level
mean using the declared priority `ThermoFresh > ComteOlden > GlobTherm`, matching
the precedence of the original analysis while now using GlobTherm correctly.

Interaction databases document reported interactions rather than complete
diets. Coverage values therefore refer to the proportion of *recorded* prey
with the corresponding thermal metric and must not be interpreted as biological
diet completeness.
