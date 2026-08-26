# Figure 3 SDM application

This pipeline tests whether resource information improves fitted SDMs in
the degree-correlation conditions where the simulation predicts larger trophic
mismatch.

For every eligible focal consumer of degree 2 and 6 in each simulated community, it:

1. draws up to 50 spatially biased presences from the recursive `AB` truth;
2. draws 180 (5% of the landscape, 3600 cells) biased background cells;
3. fits an abiotic logistic model with `E` and `E²`;
4. fits the same model with true local prey availability added;
5. evaluates both models against the complete `AB` map over all grid cells.

Consumers are first averaged within their community and communities are then
treated as replicates. This avoids treating consumers sharing the same simulated
world as independent observations and removes random focal-consumer selection.

The correlation treatments use matched seeds: within each environment, niche
regime, and replicate they begin with the same degree assignment, trophic
community, environment, and niche breadth. Consumer optima and therefore focal eligibility can
change with the correlation treatment. Figure 3 highlights the low and high
correlation treatments (`0` and `0.95`) at degrees 2 and 6.

The primary result is
`AUC(resource-informed) - AUC(abiotic)`. Brier-score improvement is retained as
a secondary output.

Run from the repository root:

```powershell
$env:JULIA_NUM_THREADS = 'auto'
julia --project=. sdm\MainSDM.jl
& 'C:\Program Files\R\R-4.3.2\bin\Rscript.exe' sdm\PlotFigure3.R
```

Outputs:

- `oracle_comparisons.tsv`: one comparison per eligible focal consumer;
- `oracle_community_results.tsv`: consumer means within each community and degree;
- `oracle_summary.tsv`: community-level summaries by environment, niche regime, degree, and
  community niche correlation;
- `figure3_selected_scenarios.tsv`: the four cases shown in Panel B;
- `Figure3.png`: the final two-panel figure;
- `run_metadata.tsv`: run settings.

The pipeline checkpoints communities and resumes safely after interruption.
Temporary checkpoints are removed on a best-effort basis after the final tables
are written.
The SDM outputs must be regenerated whenever the simulation framework changes.
