# Figure 5 SDM application

This pipeline tests whether resource information improves fitted SDMs in
the degree-correlation conditions where the simulation predicts larger trophic
mismatch.

For one eligible focal consumer per degree and simulated community, it:

1. draws up to 50 spatially biased presences from the recursive `AB` truth;
2. draws 500 biased background cells;
3. fits an abiotic logistic model with `E` and `E²`;
4. fits the same model with true local prey availability added;
5. evaluates both models against the complete `AB` map over all grid cells.

The correlation treatments use matched seeds: within each environment, niche
regime, and replicate they begin with the same degree assignment, trophic
community, environment, and niche breadth. Consumer optima and therefore focal eligibility can
change with the correlation treatment. Figure 5 highlights the low and high
correlation treatments (`0` and `0.95`) at degrees 2 and 6.

The primary result is
`AUC(resource-informed) - AUC(abiotic)`. Brier-score improvement is retained as
a secondary output.

Run from the repository root:

```powershell
$env:JULIA_NUM_THREADS = 'auto'
julia --project=. sdm\MainSDM.jl
& 'C:\Program Files\R\R-4.3.2\bin\Rscript.exe' sdm\PlotFigure5.R
```

Outputs:

- `oracle_comparisons.tsv`: one paired comparison per focal consumer;
- `oracle_summary.tsv`: summaries by environment, niche regime, degree, and
  community niche correlation;
- `figure5_selected_scenarios.tsv`: the four cases shown in Panel B;
- `figure5_sdm_application.png`: the final two-panel figure;
- `run_metadata.tsv`: run settings.

The pipeline checkpoints communities and resumes safely after interruption,
then removes those temporary checkpoints after the final tables are written.
The SDM outputs must be regenerated whenever the simulation framework changes.
