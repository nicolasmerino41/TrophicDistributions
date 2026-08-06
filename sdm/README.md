# Figure 5: translating the framework into SDM performance

This folder contains one controlled application of the theoretical framework.
It asks a single question:

> Does adding resource distributions improve an SDM most strongly where the
> framework predicts that trophic interactions alter consumer distributions?

## Method

For every simulated community and focal-consumer degree:

1. The existing simulation generates the environmental raster, food web,
   abiotic distribution (`A`), and interaction-constrained distribution (`AB`).
2. A finite, spatially biased sample of focal-consumer presences is drawn from
   the true `AB` distribution.
3. An abiotic presence-background SDM is fitted using the environment and its
   quadratic term.
4. A resource-informed SDM adds the true spatial availability of at least one
   prey species.
5. Both models are evaluated against the complete simulated `AB` distribution.

The main performance measure is the increase in conventional whole-landscape
AUC:

`AUC(resource-informed model) - AUC(abiotic model)`

Positive values mean that adding resource distributions improved discrimination.
Brier prediction error is retained in the result table as a secondary metric.

The abiotic model is expected to perform strongly because the simulated niche
is generated from one environmental variable and the fitted SDM uses the
matching linear and quadratic terms. Figure 5 therefore tests for a modest but
systematic improvement from biotic information rather than expecting a poorly
performing abiotic baseline.

The experiment uses imperfect focal-occurrence data but complete resource-range
information. It is therefore a controlled demonstration of applicability, not
an analysis of uncertainty in empirical prey distributions.

## Files

- `Parameters.jl`: SDM sampling, model, grid, and highlighted-scenario settings.
- `Functions.jl`: sampling, logistic SDMs, evaluation, sweep, and summaries.
- `MainSDM.jl`: ready-to-run analysis.
- `PlotFigure5.R`: self-contained two-panel Figure 5 script.

## Run

From the repository root:

```powershell
$env:JULIA_NUM_THREADS = 'auto'
& 'C:\Users\nicol\.julia\juliaup\julia-1.12.5+0.x64.w64.mingw32\bin\julia.exe' --project=. sdm\MainSDM.jl
& 'C:\Program Files\R\R-4.3.2\bin\Rscript.exe' sdm\PlotFigure5.R
```

The Julia run checkpoints every community and safely resumes after interruption.
The valid resource-informed results from the previously completed run have
already been preserved in the current output files, so rerunning is unnecessary
unless parameters change.

## Outputs

- `oracle_comparisons.tsv`: one paired SDM comparison per focal consumer.
- `oracle_summary.tsv`: replicate summary for every treatment, degree, and correlation.
- `figure5_selected_scenarios.tsv`: the four values plotted in Panel B.
- `figure5_sdm_application.png` and `.pdf`: final two-panel figure.
- `run_metadata.tsv`: settings used by a newly completed Julia run.

Panel A averages theoretical mismatch across both environments, all three
food-web architectures, and all four niche regimes. Panel B shows the mean and
95% interval for degree 2 and degree 6 at correlations approximately 0.1 and
0.8. The exact simulated correlation treatments are 0.0679 and 0.8143. Colors
correspond directly to the markers in Panel A.
