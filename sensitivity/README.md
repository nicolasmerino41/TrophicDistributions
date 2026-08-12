# Secondary-variable sensitivity analysis

This folder tests four modelling choices:

- grid dimension: `20`, `40`, `60`, `80`, and `100` cells per side;
- minimum viable patch: `0`, `0.5`, `1`, `1.39`, `1.5`, `2`, and `3%` of the
  landscape;
- suitability threshold: `0.15`, `0.20`, `0.25`, `0.30`, and `0.35`;
- species richness: `150`, `200`, `250`, `300`, `350`, and `400`.

Only one secondary variable changes at a time. Correlation treatments and
alternative settings use matched seeds, so differences are not driven by a new
random community for every value. When grid resolution changes, environmental
smoothing iterations scale with the square of grid dimension to preserve the
approximate spatial scale of autocorrelation relative to the landscape.
The baseline minimum-patch proportion is also preserved when grid resolution
changes, preventing grid size from silently changing the movement filter.

For each setting, the analysis retains the mean-mismatch surface across all
environment, niche-breadth, degree, and correlation cells. The final figure
shows the median change from the baseline surface, the interquartile range
(thick interval), and the 5th-95th percentile range (thin interval).

Run from the repository root:

```powershell
$env:JULIA_NUM_THREADS = 'auto'
julia --project=. sensitivity\MainSensitivity.jl
& 'C:\Program Files\R\R-4.3.2\bin\Rscript.exe' sensitivity\PlotSensitivity.R
```

`N_REPLICATES` is set in `sensitivity/Parameters.jl`. The pipeline produces:

- `Outputs/scenarios/*.tsv`: one degree-correlation surface per setting;
- `Outputs/sensitivity_summary.tsv`: compact effect summaries;
- `Outputs/sensitivity_secondary_variables.png`: the single sensitivity figure.
