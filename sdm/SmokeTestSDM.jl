include(joinpath(@__DIR__, "..", "SimulationsCode", "Functions.jl"))
include(joinpath(@__DIR__, "Parameters.jl"))
include(joinpath(@__DIR__, "Functions.jl"))

using .Functions.IO: save_table_tsv
using .SDMParameters: OUTPUT_DIR
using .SDMFunctions: run_experiment, summarize_comparisons, summarize_relationship

smoke_dir = joinpath(OUTPUT_DIR, "smoke")
checkpoint_dir = joinpath(smoke_dir, "checkpoints")
results = run_experiment(
    environments=[:autocorr],
    networks=[:random],
    regime_indices=[3],
    correlations=[0.0, 0.475, 0.95],
    replicates=1,
    checkpoint_dir=checkpoint_dir,
    resume=false
)

isempty(results.comparison_rows) && error("Smoke run produced no SDM comparisons")
summary = summarize_comparisons(results.comparison_rows)
relationship = summarize_relationship(results.comparison_rows)
mkpath(smoke_dir)
save_table_tsv(joinpath(smoke_dir, "model_comparisons.tsv"), results.comparison_rows)
save_table_tsv(joinpath(smoke_dir, "final_summary.tsv"), summary)
save_table_tsv(joinpath(smoke_dir, "relationship_summary.tsv"), relationship)

all(row -> row.abiotic_converged && row.biotic_converged,
    results.comparison_rows) || error("At least one smoke-test SDM did not converge")
println("Smoke test passed with $(length(results.comparison_rows)) model comparisons.")
