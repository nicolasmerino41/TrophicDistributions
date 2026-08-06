include(joinpath(@__DIR__, "..", "SimulationsCode", "Functions.jl"))
include(joinpath(@__DIR__, "Parameters.jl"))
include(joinpath(@__DIR__, "Functions.jl"))

using Dates
using .Functions.IO: save_table_tsv
using .SDMParameters
using .SDMFunctions: run_experiment, summarize_results

mkpath(OUTPUT_DIR)
println("Starting the Figure 5 SDM experiment")
println("Comparison: abiotic-only versus true-resource-informed SDMs")
println("Julia threads: $(Threads.nthreads())")

results = run_experiment()
summary = summarize_results(results.rows)

save_table_tsv(joinpath(OUTPUT_DIR, "oracle_comparisons.tsv"), results.rows)
save_table_tsv(joinpath(OUTPUT_DIR, "oracle_summary.tsv"), summary)
if !isempty(results.exclusions)
    save_table_tsv(joinpath(OUTPUT_DIR, "exclusions.tsv"), results.exclusions)
end

metadata = [(
    pipeline_version=PIPELINE_VERSION,
    completed_at=string(now()),
    comparison="abiotic versus true-resource-informed",
    environments=join(string.(ENVIRONMENTS), ","),
    networks=join(string.(NETWORKS), ","),
    regime_indices=join(REGIME_INDICES, ","),
    degrees=join(DEGREES, ","),
    correlations=join(CORRELATIONS, ","),
    replicates=N_REPLICATES,
    focal_presence_sample_size=FOCAL_PRESENCE_SAMPLE_SIZE,
    background_sample_size=BACKGROUND_SAMPLE_SIZE,
    spatial_sampling_bias=true,
    base_seed=BASE_SEED,
    julia_threads=Threads.nthreads()
)]
save_table_tsv(joinpath(OUTPUT_DIR, "run_metadata.tsv"), metadata)

println("Saved oracle_comparisons.tsv and oracle_summary.tsv")
println("Run sdm/PlotFigure5.R to create the figure.")
