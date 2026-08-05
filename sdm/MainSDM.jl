include(joinpath(@__DIR__, "..", "SimulationsCode", "Functions.jl"))
include(joinpath(@__DIR__, "Parameters.jl"))
include(joinpath(@__DIR__, "Functions.jl"))

using Dates
using .Functions.IO: save_table_tsv
using .SDMParameters
using .SDMFunctions: run_experiment, summarize_comparisons, summarize_relationship

function save_nonempty(path, rows)
    if isempty(rows)
        println("No rows to save for $(basename(path))")
        return nothing
    end
    save_table_tsv(path, rows)
    println("Saved $(length(rows)) rows to $path")
    return path
end

mkpath(OUTPUT_DIR)
println("Starting complete Figure 5 SDM experiment")
println("Communities: $(length(ENVIRONMENTS) * length(NETWORKS) * length(REGIME_INDICES) * length(CORRELATIONS) * N_REPLICATES)")
println("Julia threads: $(Threads.nthreads())")

results = run_experiment()
final_summary = summarize_comparisons(results.comparison_rows)
relationship_summary = summarize_relationship(results.comparison_rows)

save_nonempty(joinpath(OUTPUT_DIR, "community_diagnostics.tsv"), results.community_rows)
save_nonempty(joinpath(OUTPUT_DIR, "focal_metadata.tsv"), results.focal_rows)
save_nonempty(joinpath(OUTPUT_DIR, "model_results.tsv"), results.model_rows)
save_nonempty(joinpath(OUTPUT_DIR, "model_comparisons.tsv"), results.comparison_rows)
save_nonempty(joinpath(OUTPUT_DIR, "final_summary.tsv"), final_summary)
save_nonempty(joinpath(OUTPUT_DIR, "relationship_summary.tsv"), relationship_summary)
save_nonempty(joinpath(OUTPUT_DIR, "exclusions.tsv"), results.exclusions)

metadata = [(
    pipeline_version=PIPELINE_VERSION,
    completed_at=string(now()),
    environments=join(string.(ENVIRONMENTS), ","),
    networks=join(string.(NETWORKS), ","),
    regime_indices=join(REGIME_INDICES, ","),
    degrees=join(DEGREES, ","),
    correlations=join(CORRELATIONS, ","),
    replicates=N_REPLICATES,
    focal_presence_sample_size=FOCAL_PRESENCE_SAMPLE_SIZE,
    prey_presence_sample_size=PREY_PRESENCE_SAMPLE_SIZE,
    background_sample_size=BACKGROUND_SAMPLE_SIZE,
    link_completeness_levels=join(LINK_COMPLETENESS_LEVELS, ","),
    biased_background=BIASED_BACKGROUND,
    base_seed=BASE_SEED,
    julia_threads=Threads.nthreads()
)]
save_table_tsv(joinpath(OUTPUT_DIR, "run_metadata.tsv"), metadata)

println("Figure 5 SDM analysis complete.")
println("Run PlotFigure5.R to create the final figures.")
