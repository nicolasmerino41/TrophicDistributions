include(joinpath(@__DIR__, "..", "SimulationsCode", "Functions.jl"))
include(joinpath(@__DIR__, "Parameters.jl"))
include(joinpath(@__DIR__, "Functions.jl"))

using Dates
using .Functions.IO: save_table_tsv
using .SDMParameters
using .SDMFunctions: run_experiment, summarize_communities, summarize_results

mkpath(OUTPUT_DIR)
println("Comparison: abiotic-only versus resource-informed SDMs")
println("Julia threads: $(Threads.nthreads())")

results = run_experiment()
community_results = summarize_communities(results.rows)
summary = summarize_results(community_results)

save_table_tsv(joinpath(OUTPUT_DIR, "oracle_comparisons.tsv"), results.rows)
save_table_tsv(joinpath(OUTPUT_DIR, "oracle_community_results.tsv"), community_results)
save_table_tsv(joinpath(OUTPUT_DIR, "oracle_summary.tsv"), summary)
if !isempty(results.exclusions)
    save_table_tsv(joinpath(OUTPUT_DIR, "exclusions.tsv"), results.exclusions)
end

metadata = [(
    pipeline_version=PIPELINE_VERSION,
    completed_at=string(now()),
    comparison="abiotic versus true-resource-informed",
    environments=join(string.(ENVIRONMENTS), ","),
    regime_indices=join(REGIME_INDICES, ","),
    degrees=join(DEGREES, ","),
    correlations=join(CORRELATIONS, ","),
    replicates=N_REPLICATES,
    focal_design="all eligible consumers; community-level means",
    focal_presence_sample_size=FOCAL_PRESENCE_SAMPLE_SIZE,
    background_sample_size=BACKGROUND_SAMPLE_SIZE,
    spatial_sampling_bias=true,
    base_seed=BASE_SEED,
    julia_threads=Threads.nthreads()
)]
save_table_tsv(joinpath(OUTPUT_DIR, "run_metadata.tsv"), metadata)

# Checkpoints are only needed while a run is incomplete. Once the final tables
# are safely written, remove them so the published output stays compact.
function remove_checkpoints(path; attempts=5)
    !isdir(path) && return true
    GC.gc()
    for attempt in 1:attempts
        try
            rm(path; recursive=true, force=true, allow_delayed_delete=true)
            return true
        catch err
            if attempt < attempts
                sleep(0.5 * attempt)
            else
                return false
            end
        end
    end
end

remove_checkpoints(CHECKPOINT_DIR)

println("Saved consumer-level, community-level, and summary SDM results")
println("Run sdm/PlotFigure3.R to create the figure.")
