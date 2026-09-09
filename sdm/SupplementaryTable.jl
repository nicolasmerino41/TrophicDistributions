using Statistics
using Printf

function supplementary_rows(community_rows;
    degrees=SDMParameters.DEGREES, correlations=SDMParameters.CORRELATIONS,
    attempted=length(SDMParameters.ENVIRONMENTS) *
              length(SDMParameters.REGIME_INDICES) * SDMParameters.N_REPLICATES)
    rows = NamedTuple[]
    for degree in degrees, target_r in correlations
        group = filter(r -> r.degree == degree && r.target_r == target_r, community_rows)
        auc = filter(isfinite, Float64[r.mean_delta_auc for r in group])
        brier = filter(isfinite, Float64[r.mean_delta_brier for r in group])
        mismatch = filter(isfinite, Float64[r.mean_true_mismatch for r in group])
        push!(rows, (
            degree=degree, target_r=target_r,
            attempted_communities=attempted, n_auc=length(auc),
            n_brier=length(brier), n_consumers=sum((r.n_consumers for r in group); init=0),
            contributing_strata=length(unique([(r.environment, r.regime) for r in group])),
            mean_true_mismatch=isempty(mismatch) ? NaN : mean(mismatch),
            mean_delta_auc=isempty(auc) ? NaN : mean(auc),
            se_delta_auc=length(auc) < 2 ? NaN : std(auc)/sqrt(length(auc)),
            mean_delta_brier=isempty(brier) ? NaN : mean(brier),
            se_delta_brier=length(brier) < 2 ? NaN : std(brier)/sqrt(length(brier))
        ))
    end
    return rows
end

function write_supplementary_table(community_rows, output_dir; comparisons=nothing)
    rows = supplementary_rows(community_rows)
    Functions.IO.save_table_tsv(joinpath(output_dir, "supplementary_sdm_continuum.tsv"), rows)
    fmt(x) = isfinite(x) ? @sprintf("%.4f", x) : "NA"
    open(joinpath(output_dir, "supplementary_sdm.md"), "w") do io
        println(io, "# Supplementary table. Resource information across the interaction-relevance grid\n")
        println(io, "Finite, spatially biased presence-background SDMs across all 15 consumer degrees and 15 target niche correlations. Each cell attempts 800 communities (100 replicates in each of two environments and four niche regimes). Consumers are averaged within communities before caculating community means. Entries are mean ± standard error across communities. ΔAUC = resource-informed minus abiotic AUC; ΔBrier = abiotic minus resource-informed Brier score; positive values indicate improvement. Mismatch is the mean true trophic mismatch for eligible consumers.\n")
        if comparisons !== nothing
            flagged = count(r -> !(r.abiotic_converged && r.biotic_converged), comparisons)
            println(io, "Convergence: $flagged of $(length(comparisons)) consumer comparisons had at least one model that did not meet the convergence criterion. Their finite scores are retained, matching the existing fitting pipeline; convergence flags are available in `oracle_comparisons.tsv`.\n")
        end
        println(io, "| Degree | Target r | N (AUC/Brier) | S | Mismatch | ΔAUC ± SE | ΔBrier ± SE |")
        println(io, "|---:|---:|---:|---:|---:|---:|---:|")
        for r in rows
            println(io, "| $(r.degree) | $(@sprintf("%.6f", r.target_r)) | $(r.n_auc)/$(r.n_brier) | $(r.contributing_strata) | $(fmt(r.mean_true_mismatch)) | $(fmt(r.mean_delta_auc)) ± $(fmt(r.se_delta_auc)) | $(fmt(r.mean_delta_brier)) ± $(fmt(r.se_delta_brier)) |")
        end
    end
    return rows
end
