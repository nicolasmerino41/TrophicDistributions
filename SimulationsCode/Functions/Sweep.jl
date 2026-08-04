module Sweep

using Random
using Statistics
using ..Parameters: CORR_RANGE, N_CORR, NREP, BASE_SEED, NX, NY, NCELLS,
                    S, BASAL_FRAC, DEGREE_CLASSES, Emin_patch,
                    USE_CONNECTIVITY_FILTER, OUTDIR, TAIL_THRESH
using ..Connectivity: make_workspaces
using ..Niches: regimes
using ..Simulation: simulate_one!
using ..Metrics: mismatch_q90

export sweep_all, ENVKINDS, NETFAMS, NETNAMES, regime_name

const ENVKINDS = [:random, :autocorr]
const NETFAMS = [:random, :modular, :cascade]
const NETNAMES = Dict(
    :random => "Random",
    :modular => "Modular",
    :cascade => "Cascade"
)

regime_name(reg) = reg.name

function degree_replicate_rows(simulation, metadata)
    rows = NamedTuple[]
    for degree in DEGREE_CLASSES
        assigned = filter(row -> row.target_degree == degree, simulation.consumers)
        eligible = filter(row -> row.eligible && isfinite(row.mismatch), assigned)
        mismatch = [row.mismatch for row in eligible]
        push!(rows, merge(metadata, (
            community_id=simulation.community.community_id,
            degree=degree,
            n_assigned=length(assigned),
            n_eligible=length(eligible),
            mean_mismatch=isempty(mismatch) ? NaN : mean(mismatch),
            median_mismatch=isempty(mismatch) ? NaN : median(mismatch),
            mismatch_q90=mismatch_q90(mismatch),
            mismatch_frac_gt=isempty(mismatch) ? NaN : mean(mismatch .> TAIL_THRESH),
            degree_correlation=simulation.degree_correlations[degree]
        )))
    end
    return rows
end

function aggregate_degree_rows(replicate_rows)
    grouped = Dict{Tuple{Symbol,Symbol,Int,Float64,Int},Vector{NamedTuple}}()
    for row in replicate_rows
        key = (row.environment, row.network, row.regime_id, row.target_r, row.degree)
        push!(get!(grouped, key, NamedTuple[]), row)
    end

    summary = NamedTuple[]
    for (key, rows) in grouped
        environment, network, regime_id, target_r, degree = key
        means = filter(isfinite, [row.mean_mismatch for row in rows])
        medians = filter(isfinite, [row.median_mismatch for row in rows])
        q90s = filter(isfinite, [row.mismatch_q90 for row in rows])
        fractions = filter(isfinite, [row.mismatch_frac_gt for row in rows])
        degree_corrs = filter(isfinite, [row.degree_correlation for row in rows])
        push!(summary, (
            environment=environment,
            network=network,
            regime_id=regime_id,
            regime=rows[1].regime,
            target_r=target_r,
            degree=degree,
            n_communities=length(rows),
            n_communities_with_data=length(means),
            n_communities_correlation_ok=count(row -> row.correlation_ok, rows),
            n_consumers_assigned=sum(row.n_assigned for row in rows),
            n_consumers_eligible=sum(row.n_eligible for row in rows),
            mean_mismatch=isempty(means) ? NaN : mean(means),
            sd_replicate_mean=length(means) < 2 ? NaN : std(means),
            mean_median_mismatch=isempty(medians) ? NaN : mean(medians),
            mean_q90_mismatch=isempty(q90s) ? NaN : mean(q90s),
            mean_frac_gt=isempty(fractions) ? NaN : mean(fractions),
            mean_degree_correlation=isempty(degree_corrs) ? NaN : mean(degree_corrs),
            mean_max_degree_correlation_error=mean(
                row.max_degree_correlation_error for row in rows
            )
        ))
    end
    sort!(summary, by=row -> (
        findfirst(==(row.environment), ENVKINDS),
        findfirst(==(row.network), NETFAMS),
        row.regime_id, row.target_r, row.degree
    ))
    return summary
end

function sweep_all()
    Rvals = collect(range(CORR_RANGE[1], CORR_RANGE[2], length=N_CORR))
    workspaces = make_workspaces()

    combinations = NamedTuple[]
    community_id = 0
    for (env_idx, env) in enumerate(ENVKINDS),
        (net_idx, net) in enumerate(NETFAMS),
        (regime_id, regime) in enumerate(regimes),
        (corr_idx, target_r) in enumerate(Rvals),
        replicate in 1:NREP
        community_id += 1
        push!(combinations, (
            community_id=community_id,
            env_idx=env_idx,
            environment=env,
            net_idx=net_idx,
            network=net,
            regime_id=regime_id,
            regime=regime,
            corr_idx=corr_idx,
            target_r=target_r,
            replicate=replicate
        ))
    end

    println("Threads: ", Threads.nthreads())
    println("Grid: $(NX)x$(NY) cells=$(NCELLS), S=$(S), basal_frac=$(BASAL_FRAC), Emin=$(Emin_patch), connectivity_filter=$(USE_CONNECTIVITY_FILTER)")
    println("OUTDIR: ", OUTDIR)
    println("Degree classes: $(join(DEGREE_CLASSES, ", "))")
    println("Sweep: env=$(length(ENVKINDS)) × net=$(length(NETFAMS)) × regimes=$(length(regimes)) × r=$(N_CORR) × reps=$(NREP)")

    simulations = Vector{Any}(undef, length(combinations))
    Threads.@threads for idx in eachindex(combinations)
        combination = combinations[idx]
        thread_id = Threads.threadid()
        rng = MersenneTwister(
            BASE_SEED +
            10_000_000 * combination.env_idx +
            1_000_000 * combination.net_idx +
            100_000 * combination.regime_id +
            10_000 * combination.corr_idx +
            combination.replicate
        )
        simulations[idx] = simulate_one!(
            rng,
            workspaces[thread_id],
            combination.environment,
            combination.network,
            combination.regime,
            combination.target_r,
            combination.community_id
        )
    end

    consumer_results = NamedTuple[]
    community_results = NamedTuple[]
    degree_replicate = NamedTuple[]
    degree_correlations = NamedTuple[]

    for (combination, simulation) in zip(combinations, simulations)
        metadata = (
            replicate=combination.replicate,
            environment=combination.environment,
            network=combination.network,
            regime_id=combination.regime_id,
            regime=combination.regime.name,
            target_r=combination.target_r,
            achieved_r=simulation.community.achieved_r,
            correlation_ok=simulation.community.correlation_ok,
            overall_correlation_ok=simulation.community.overall_correlation_ok,
            degree_correlation_ok=simulation.community.degree_correlation_ok,
            max_degree_correlation_error=simulation.community.max_degree_correlation_error
        )
        append!(consumer_results, [merge(metadata, row) for row in simulation.consumers])
        push!(community_results, merge(metadata, simulation.community))
        append!(degree_replicate, degree_replicate_rows(simulation, metadata))
        for degree in DEGREE_CLASSES
            push!(degree_correlations, merge(metadata, (
                community_id=combination.community_id,
                degree=degree,
                degree_correlation=simulation.degree_correlations[degree]
            )))
        end
    end

    degree_summary = aggregate_degree_rows(degree_replicate)
    return (
        consumer_results=consumer_results,
        degree_replicate=degree_replicate,
        degree_summary=degree_summary,
        community_results=community_results,
        degree_correlations=degree_correlations,
        degree_values=copy(DEGREE_CLASSES),
        correlation_values=Rvals
    )
end

end
