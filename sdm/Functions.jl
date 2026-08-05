module SDMFunctions

using LinearAlgebra
using Random
using Serialization
using Statistics
using ..Functions.Parameters: NX, NY, NCELLS
using ..Functions.Grid: x_of, y_of
using ..Functions.Connectivity: make_workspaces
using ..Functions.Niches: regimes
using ..Functions.Simulation: simulate_world!
using ..SDMParameters

export LogisticModel, fit_ridge_logistic, predict_probability,
       sampling_bias_surface, sample_presence_background,
       oracle_resource_availability, estimated_resource_availability,
       choose_focal_consumers, run_world_sdms, build_jobs,
       run_experiment, summarize_comparisons, summarize_relationship

struct LogisticModel
    coefficients::Vector{Float64}
    converged::Bool
    iterations::Int
end

@inline function logistic(value::Float64)
    if value >= 0
        return 1.0 / (1.0 + exp(-value))
    end
    exp_value = exp(value)
    return exp_value / (1.0 + exp_value)
end

predict_probability(model::LogisticModel, design::Matrix{Float64}) =
    logistic.(design * model.coefficients)

function fit_ridge_logistic(
    design::Matrix{Float64},
    response::Vector{Float64};
    ridge::Float64=LOGISTIC_RIDGE,
    max_iters::Int=LOGISTIC_MAX_ITERS,
    tolerance::Float64=LOGISTIC_TOL
)
    nrow, ncol = size(design)
    length(response) == nrow || error("Response and design have different lengths")
    all(y -> y == 0.0 || y == 1.0, response) || error("Response must be binary")

    coefficients = zeros(Float64, ncol)
    penalty = Diagonal(vcat(0.0, fill(ridge, ncol - 1)))
    for iteration in 1:max_iters
        linear_predictor = design * coefficients
        probability = logistic.(linear_predictor)
        weight = clamp.(probability .* (1.0 .- probability), 1e-7, Inf)
        working_response = linear_predictor .+ (response .- probability) ./ weight
        weighted_design = design .* reshape(weight, :, 1)
        updated = (design' * weighted_design + penalty) \
                  (design' * (weight .* working_response))
        if maximum(abs.(updated .- coefficients)) <= tolerance
            return LogisticModel(updated, true, iteration)
        end
        coefficients = updated
    end
    return LogisticModel(coefficients, false, max_iters)
end

function standardized_environment(environment::Vector{Float64})
    scale = std(environment)
    scale > 0 || error("Environmental raster has zero variance")
    return (environment .- mean(environment)) ./ scale
end

function design_matrix(
    environment_z::Vector{Float64};
    biotic_predictor::Union{Nothing,Vector{Float64}}=nothing
)
    base = hcat(ones(length(environment_z)), environment_z, environment_z.^2)
    biotic_predictor === nothing && return base
    length(biotic_predictor) == length(environment_z) ||
        error("Biotic predictor has the wrong number of cells")
    return hcat(base, biotic_predictor)
end

function sampling_bias_surface(
    strength::Float64=SAMPLING_BIAS_STRENGTH,
    floor_weight::Float64=SAMPLING_BIAS_FLOOR
)
    center_x = 0.25 * NX
    center_y = 0.75 * NY
    weights = Vector{Float64}(undef, NCELLS)
    for cell in 1:NCELLS
        dx = (x_of(cell) - center_x) / NX
        dy = (y_of(cell) - center_y) / NY
        weights[cell] = floor_weight + exp(-strength * (dx^2 + dy^2))
    end
    return weights
end

function weighted_sample_without_replacement(
    rng::AbstractRNG,
    candidates::Vector{Int},
    weights::Vector{Float64},
    requested::Int
)
    requested > 0 || return Int[]
    isempty(candidates) && return Int[]
    sample_size = min(requested, length(candidates))
    pool = copy(candidates)
    pool_weights = Float64[max(weights[cell], 0.0) for cell in pool]
    selected = Int[]
    sizehint!(selected, sample_size)
    for _ in 1:sample_size
        total_weight = sum(pool_weights)
        chosen = if total_weight <= 0
            rand(rng, eachindex(pool))
        else
            draw = rand(rng) * total_weight
            cumulative = 0.0
            found = lastindex(pool)
            for index in eachindex(pool)
                cumulative += pool_weights[index]
                if draw <= cumulative
                    found = index
                    break
                end
            end
            found
        end
        push!(selected, pool[chosen])
        deleteat!(pool, chosen)
        deleteat!(pool_weights, chosen)
    end
    return selected
end

function sample_presence_background(
    rng::AbstractRNG,
    truth::BitVector,
    accessibility::Vector{Float64};
    n_presence::Int,
    n_background::Int=BACKGROUND_SAMPLE_SIZE,
    biased_background::Bool=BIASED_BACKGROUND
)
    presences = weighted_sample_without_replacement(
        rng, findall(identity, truth), accessibility, n_presence
    )
    observed = falses(length(truth))
    observed[presences] .= true
    background_candidates = findall(!, observed)
    background_weights = biased_background ? accessibility : ones(Float64, length(truth))
    background = weighted_sample_without_replacement(
        rng, background_candidates, background_weights, n_background
    )
    return presences, background
end

function oracle_resource_availability(world, focal_consumer::Int; resources=nothing)
    selected = resources === nothing ? world.prey[focal_consumer] : resources
    availability = BitVector(falses(NCELLS))
    for resource in selected
        availability .|= world.AB[resource]
    end
    return Float64.(availability)
end

function estimated_resource_availability(
    prey_predictions::Dict{Int,Vector{Float64}},
    resources::Vector{Int}
)
    availability = zeros(Float64, NCELLS)
    valid_resources = [resource for resource in resources if haskey(prey_predictions, resource)]
    isempty(valid_resources) && return availability
    log_absence = zeros(Float64, NCELLS)
    for resource in valid_resources
        probability = clamp.(prey_predictions[resource], 0.0, 1.0 - 1e-10)
        log_absence .+= log1p.(-probability)
    end
    return 1.0 .- exp.(log_absence)
end

function jaccard_similarity(predicted::BitVector, truth::BitVector)
    union_size = count(predicted .| truth)
    union_size == 0 && return 1.0
    return count(predicted .& truth) / union_size
end

function auc_score(probability::Vector{Float64}, truth::BitVector)
    positives = count(truth)
    negatives = length(truth) - positives
    (positives == 0 || negatives == 0) && return NaN
    order = sortperm(probability)
    ranks = zeros(Float64, length(probability))
    position = 1
    while position <= length(order)
        final_position = position
        value = probability[order[position]]
        while final_position < length(order) &&
              probability[order[final_position + 1]] == value
            final_position += 1
        end
        average_rank = 0.5 * (position + final_position)
        for rank_position in position:final_position
            ranks[order[rank_position]] = average_rank
        end
        position = final_position + 1
    end
    positive_rank_sum = sum(ranks[truth])
    return (positive_rank_sum - positives * (positives + 1) / 2) /
           (positives * negatives)
end

function best_tss_threshold(probability::Vector{Float64}, response::Vector{Float64})
    positive = response .== 1.0
    negative = .!positive
    n_positive = count(positive)
    n_negative = count(negative)
    (n_positive == 0 || n_negative == 0) && return 0.5
    best_threshold = 0.5
    best_tss = -Inf
    for threshold in range(0.0, 1.0, length=201)
        predicted = probability .>= threshold
        sensitivity = count(predicted .& positive) / n_positive
        specificity = count((.!predicted) .& negative) / n_negative
        tss = sensitivity + specificity - 1.0
        if tss > best_tss
            best_tss = tss
            best_threshold = threshold
        end
    end
    return best_threshold
end

function binary_log_loss(probability::Vector{Float64}, truth::BitVector)
    p = clamp.(probability, 1e-8, 1.0 - 1e-8)
    y = Float64.(truth)
    return -mean(y .* log.(p) .+ (1.0 .- y) .* log.(1.0 .- p))
end

function focal_truth_mismatch(world, focal_consumer::Int)
    abiotic = world.A[focal_consumer]
    realized = world.AB[focal_consumer]
    union_size = count(abiotic .| realized)
    union_size == 0 && return NaN
    return 1.0 - count(abiotic .& realized) / union_size
end

function choose_focal_consumers(
    rng::AbstractRNG,
    world,
    degree::Int;
    number::Int=FOCALS_PER_DEGREE,
    minimum_presence_cells::Int=MIN_FOCAL_PRESENCE_CELLS
)
    candidates = [
        consumer for consumer in world.consumers
        if world.realized_degree[consumer] == degree &&
           count(world.A[consumer]) > 0 &&
           count(world.AB[consumer]) >= minimum_presence_cells
    ]
    shuffle!(rng, candidates)
    return candidates[1:min(number, length(candidates))]
end

function fit_sampled_distribution(
    rng::AbstractRNG,
    truth::BitVector,
    abiotic_design::Matrix{Float64},
    accessibility::Vector{Float64};
    n_presence::Int,
    minimum_presence_cells::Int
)
    true_cells = count(truth)
    if true_cells < minimum_presence_cells
        return (success=false, reason="insufficient_true_cells", probability=nothing,
                converged=false, sampled_presences=0, sampled_background=0)
    end
    presences, background = sample_presence_background(
        rng, truth, accessibility;
        n_presence=n_presence,
        n_background=BACKGROUND_SAMPLE_SIZE
    )
    if isempty(presences) || isempty(background)
        return (success=false, reason="sampling_failed", probability=nothing,
                converged=false, sampled_presences=length(presences),
                sampled_background=length(background))
    end
    training_cells = vcat(presences, background)
    response = vcat(ones(Float64, length(presences)), zeros(Float64, length(background)))
    model = fit_ridge_logistic(abiotic_design[training_cells, :], response)
    return (success=true, reason="ok",
            probability=predict_probability(model, abiotic_design),
            converged=model.converged, sampled_presences=length(presences),
            sampled_background=length(background))
end

function evaluate_model(
    design::Matrix{Float64},
    training_cells::Vector{Int},
    training_response::Vector{Float64},
    truth::BitVector
)
    model = fit_ridge_logistic(design[training_cells, :], training_response)
    probability = predict_probability(model, design)
    threshold = best_tss_threshold(probability[training_cells], training_response)
    binary = BitVector(probability .>= threshold)
    return (
        converged=model.converged,
        iterations=model.iterations,
        auc=auc_score(probability, truth),
        brier=mean((probability .- Float64.(truth)).^2),
        logloss=binary_log_loss(probability, truth),
        jaccard=jaccard_similarity(binary, truth),
        threshold=threshold
    )
end

function fit_prey_models(
    world,
    resources::Vector{Int},
    abiotic_design::Matrix{Float64},
    accessibility::Vector{Float64},
    seed::Int
)
    predictions = Dict{Int,Vector{Float64}}()
    diagnostics = NamedTuple[]
    for resource in sort(unique(resources))
        prey_rng = MersenneTwister(seed + 104729 * resource)
        fit = fit_sampled_distribution(
            prey_rng, world.AB[resource], abiotic_design, accessibility;
            n_presence=PREY_PRESENCE_SAMPLE_SIZE,
            minimum_presence_cells=MIN_PREY_PRESENCE_CELLS
        )
        if fit.success
            predictions[resource] = fit.probability
        end
        push!(diagnostics, (
            resource_id=resource,
            true_presence_cells=count(world.AB[resource]),
            sampled_presences=fit.sampled_presences,
            sampled_background=fit.sampled_background,
            model_available=fit.success,
            model_converged=fit.converged,
            failure_reason=fit.reason
        ))
    end
    return predictions, diagnostics
end

function model_record(
    model_name::String,
    link_completeness::Float64,
    known_links::Int,
    modeled_links::Int,
    predictor::Vector{Float64},
    metrics
)
    return (
        model=model_name,
        link_completeness=link_completeness,
        known_links=known_links,
        modeled_links=modeled_links,
        predictor_prevalence=mean(predictor),
        converged=metrics.converged,
        iterations=metrics.iterations,
        auc=metrics.auc,
        brier=metrics.brier,
        logloss=metrics.logloss,
        jaccard=metrics.jaccard,
        threshold=metrics.threshold
    )
end

function comparison_record(candidate, baseline)
    return (
        information_level=candidate.model,
        link_completeness=candidate.link_completeness,
        known_links=candidate.known_links,
        modeled_links=candidate.modeled_links,
        modeled_link_fraction=candidate.known_links == 0 ? 0.0 :
            candidate.modeled_links / candidate.known_links,
        predictor_prevalence=candidate.predictor_prevalence,
        abiotic_auc=baseline.auc,
        biotic_auc=candidate.auc,
        delta_auc=candidate.auc - baseline.auc,
        abiotic_brier=baseline.brier,
        biotic_brier=candidate.brier,
        delta_brier=baseline.brier - candidate.brier,
        abiotic_logloss=baseline.logloss,
        biotic_logloss=candidate.logloss,
        delta_logloss=baseline.logloss - candidate.logloss,
        abiotic_jaccard=baseline.jaccard,
        biotic_jaccard=candidate.jaccard,
        delta_jaccard=candidate.jaccard - baseline.jaccard,
        abiotic_converged=baseline.converged,
        biotic_converged=candidate.converged
    )
end

function run_focal_models(
    rng::AbstractRNG,
    world,
    focal_consumer::Int,
    abiotic_design::Matrix{Float64},
    accessibility::Vector{Float64},
    prey_predictions::Dict{Int,Vector{Float64}}
)
    truth = world.AB[focal_consumer]
    presences, background = sample_presence_background(
        rng, truth, accessibility;
        n_presence=FOCAL_PRESENCE_SAMPLE_SIZE,
        n_background=BACKGROUND_SAMPLE_SIZE
    )
    isempty(presences) && error("Focal consumer has no sampled presences")
    training_cells = vcat(presences, background)
    training_response = vcat(
        ones(Float64, length(presences)), zeros(Float64, length(background))
    )

    baseline_metrics = evaluate_model(
        abiotic_design, training_cells, training_response, truth
    )
    baseline = model_record(
        "abiotic", NaN, 0, 0, zeros(Float64, NCELLS), baseline_metrics
    )
    models = [baseline]

    resources = copy(world.prey[focal_consumer])
    shuffled_resources = shuffle(rng, resources)
    oracle_predictor = oracle_resource_availability(world, focal_consumer)
    oracle_metrics = evaluate_model(
        design_matrix(abiotic_design[:, 2]; biotic_predictor=oracle_predictor),
        training_cells, training_response, truth
    )
    oracle = model_record(
        "oracle", 1.0, length(resources), length(resources),
        oracle_predictor, oracle_metrics
    )
    push!(models, oracle)

    for completeness in LINK_COMPLETENESS_LEVELS
        number_known = max(1, min(length(resources), round(Int, completeness * length(resources))))
        known_resources = shuffled_resources[1:number_known]
        number_modeled = count(resource -> haskey(prey_predictions, resource), known_resources)
        estimated_predictor = estimated_resource_availability(
            prey_predictions, known_resources
        )
        metrics = evaluate_model(
            design_matrix(abiotic_design[:, 2]; biotic_predictor=estimated_predictor),
            training_cells, training_response, truth
        )
        label = "estimated_" * string(round(Int, 100 * completeness))
        push!(models, model_record(
            label, completeness, number_known, number_modeled,
            estimated_predictor, metrics
        ))
    end
    comparisons = [comparison_record(model, baseline) for model in models[2:end]]
    return (
        models=models,
        comparisons=comparisons,
        sampled_presences=length(presences),
        sampled_background=length(background)
    )
end

function build_jobs(;
    environments=ENVIRONMENTS,
    networks=NETWORKS,
    regime_indices=REGIME_INDICES,
    correlations=CORRELATIONS,
    replicates=N_REPLICATES
)
    jobs = NamedTuple[]
    community_id = 0
    for (environment_index, environment) in enumerate(environments),
        (network_index, network) in enumerate(networks),
        regime_index in regime_indices,
        (correlation_index, target_r) in enumerate(correlations),
        replicate in 1:replicates
        community_id += 1
        push!(jobs, (
            community_id=community_id,
            environment=environment,
            environment_index=environment_index,
            network=network,
            network_index=network_index,
            regime_index=regime_index,
            correlation_index=correlation_index,
            target_r=Float64(target_r),
            replicate=replicate,
            seed=BASE_SEED + 1_000_003 * community_id
        ))
    end
    return jobs
end

function common_metadata(job, world, focal_consumer::Int)
    return (
        community_id=job.community_id,
        replicate=job.replicate,
        environment=String(job.environment),
        network=String(job.network),
        regime=String(world.regime_name),
        target_r=world.target_r,
        achieved_r=world.achieved_r,
        focal_consumer=focal_consumer,
        degree=world.realized_degree[focal_consumer],
        true_mismatch=focal_truth_mismatch(world, focal_consumer),
        true_presence_cells=count(world.AB[focal_consumer]),
        true_resource_links=length(world.prey[focal_consumer])
    )
end

function run_world_sdms(job, workspace)
    simulation_rng = MersenneTwister(job.seed)
    world = simulate_world!(
        simulation_rng, workspace, job.environment, job.network,
        regimes[job.regime_index], job.target_r, job.community_id
    )
    accessibility = sampling_bias_surface()
    environment_z = standardized_environment(world.environment)
    abiotic_design = design_matrix(environment_z)

    selected_focals = Int[]
    exclusions = NamedTuple[]
    selection_rng = MersenneTwister(job.seed + 17)
    if world.correlation_ok
        for degree in DEGREES
            focals = choose_focal_consumers(selection_rng, world, degree)
            if isempty(focals)
                push!(exclusions, (
                    community_id=job.community_id,
                    replicate=job.replicate,
                    environment=String(job.environment),
                    network=String(job.network),
                    regime=String(world.regime_name),
                    target_r=world.target_r,
                    degree=degree,
                    reason="no_eligible_focal"
                ))
            else
                append!(selected_focals, focals)
            end
        end
    else
        for degree in DEGREES
            push!(exclusions, (
                community_id=job.community_id,
                replicate=job.replicate,
                environment=String(job.environment),
                network=String(job.network),
                regime=String(world.regime_name),
                target_r=world.target_r,
                degree=degree,
                reason="correlation_calibration_failed"
            ))
        end
    end

    required_resources = isempty(selected_focals) ? Int[] :
        unique(vcat((world.prey[focal] for focal in selected_focals)...))
    prey_predictions, prey_diagnostics = fit_prey_models(
        world, required_resources, abiotic_design, accessibility, job.seed + 37
    )

    model_rows = NamedTuple[]
    comparison_rows = NamedTuple[]
    focal_rows = NamedTuple[]
    for focal in selected_focals
        focal_rng = MersenneTwister(job.seed + 1_000_033 * focal)
        result = run_focal_models(
            focal_rng, world, focal, abiotic_design, accessibility, prey_predictions
        )
        metadata = common_metadata(job, world, focal)
        append!(model_rows, [merge(metadata, row) for row in result.models])
        append!(comparison_rows, [merge(metadata, row) for row in result.comparisons])
        push!(focal_rows, merge(metadata, (
            sampled_presences=result.sampled_presences,
            sampled_background=result.sampled_background,
            prey_models_available=count(resource -> haskey(prey_predictions, resource),
                                        world.prey[focal])
        )))
    end

    community_row = (
        community_id=job.community_id,
        replicate=job.replicate,
        environment=String(job.environment),
        network=String(job.network),
        regime=String(world.regime_name),
        target_r=world.target_r,
        achieved_r=world.achieved_r,
        correlation_ok=world.correlation_ok,
        max_degree_correlation_error=world.max_degree_correlation_error,
        eligible_focals=length(selected_focals),
        unique_required_resources=length(required_resources),
        prey_models_available=length(prey_predictions),
        prey_models_converged=count(row -> row.model_converged, prey_diagnostics)
    )
    return (
        community_rows=[community_row],
        focal_rows=focal_rows,
        model_rows=model_rows,
        comparison_rows=comparison_rows,
        exclusions=exclusions
    )
end

function checkpoint_path(job, checkpoint_dir::String)
    return joinpath(checkpoint_dir, "community_" * lpad(job.community_id, 6, '0') * ".jls")
end

function combine_world_results(world_results)
    combined = (
        community_rows=NamedTuple[], focal_rows=NamedTuple[],
        model_rows=NamedTuple[], comparison_rows=NamedTuple[],
        exclusions=NamedTuple[]
    )
    for result in world_results
        append!(combined.community_rows, result.community_rows)
        append!(combined.focal_rows, result.focal_rows)
        append!(combined.model_rows, result.model_rows)
        append!(combined.comparison_rows, result.comparison_rows)
        append!(combined.exclusions, result.exclusions)
    end
    return combined
end

function run_experiment(;
    environments=ENVIRONMENTS,
    networks=NETWORKS,
    regime_indices=REGIME_INDICES,
    correlations=CORRELATIONS,
    replicates=N_REPLICATES,
    checkpoint_dir=CHECKPOINT_DIR,
    resume::Bool=true
)
    jobs = build_jobs(
        environments=environments, networks=networks,
        regime_indices=regime_indices, correlations=correlations,
        replicates=replicates
    )
    mkpath(checkpoint_dir)
    workspaces = make_workspaces()
    results = Vector{Any}(undef, length(jobs))
    completed = Threads.Atomic{Int}(0)
    Threads.@threads for index in eachindex(jobs)
        job = jobs[index]
        path = checkpoint_path(job, checkpoint_dir)
        results[index] = if resume && isfile(path)
            deserialize(path)
        else
            result = run_world_sdms(job, workspaces[Threads.threadid()])
            serialize(path, result)
            result
        end
        current = Threads.atomic_add!(completed, 1) + 1
        if current == length(jobs) || current % max(1, length(jobs) ÷ 100) == 0
            println("Completed $current / $(length(jobs)) communities")
        end
    end
    return combine_world_results(results)
end

function finite_mean(values)
    valid = filter(isfinite, Float64.(values))
    return isempty(valid) ? NaN : mean(valid)
end

function finite_se(values)
    valid = filter(isfinite, Float64.(values))
    return length(valid) <= 1 ? NaN : std(valid) / sqrt(length(valid))
end

function summarize_comparisons(rows)
    groups = Dict{Tuple,Vector{NamedTuple}}()
    for row in rows
        key = (row.environment, row.network, row.regime, row.target_r,
               row.degree, row.information_level, row.link_completeness)
        push!(get!(groups, key, NamedTuple[]), row)
    end
    summary = NamedTuple[]
    for (key, group) in groups
        environment, network, regime, target_r, degree,
            information_level, link_completeness = key
        push!(summary, (
            environment=environment, network=network, regime=regime,
            target_r=target_r, degree=degree,
            information_level=information_level,
            link_completeness=link_completeness,
            n=length(group),
            mean_achieved_r=finite_mean(getproperty.(group, :achieved_r)),
            mean_true_mismatch=finite_mean(getproperty.(group, :true_mismatch)),
            mean_modeled_link_fraction=finite_mean(getproperty.(group, :modeled_link_fraction)),
            mean_delta_auc=finite_mean(getproperty.(group, :delta_auc)),
            se_delta_auc=finite_se(getproperty.(group, :delta_auc)),
            mean_delta_brier=finite_mean(getproperty.(group, :delta_brier)),
            se_delta_brier=finite_se(getproperty.(group, :delta_brier)),
            mean_delta_logloss=finite_mean(getproperty.(group, :delta_logloss)),
            se_delta_logloss=finite_se(getproperty.(group, :delta_logloss)),
            mean_delta_jaccard=finite_mean(getproperty.(group, :delta_jaccard)),
            se_delta_jaccard=finite_se(getproperty.(group, :delta_jaccard))
        ))
    end
    sort!(summary, by=row -> (row.environment, row.network, row.regime,
                              row.information_level, row.target_r, row.degree))
    return summary
end

function summarize_relationship(rows; bin_width::Float64=MISMATCH_BIN_WIDTH)
    groups = Dict{Tuple,Vector{NamedTuple}}()
    for row in rows
        isfinite(row.true_mismatch) || continue
        bin = clamp(floor(Int, row.true_mismatch / bin_width), 0, floor(Int, 1 / bin_width))
        midpoint = min(1.0, (bin + 0.5) * bin_width)
        key = (row.information_level, row.link_completeness, midpoint)
        push!(get!(groups, key, NamedTuple[]), row)
    end
    summary = NamedTuple[]
    for (key, group) in groups
        information_level, link_completeness, mismatch_midpoint = key
        push!(summary, (
            information_level=information_level,
            link_completeness=link_completeness,
            mismatch_midpoint=mismatch_midpoint,
            n=length(group),
            mean_true_mismatch=finite_mean(getproperty.(group, :true_mismatch)),
            mean_delta_auc=finite_mean(getproperty.(group, :delta_auc)),
            se_delta_auc=finite_se(getproperty.(group, :delta_auc)),
            mean_delta_brier=finite_mean(getproperty.(group, :delta_brier)),
            se_delta_brier=finite_se(getproperty.(group, :delta_brier)),
            mean_delta_jaccard=finite_mean(getproperty.(group, :delta_jaccard)),
            se_delta_jaccard=finite_se(getproperty.(group, :delta_jaccard))
        ))
    end
    sort!(summary, by=row -> (row.information_level, row.mismatch_midpoint))
    return summary
end

end
