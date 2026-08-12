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
       sampling_bias_surface, true_resource_availability,
       choose_focal_consumer, run_focal_sdm, build_jobs,
       run_world_sdms, run_experiment, summarize_results

struct LogisticModel
    coefficients::Vector{Float64}
    converged::Bool
    iterations::Int
end

@inline function logistic(value::Float64)
    value >= 0 && return 1.0 / (1.0 + exp(-value))
    exp_value = exp(value)
    return exp_value / (1.0 + exp_value)
end

predict_probability(model::LogisticModel, design::Matrix{Float64}) =
    logistic.(design * model.coefficients)

function fit_ridge_logistic(
    design::Matrix{Float64}, response::Vector{Float64};
    ridge::Float64=LOGISTIC_RIDGE,
    max_iters::Int=LOGISTIC_MAX_ITERS,
    tolerance::Float64=LOGISTIC_TOL
)
    size(design, 1) == length(response) || error("Response and design lengths differ")
    ncol = size(design, 2)
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

function design_matrix(environment_z::Vector{Float64}; resource=nothing)
    abiotic = hcat(ones(length(environment_z)), environment_z, environment_z.^2)
    resource === nothing && return abiotic
    return hcat(abiotic, resource)
end

function sampling_bias_surface()
    center_x = 0.25 * NX
    center_y = 0.75 * NY
    return [
        SAMPLING_BIAS_FLOOR + exp(-SAMPLING_BIAS_STRENGTH *
            (((x_of(cell) - center_x) / NX)^2 + ((y_of(cell) - center_y) / NY)^2))
        for cell in 1:NCELLS
    ]
end

function weighted_sample_without_replacement(
    rng::AbstractRNG, candidates::Vector{Int}, weights::Vector{Float64}, requested::Int
)
    sample_size = min(requested, length(candidates))
    pool = copy(candidates)
    pool_weights = [max(weights[cell], 0.0) for cell in pool]
    selected = Int[]
    for _ in 1:sample_size
        total = sum(pool_weights)
        chosen = if total <= 0
            rand(rng, eachindex(pool))
        else
            draw = rand(rng) * total
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
    rng::AbstractRNG, truth::BitVector, accessibility::Vector{Float64}
)
    presences = weighted_sample_without_replacement(
        rng, findall(identity, truth), accessibility, FOCAL_PRESENCE_SAMPLE_SIZE
    )
    sampled_presence = falses(length(truth))
    sampled_presence[presences] .= true
    background = weighted_sample_without_replacement(
        rng, findall(!, sampled_presence), accessibility, BACKGROUND_SAMPLE_SIZE
    )
    return presences, background
end

function true_resource_availability(world, focal_consumer::Int)
    availability = falses(NCELLS)
    for resource in world.prey[focal_consumer]
        availability .|= world.AB[resource]
    end
    return Float64.(availability)
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
        while final_position < length(order) &&
              probability[order[final_position + 1]] == probability[order[position]]
            final_position += 1
        end
        average_rank = 0.5 * (position + final_position)
        for rank_position in position:final_position
            ranks[order[rank_position]] = average_rank
        end
        position = final_position + 1
    end
    return (sum(ranks[truth]) - positives * (positives + 1) / 2) /
           (positives * negatives)
end

function focal_truth_mismatch(world, focal_consumer::Int)
    abiotic = world.A[focal_consumer]
    realized = world.AB[focal_consumer]
    union_size = count(abiotic .| realized)
    return union_size == 0 ? NaN : 1.0 - count(abiotic .& realized) / union_size
end

function choose_focal_consumer(rng::AbstractRNG, world, degree::Int)
    candidates = [
        consumer for consumer in world.consumers
        if world.realized_degree[consumer] == degree &&
           count(world.AB[consumer]) >= MIN_FOCAL_PRESENCE_CELLS
    ]
    isempty(candidates) && return nothing
    return rand(rng, candidates)
end

function run_focal_sdm(
    rng::AbstractRNG, world, focal_consumer::Int,
    environment_z::Vector{Float64}, accessibility::Vector{Float64}
)
    truth = world.AB[focal_consumer]
    presences, background = sample_presence_background(rng, truth, accessibility)
    training_cells = vcat(presences, background)
    response = vcat(ones(Float64, length(presences)), zeros(Float64, length(background)))
    resource = true_resource_availability(world, focal_consumer)

    abiotic_design = design_matrix(environment_z)
    biotic_design = design_matrix(environment_z; resource=resource)
    abiotic_model = fit_ridge_logistic(abiotic_design[training_cells, :], response)
    biotic_model = fit_ridge_logistic(biotic_design[training_cells, :], response)
    abiotic_probability = predict_probability(abiotic_model, abiotic_design)
    biotic_probability = predict_probability(biotic_model, biotic_design)
    abiotic_auc = auc_score(abiotic_probability, truth)
    biotic_auc = auc_score(biotic_probability, truth)
    numeric_truth = Float64.(truth)
    abiotic_brier = mean((abiotic_probability .- numeric_truth).^2)
    biotic_brier = mean((biotic_probability .- numeric_truth).^2)

    return (
        true_mismatch=focal_truth_mismatch(world, focal_consumer),
        true_presence_cells=count(truth),
        true_resource_links=length(world.prey[focal_consumer]),
        predictor_prevalence=mean(resource),
        abiotic_auc=abiotic_auc,
        biotic_auc=biotic_auc,
        delta_auc=biotic_auc - abiotic_auc,
        abiotic_brier=abiotic_brier,
        biotic_brier=biotic_brier,
        delta_brier=abiotic_brier - biotic_brier,
        abiotic_converged=abiotic_model.converged,
        biotic_converged=biotic_model.converged
    )
end

function build_jobs(;
    environments=ENVIRONMENTS, regime_indices=REGIME_INDICES,
    correlations=CORRELATIONS,
    replicates=N_REPLICATES
)
    jobs = NamedTuple[]
    community_id = 0
    for (environment_index, environment) in enumerate(environments),
        regime_index in regime_indices,
        target_r in correlations,
        replicate in 1:replicates
        community_id += 1
        push!(jobs, (
            community_id=community_id, environment=environment,
            regime_index=regime_index, target_r=Float64(target_r),
            replicate=replicate,
            seed=BASE_SEED + 10_000_000 * environment_index +
                 100_000 * regime_index + replicate
        ))
    end
    return jobs
end

function run_world_sdms(job, workspace)
    world = simulate_world!(
        MersenneTwister(job.seed), workspace, job.environment,
        regimes[job.regime_index], job.target_r, job.community_id
    )
    rows = NamedTuple[]
    exclusions = NamedTuple[]
    world.correlation_ok || return (rows=rows, exclusions=[(
        community_id=job.community_id, environment=String(job.environment),
        regime=String(world.regime_name),
        target_r=world.target_r, degree=0, reason="correlation_calibration_failed"
    )])

    environment_z = standardized_environment(world.environment)
    accessibility = sampling_bias_surface()
    selection_rng = MersenneTwister(job.seed + 17)
    for degree in DEGREES
        focal = choose_focal_consumer(selection_rng, world, degree)
        if focal === nothing
            push!(exclusions, (
                community_id=job.community_id, environment=String(job.environment),
                regime=String(world.regime_name),
                target_r=world.target_r, degree=degree, reason="no_eligible_focal"
            ))
            continue
        end
        result = run_focal_sdm(
            MersenneTwister(job.seed + 1_000_033 * focal), world, focal,
            environment_z, accessibility
        )
        push!(rows, merge((
            community_id=job.community_id, replicate=job.replicate,
            environment=String(job.environment),
            regime=String(world.regime_name), target_r=world.target_r,
            achieved_r=world.achieved_r, focal_consumer=focal,
            degree=world.realized_degree[focal]
        ), result))
    end
    return (rows=rows, exclusions=exclusions)
end

function checkpoint_path(job, directory)
    joinpath(directory, "community_" * lpad(job.community_id, 6, '0') * ".jls")
end

function run_experiment(;
    environments=ENVIRONMENTS, regime_indices=REGIME_INDICES,
    correlations=CORRELATIONS,
    replicates=N_REPLICATES, checkpoint_dir=CHECKPOINT_DIR, resume=true
)
    jobs = build_jobs(
        environments=environments, regime_indices=regime_indices,
        correlations=correlations,
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
        current % max(1, length(jobs) ÷ 100) == 0 &&
            println("Completed $current / $(length(jobs)) communities")
    end
    return (
        rows=vcat((result.rows for result in results)...),
        exclusions=vcat((result.exclusions for result in results)...)
    )
end

finite_mean(values) = mean(filter(isfinite, Float64.(values)))

function finite_se(values)
    valid = filter(isfinite, Float64.(values))
    length(valid) <= 1 && return NaN
    return std(valid) / sqrt(length(valid))
end

function summarize_results(rows)
    groups = Dict{Tuple,Vector{NamedTuple}}()
    for row in rows
        key = (row.environment, row.regime, row.target_r, row.degree)
        push!(get!(groups, key, NamedTuple[]), row)
    end
    summary = NamedTuple[]
    for (key, group) in groups
        environment, regime, target_r, degree = key
        push!(summary, (
            environment=environment, regime=regime,
            target_r=target_r, degree=degree, n=length(group),
            mean_achieved_r=finite_mean(getproperty.(group, :achieved_r)),
            mean_true_mismatch=finite_mean(getproperty.(group, :true_mismatch)),
            mean_delta_auc=finite_mean(getproperty.(group, :delta_auc)),
            se_delta_auc=finite_se(getproperty.(group, :delta_auc)),
            mean_delta_brier=finite_mean(getproperty.(group, :delta_brier)),
            se_delta_brier=finite_se(getproperty.(group, :delta_brier))
        ))
    end
    sort!(summary, by=row -> (row.environment, row.regime,
                              row.target_r, row.degree))
    return summary
end

end
