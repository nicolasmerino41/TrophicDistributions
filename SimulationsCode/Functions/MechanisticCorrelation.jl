module MechanisticCorrelation

using LinearAlgebra
using Random
using Statistics
using ..Parameters: S, E_MIN, E_MAX, DEGREE_CLASSES,
                    CORR_CALIBRATION_ITERS, CORR_CALIBRATION_RESTARTS,
                    CORR_CALIBRATION_DAMPING, CORR_WITHIN_DEGREE_SD,
                    TARGET_R_TOL,
                    DEGREE_TARGET_R_TOL

export pearson_r, prey_means, mechanistic_corr, degree_corr_diagnostics,
       target_within_tolerance, degree_targets_within_tolerance,
       max_degree_target_error, assign_mus_with_target_corr!

function pearson_r(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    length(a) == length(b) || error("Correlation vectors must have equal length")
    length(a) < 2 && return NaN
    ac = Float64.(a) .- mean(a)
    bc = Float64.(b) .- mean(b)
    denominator = sqrt(sum(abs2, ac) * sum(abs2, bc))
    denominator <= eps(Float64) && return 0.0
    return clamp(dot(ac, bc) / denominator, -1.0, 1.0)
end

function prey_means(mu::Vector{Float64}, prey::Vector{Vector{Int}}, basal_mask::BitVector)
    pm = fill(NaN, S)
    for i in 1:S
        if basal_mask[i] || isempty(prey[i])
            continue
        end
        pm[i] = mean(mu[prey[i]])
    end
    return pm
end

function degree_groups(
    prey::Vector{Vector{Int}},
    basal_mask::BitVector;
    degree_classes::Vector{Int}=DEGREE_CLASSES
)
    return Dict(
        degree => findall(
            i -> !basal_mask[i] && length(prey[i]) == degree,
            1:S
        )
        for degree in degree_classes
    )
end

function correlations_from_prey_means(
    mu::Vector{Float64},
    pm::Vector{Float64},
    basal_mask::BitVector,
    groups::Dict{Int,Vector{Int}}
)
    consumers = findall(i -> !basal_mask[i] && isfinite(pm[i]), 1:S)
    overall = length(consumers) < 2 ? NaN : pearson_r(mu[consumers], pm[consumers])
    by_degree = Dict{Int,Float64}()
    for (degree, idx) in groups
        valid = filter(i -> isfinite(pm[i]), idx)
        by_degree[degree] = length(valid) < 3 ? NaN : pearson_r(mu[valid], pm[valid])
    end
    return overall, by_degree
end

function mechanistic_corr(mu::Vector{Float64}, prey::Vector{Vector{Int}}, basal_mask::BitVector)
    pm = prey_means(mu, prey, basal_mask)
    groups = degree_groups(prey, basal_mask)
    overall, _ = correlations_from_prey_means(mu, pm, basal_mask, groups)
    return overall
end

function degree_corr_diagnostics(
    mu::Vector{Float64},
    prey::Vector{Vector{Int}},
    basal_mask::BitVector;
    degree_classes::Vector{Int}=DEGREE_CLASSES
)
    pm = prey_means(mu, prey, basal_mask)
    groups = degree_groups(prey, basal_mask; degree_classes=degree_classes)
    _, diagnostics = correlations_from_prey_means(mu, pm, basal_mask, groups)
    return diagnostics
end

target_within_tolerance(achieved_r::Float64, target_r::Float64) =
    isfinite(achieved_r) && abs(achieved_r - target_r) <= TARGET_R_TOL

function max_degree_target_error(diagnostics::Dict{Int,Float64}, target_r::Float64)
    values_r = collect(values(diagnostics))
    any(r -> !isfinite(r), values_r) && return Inf
    return maximum(abs(r - target_r) for r in values_r)
end

degree_targets_within_tolerance(diagnostics::Dict{Int,Float64}, target_r::Float64) =
    max_degree_target_error(diagnostics, target_r) <= DEGREE_TARGET_R_TOL

function orthogonal_residual(raw::Vector{Float64}, predictor_unit::Vector{Float64})
    residual = raw .- mean(raw)
    residual .-= dot(residual, predictor_unit) .* predictor_unit
    residual_norm = norm(residual)
    if residual_norm <= 1e-10
        residual = collect(range(-1.0, 1.0, length=length(raw)))
        residual .-= mean(residual)
        residual .-= dot(residual, predictor_unit) .* predictor_unit
        residual_norm = norm(residual)
    end
    residual_norm <= 1e-10 && error("Could not construct an orthogonal niche residual")
    return residual ./ residual_norm
end

function weighted_orthogonal_group_residual(
    raw::Vector{Float64},
    weights::Vector{Float64},
    predictor::Vector{Float64}
)
    residual = raw .- sum(weights .* raw) / sum(weights)
    predictor_ss = sum(weights .* predictor.^2)
    if predictor_ss > 1e-12
        residual .-= (sum(weights .* residual .* predictor) / predictor_ss) .* predictor
    end
    residual_norm = sqrt(sum(weights .* residual.^2))
    if residual_norm <= 1e-10
        residual = collect(range(-1.0, 1.0, length=length(raw)))
        residual .-= sum(weights .* residual) / sum(weights)
        if predictor_ss > 1e-12
            residual .-= (sum(weights .* residual .* predictor) / predictor_ss) .* predictor
        end
        residual_norm = sqrt(sum(weights .* residual.^2))
    end
    residual_norm <= 1e-10 && error("Could not construct a group-level niche residual")
    return residual ./ residual_norm
end

function joint_projected_consumer_optima(
    pm::Vector{Float64},
    consumers::Vector{Int},
    groups::Dict{Int,Vector{Int}},
    residuals::Dict{Int,Vector{Float64}},
    group_residual::Dict{Int,Float64},
    target_r::Float64,
    desired_sd::Float64
)
    global_prey_mean = mean(pm[consumers])
    signal = pm[consumers] .- global_prey_mean
    residual = zeros(Float64, length(consumers))
    consumer_position = Dict(species => position for (position, species) in enumerate(consumers))

    if abs(target_r) <= 1e-10
        for (degree, idx) in groups
            predictor = pm[idx] .- mean(pm[idx])
            predictor_norm = norm(predictor)
            predictor_norm <= 1e-10 && error("Prey means have no variance in degree class $degree")
            within = orthogonal_residual(residuals[degree], predictor ./ predictor_norm)
            for (local_position, species) in enumerate(idx)
                residual[consumer_position[species]] = within[local_position]
            end
        end
        combined = residual
    else
        ratio = sqrt(1.0 - target_r^2) / abs(target_r)
        for (degree, idx) in groups
            predictor = pm[idx] .- mean(pm[idx])
            predictor_norm = norm(predictor)
            predictor_norm <= 1e-10 && error("Prey means have no variance in degree class $degree")
            within = orthogonal_residual(residuals[degree], predictor ./ predictor_norm)
            within .*= ratio * predictor_norm
            for (local_position, species) in enumerate(idx)
                residual[consumer_position[species]] = within[local_position]
            end
        end

        degrees = collect(DEGREE_CLASSES)
        weights = Float64[length(groups[degree]) for degree in degrees]
        group_predictor = Float64[
            mean(pm[groups[degree]]) - global_prey_mean for degree in degrees
        ]
        between_ss = sum(weights .* group_predictor.^2)
        if between_ss > 1e-12
            raw_group = Float64[group_residual[degree] for degree in degrees]
            group_unit = weighted_orthogonal_group_residual(
                raw_group, weights, group_predictor
            )
            group_unit .*= ratio * sqrt(between_ss)
            for (group_position, degree) in enumerate(degrees), species in groups[degree]
                residual[consumer_position[species]] += group_unit[group_position]
            end
        end
        combined = sign(target_r) .* signal .+ residual
    end

    combined .-= mean(combined)
    combined_sd = std(combined)
    combined_sd <= 1e-10 && error("Projected consumer optima have no variance")
    combined ./= combined_sd

    center = 0.5 * (E_MIN + E_MAX)
    feasible_sd = Inf
    for value in combined
        if value > 0
            feasible_sd = min(feasible_sd, (E_MAX - center) / value)
        elseif value < 0
            feasible_sd = min(feasible_sd, (E_MIN - center) / value)
        end
    end
    applied_sd = min(desired_sd, 0.98 * feasible_sd)
    return center .+ applied_sd .* combined
end

"""
Calibrate consumer optima so every degree class and the pooled community target
the same consumer-prey niche correlation. Degree classes are projected
simultaneously, then prey means are recomputed until the coupled food web
converges. The best solution across restarts is retained.
"""
function assign_mus_with_target_corr!(
    rng::AbstractRNG,
    mu::Vector{Float64},
    prey::Vector{Vector{Int}},
    basal_mask::BitVector,
    target_r::Float64;
    calibration_iters::Int=CORR_CALIBRATION_ITERS,
    restarts::Int=CORR_CALIBRATION_RESTARTS,
    damping::Float64=CORR_CALIBRATION_DAMPING
)
    -1.0 <= target_r <= 1.0 || error("target_r must be between -1 and 1")
    0.0 < damping <= 1.0 || error("Calibration damping must be in (0, 1]")

    consumers = findall(!, basal_mask)
    groups = degree_groups(prey, basal_mask)
    any(isempty, values(groups)) && error("Every configured degree class must be represented")

    original_mu = copy(mu)
    best_mu = copy(mu)
    best_overall = mechanistic_corr(mu, prey, basal_mask)
    best_score = Inf

    for restart in 1:restarts
        current = copy(original_mu)
        if restart > 1
            for i in consumers
                current[i] = rand(rng) * (E_MAX - E_MIN) + E_MIN
            end
        end

        residuals = Dict(
            degree => randn(rng, length(idx))
            for (degree, idx) in groups
        )
        group_residual = Dict(degree => randn(rng) for degree in keys(groups))
        desired_sd = CORR_WITHIN_DEGREE_SD

        for _ in 1:calibration_iters
            pm = prey_means(current, prey, basal_mask)
            projected = joint_projected_consumer_optima(
                pm, consumers, groups, residuals, group_residual,
                target_r, desired_sd
            )
            current[consumers] .= (1.0 - damping) .* current[consumers] .+
                                  damping .* projected

            updated_pm = prey_means(current, prey, basal_mask)
            overall, by_degree = correlations_from_prey_means(
                current, updated_pm, basal_mask, groups
            )
            degree_error = max_degree_target_error(by_degree, target_r)
            overall_error = isfinite(overall) ? abs(overall - target_r) : Inf
            score = max(degree_error, overall_error)

            if score < best_score
                best_score = score
                best_mu = copy(current)
                best_overall = overall
            end
            if degree_error <= DEGREE_TARGET_R_TOL && overall_error <= TARGET_R_TOL
                break
            end
        end
        best_score <= min(DEGREE_TARGET_R_TOL, TARGET_R_TOL) && break
    end

    copyto!(mu, best_mu)
    return best_overall
end

end
