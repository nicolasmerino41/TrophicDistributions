module Simulation

using Random
using Statistics
using ..Parameters: S, SUIT_THRESH, Emin_patch, E_MIN, E_MAX, TAIL_THRESH
using ..Environment: make_environment
using ..Niches: suitability_mask_1d, draw_sigmas, BreadthRegime
using ..Networks: consumers_and_basal, assign_balanced_degrees,
                  realized_degrees, realized_connectance,
                  build_metaweb_random, build_metaweb_modular,
                  build_metaweb_cascade
using ..MechanisticCorrelation: assign_mus_with_target_corr!,
                                degree_corr_diagnostics,
                                target_within_tolerance,
                                degree_targets_within_tolerance,
                                max_degree_target_error
using ..Dynamics: fixed_point_AB
using ..Connectivity: apply_connectivity_filter, CCWorkspace
using ..Metrics: gamma_richness_cons, frac_affected, realized_overlap,
                 realized_overlap_by_species, jaccard_mismatch_by_species,
                 mismatch_q90, mismatch_frac_gt

export simulate_world!, summarize_world, simulate_one!, count_links

@inline function count_links(prey::Vector{Vector{Int}}, basal_mask::BitVector)
    return sum(length(prey[i]) for i in 1:S if !basal_mask[i])
end

"""
Generate one complete virtual community and retain the spatial truth required by
downstream applications. The standard sweep calls `summarize_world` immediately;
the SDM application consumes the returned maps directly without serializing the
full set of species distributions.
"""
function simulate_world!(
    rng::AbstractRNG,
    ws::CCWorkspace,
    envkind::Symbol,
    netfamily::Symbol,
    regime::BreadthRegime,
    target_r::Float64,
    community_id::Int
)
    _, basal_mask, consumers = consumers_and_basal()
    target_degree = assign_balanced_degrees(rng, basal_mask)
    environment = make_environment(rng, envkind)

    prey = if netfamily == :random
        build_metaweb_random(rng, basal_mask, target_degree)
    elseif netfamily == :modular
        first(build_metaweb_modular(rng, basal_mask, target_degree))
    elseif netfamily == :cascade
        first(build_metaweb_cascade(rng, basal_mask, target_degree))
    else
        error("Unknown netfamily: $netfamily")
    end
    realized_degree = realized_degrees(prey)

    sigma = draw_sigmas(rng, regime)
    mu = [rand(rng) * (E_MAX - E_MIN) + E_MIN for _ in 1:S]
    achieved_r = assign_mus_with_target_corr!(rng, mu, prey, basal_mask, target_r)
    corr_by_degree = degree_corr_diagnostics(mu, prey, basal_mask)
    overall_correlation_ok = target_within_tolerance(achieved_r, target_r)
    degree_correlation_ok = degree_targets_within_tolerance(corr_by_degree, target_r)

    A_raw = [
        suitability_mask_1d(environment, mu[i], sigma[i], SUIT_THRESH)
        for i in 1:S
    ]
    AB_raw = fixed_point_AB(A_raw, prey, basal_mask)

    A = Vector{BitVector}(undef, S)
    AB = Vector{BitVector}(undef, S)
    for i in 1:S
        A[i] = apply_connectivity_filter(ws, A_raw[i], Emin_patch)
        AB[i] = apply_connectivity_filter(ws, AB_raw[i], Emin_patch)
    end

    return (
        community_id=community_id,
        envkind=envkind,
        netfamily=netfamily,
        regime_name=regime.name,
        target_r=target_r,
        achieved_r=achieved_r,
        correlation_ok=overall_correlation_ok && degree_correlation_ok,
        overall_correlation_ok=overall_correlation_ok,
        degree_correlation_ok=degree_correlation_ok,
        max_degree_correlation_error=max_degree_target_error(corr_by_degree, target_r),
        degree_correlations=corr_by_degree,
        environment=environment,
        basal_mask=basal_mask,
        consumers=consumers,
        prey=prey,
        target_degree=target_degree,
        realized_degree=realized_degree,
        mu=mu,
        sigma=sigma,
        A_raw=A_raw,
        AB_raw=AB_raw,
        A=A,
        AB=AB
    )
end

function summarize_world(world)
    mismatch, eligible = jaccard_mismatch_by_species(
        world.A, world.AB, world.basal_mask
    )
    overlap_by_species = realized_overlap_by_species(
        world.A_raw, world.prey, world.basal_mask
    )
    eligible_mismatch = mismatch[eligible]

    SA = gamma_richness_cons(world.A, world.basal_mask)
    SAB = gamma_richness_cons(world.AB, world.basal_mask)
    dSrel = SA == 0 ? NaN : 1.0 - SAB / SA

    consumer_results = [(
        community_id=world.community_id,
        consumer_id=i,
        target_degree=world.target_degree[i],
        realized_degree=world.realized_degree[i],
        eligible=eligible[i],
        mismatch=mismatch[i],
        realized_overlap=overlap_by_species[i]
    ) for i in world.consumers]

    community_metrics = (
        community_id=world.community_id,
        target_r=world.target_r,
        achieved_r=world.achieved_r,
        correlation_ok=world.correlation_ok,
        overall_correlation_ok=world.overall_correlation_ok,
        degree_correlation_ok=world.degree_correlation_ok,
        max_degree_correlation_error=world.max_degree_correlation_error,
        dSrel=dSrel,
        mean_jaccard_mismatch=isempty(eligible_mismatch) ? NaN : mean(eligible_mismatch),
        frac_affected=frac_affected(world.A, world.AB, world.basal_mask),
        realized_overlap=realized_overlap(world.A_raw, world.prey, world.basal_mask),
        realized_connectance=realized_connectance(world.prey, world.basal_mask),
        mismatch_q90=mismatch_q90(eligible_mismatch),
        mismatch_frac_gt=mismatch_frac_gt(eligible_mismatch, TAIL_THRESH),
        n_eligible=count(eligible)
    )

    return (
        consumers=consumer_results,
        community=community_metrics,
        degree_correlations=world.degree_correlations
    )
end

function simulate_one!(
    rng::AbstractRNG,
    ws::CCWorkspace,
    envkind::Symbol,
    netfamily::Symbol,
    regime::BreadthRegime,
    target_r::Float64,
    community_id::Int
)
    world = simulate_world!(
        rng, ws, envkind, netfamily, regime, target_r, community_id
    )
    return summarize_world(world)
end

end
