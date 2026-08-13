module Simulation

using Random
using Statistics
using ..Parameters: S, SUIT_THRESH, Emin_patch, E_MIN, E_MAX, TAIL_THRESH
using ..Environment: make_environment
using ..Niches: suitability_mask_1d, draw_sigmas, BreadthRegime
using ..Networks: consumers_and_basal, assign_balanced_degrees,
                  realized_degrees, build_trophic_community
using ..MechanisticCorrelation: assign_mus_with_target_corr!,
                                degree_corr_diagnostics,
                                target_within_tolerance,
                                degree_targets_within_tolerance,
                                max_degree_target_error
using ..Dynamics: trophic_distributions
using ..Connectivity: apply_connectivity_filter, CCWorkspace
using ..Metrics: gamma_richness_cons, frac_affected, realized_overlap,
                 realized_overlap_by_species, jaccard_mismatch_by_species,
                 mismatch_q90, mismatch_frac_gt

export simulate_world!, summarize_world, simulate_one!

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
    regime::BreadthRegime,
    target_r::Float64,
    community_id::Int
)
    _, basal_mask, consumers = consumers_and_basal()
    target_degree = assign_balanced_degrees(rng, basal_mask)
    environment = make_environment(rng, envkind)

    prey = build_trophic_community(rng, basal_mask, target_degree)
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
    AB_direct_raw, AB_raw = trophic_distributions(A_raw, prey, basal_mask)

    A = Vector{BitVector}(undef, S)
    AB_direct = Vector{BitVector}(undef, S)
    AB = Vector{BitVector}(undef, S)
    for i in 1:S
        A[i] = apply_connectivity_filter(ws, A_raw[i], Emin_patch)
        AB_direct[i] = apply_connectivity_filter(ws, AB_direct_raw[i], Emin_patch)
        AB[i] = apply_connectivity_filter(ws, AB_raw[i], Emin_patch)
    end

    return (
        community_id=community_id,
        envkind=envkind,
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
        AB_direct_raw=AB_direct_raw,
        AB_raw=AB_raw,
        A=A,
        AB_direct=AB_direct,
        AB=AB
    )
end

function summarize_world(world)
    mismatch, eligible = jaccard_mismatch_by_species(
        world.A, world.AB, world.basal_mask
    )
    direct_mismatch, direct_eligible = jaccard_mismatch_by_species(
        world.A, world.AB_direct, world.basal_mask
    )
    eligible == direct_eligible || error("Direct and recursive eligibility differ")
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
        direct_mismatch=direct_mismatch[i],
        indirect_increment=eligible[i] ? mismatch[i] - direct_mismatch[i] : NaN,
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
        mean_direct_mismatch=isempty(eligible_mismatch) ? NaN : mean(direct_mismatch[eligible]),
        mean_indirect_increment=isempty(eligible_mismatch) ? NaN :
            mean(mismatch[eligible] .- direct_mismatch[eligible]),
        frac_affected=frac_affected(world.A, world.AB, world.basal_mask),
        realized_overlap=realized_overlap(world.A_raw, world.prey, world.basal_mask),
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
    regime::BreadthRegime,
    target_r::Float64,
    community_id::Int
)
    world = simulate_world!(
        rng, ws, envkind, regime, target_r, community_id
    )
    return summarize_world(world)
end

end
