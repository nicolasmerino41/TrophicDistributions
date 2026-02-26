module Simulation

using Random
using ..Parameters: S, SUIT_THRESH, Emin_patch, E_MIN, E_MAX, TAIL_THRESH
using ..Environment: make_environment
using ..Niches: suitability_mask_1d, draw_sigmas, BreadthRegime
using ..Networks: consumers_and_basal, build_metaweb_random, build_metaweb_modular, build_metaweb_heavytail, build_metaweb_cascade
using ..MechanisticCorrelation: assign_mus_with_target_corr!
using ..Dynamics: fixed_point_AB
using ..Connectivity: apply_connectivity_filter, CCWorkspace
using ..Metrics: gamma_richness_cons, mean_jaccard_mismatch, frac_affected, realized_overlap, jaccard_mismatch_vec, mismatch_q90, mismatch_frac_gt

export simulate_one!, count_links

@inline function count_links(prey::Vector{Vector{Int}}, basal_mask::BitVector)
    L = 0
    for i in 1:S
        basal_mask[i] && continue
        L += length(prey[i])
    end
    return L
end

function simulate_one!(
    rng::AbstractRNG,
    ws::CCWorkspace,
    envkind::Symbol,
    netfamily::Symbol,
    regime::BreadthRegime,
    C::Float64,
    target_r::Float64
)
    nb, basal_mask, consumers = consumers_and_basal()

    E = make_environment(rng, envkind)

    MODULE = nothing
    ranks = nothing
    prey = nothing

    if netfamily == :random
        prey = build_metaweb_random(rng, C, basal_mask)
    elseif netfamily == :modular
        prey, MODULE = build_metaweb_modular(rng, C, basal_mask)
    elseif netfamily == :heavytail
        prey = build_metaweb_heavytail(rng, C, basal_mask)
    elseif netfamily == :cascade
        prey, ranks = build_metaweb_cascade(rng, C, basal_mask)
    else
        error("Unknown netfamily: $netfamily")
    end

    σ = draw_sigmas(rng, regime)

    μ = Vector{Float64}(undef, S)
    for i in 1:nb
        μ[i] = rand(rng) * (E_MAX - E_MIN) + E_MIN
    end
    for i in (nb+1):S
        μ[i] = rand(rng) * (E_MAX - E_MIN) + E_MIN
    end

    achieved_r = assign_mus_with_target_corr!(rng, μ, prey, basal_mask, target_r)

    A_raw = Vector{BitVector}(undef, S)
    for i in 1:S
        A_raw[i] = suitability_mask_1d(E, μ[i], σ[i], SUIT_THRESH)
    end

    AB_raw = fixed_point_AB(A_raw, prey, basal_mask)

    A = Vector{BitVector}(undef, S)
    AB = Vector{BitVector}(undef, S)
    for i in 1:S
        A[i]  = apply_connectivity_filter(ws, A_raw[i], Emin_patch)
        AB[i] = apply_connectivity_filter(ws, AB_raw[i], Emin_patch)
    end

    SA  = gamma_richness_cons(A, basal_mask)
    SAB = gamma_richness_cons(AB, basal_mask)
    dSrel = (SA == 0) ? NaN : (1.0 - (SAB / SA))

    mjm = mean_jaccard_mismatch(A, AB, basal_mask)
    fa  = frac_affected(A, AB, basal_mask)
    ov  = realized_overlap(A_raw, prey, basal_mask)

    mvec = jaccard_mismatch_vec(A, AB, basal_mask)
    q90  = mismatch_q90(mvec)
    frac_gt = mismatch_frac_gt(mvec, TAIL_THRESH)

    L = count_links(prey, basal_mask)
    Creal = L / (S^2)

    return (
        dSrel=dSrel,
        mean_jaccard_mismatch=mjm,
        frac_affected=fa,
        realized_overlap=ov,
        achieved_r=achieved_r,
        Creal=Creal,
        mismatch_q90=q90,
        mismatch_frac_gt=frac_gt
    )
end

end