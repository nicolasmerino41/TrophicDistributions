module MechanisticCorrelation

using Random
using ..Parameters: S, E_MIN, E_MAX, RELAX_ITERS, MU_NOISE_SD
using Statistics

export pearson_r, prey_means, mechanistic_corr, assign_mus_with_target_corr!

function pearson_r(a::Vector{Float64}, b::Vector{Float64})
    ma = mean(a); mb = mean(b)
    sa = std(a); sb = std(b)
    (sa == 0 || sb == 0) && return 0.0
    return mean((a .- ma) .* (b .- mb)) / (sa * sb)
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

function mechanistic_corr(mu::Vector{Float64}, prey::Vector{Vector{Int}}, basal_mask::BitVector)
    pm = prey_means(mu, prey, basal_mask)
    idx = findall(i -> !basal_mask[i] && !isnan(pm[i]), 1:S)
    length(idx) < 6 && return 0.0
    return pearson_r(mu[idx], pm[idx])
end

function assign_mus_with_target_corr!(
    rng::AbstractRNG,
    mu::Vector{Float64},
    prey::Vector{Vector{Int}},
    basal_mask::BitVector,
    target_r::Float64;
    relax_iters::Int = RELAX_ITERS
)
    # Basal μ already set; this updates consumer μ using relaxation with α.
    consumers = findall(!, basal_mask)

    function relax(alpha::Float64)
        m = copy(mu)
        for _ in 1:relax_iters
            pm = prey_means(m, prey, basal_mask)
            @inbounds for i in consumers
                if isnan(pm[i]); continue; end
                m[i] = clamp((1-alpha)*m[i] + alpha*pm[i] + MU_NOISE_SD*randn(rng), E_MIN, E_MAX)
            end
        end
        return mechanistic_corr(m, prey, basal_mask), m
    end

    lo, hi = -0.98, 0.98
    best_err = Inf
    best_mu = copy(mu)
    best_r = 0.0

    for _ in 1:26
        mid = 0.5*(lo + hi)
        rmid, mmid = relax(mid)
        err = abs(rmid - target_r)
        if err < best_err
            best_err = err
            best_mu = mmid
            best_r = rmid
        end
        if rmid < target_r
            lo = mid
        else
            hi = mid
        end
    end

    copyto!(mu, best_mu)
    return best_r
end

end