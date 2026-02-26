module Networks

using Random
using ..Parameters: S, BASAL_FRAC, N_MODULES, MODULAR_IN_BIAS,
                    HEAVYTAIL_GAMMA, HEAVYTAIL_KMAX_FRAC, CASCADE_LAMBDA
using Random

export consumers_and_basal,
       realized_connectance,
       ensure_min1_prey!,
       build_metaweb_random,
       build_metaweb_modular,
       build_metaweb_heavytail,
       build_metaweb_cascade

function consumers_and_basal()
    nb = round(Int, BASAL_FRAC * S)
    basal_mask = BitVector(falses(S))
    basal_mask[1:nb] .= true
    consumers = collect((nb+1):S)
    return nb, basal_mask, consumers
end

function realized_connectance(prey::Vector{Vector{Int}}, basal_mask::BitVector)
    L = 0
    for i in 1:S
        basal_mask[i] && continue
        L += length(prey[i])
    end
    return L / (S^2)
end

function ensure_min1_prey!(rng::AbstractRNG, prey::Vector{Vector{Int}}, basal_mask::BitVector)
    # Guarantee each consumer has at least one prey.
    for i in 1:S
        basal_mask[i] && continue
        if isempty(prey[i])
            candidates = findall(basal_mask)
            if isempty(candidates)
                explained = [j for j in 1:S if j != i]
                push!(prey[i], explained[rand(rng, 1:length(explained))])
            else
                push!(prey[i], candidates[rand(rng, 1:length(candidates))])
            end
        end
    end
end

function build_metaweb_random(rng::AbstractRNG, C::Float64, basal_mask::BitVector)
    # Random consumer->prey edges; exact Ltarget achieved (after ensuring 1 prey each).
    prey = [Int[] for _ in 1:S]
    consumers = findall(!, basal_mask)
    Ltarget = round(Int, C * S^2)
    # Step 1: 1 prey each
    for i in consumers
        cand = [j for j in 1:S if j != i]
        push!(prey[i], cand[rand(rng, 1:length(cand))])
    end
    # Step 2: fill remaining edges uniformly without duplicates
    L = length(consumers)
    while L < Ltarget
        i = consumers[rand(rng, 1:length(consumers))]
        cand = rand(rng, 1:S)
        if cand == i; continue; end
        if cand ∉ prey[i]
            push!(prey[i], cand)
            L += 1
        end
    end
    return prey
end

function build_metaweb_modular(rng::AbstractRNG, C::Float64, basal_mask::BitVector)
    # Modular SBM-like: higher probability within module.
    prey = [Int[] for _ in 1:S]
    consumers = findall(!, basal_mask)
    # module assignment for all species
    MODULE = Vector{Int}(undef, S)
    for i in 1:S
        MODULE[i] = 1 + (i - 1) % N_MODULES
    end
    Ltarget = round(Int, C * S^2)

    function sample_prey(i::Int)
        inmod = Int[]
        outmod = Int[]
        mi = MODULE[i]
        for j in 1:S
            if j == i; continue; end
            if MODULE[j] == mi
                push!(inmod, j)
            else
                push!(outmod, j)
            end
        end
        # weight within vs outside
        if isempty(inmod)
            return outmod[rand(rng, 1:length(outmod))]
        elseif isempty(outmod)
            return inmod[rand(rng, 1:length(inmod))]
        else
            if rand(rng) < MODULAR_IN_BIAS / (MODULAR_IN_BIAS + 1.0)
                return inmod[rand(rng, 1:length(inmod))]
            else
                return outmod[rand(rng, 1:length(outmod))]
            end
        end
    end

    # Step 1: 1 prey each
    for i in consumers
        push!(prey[i], sample_prey(i))
    end
    # Step 2: fill remaining
    L = length(consumers)
    while L < Ltarget
        i = consumers[rand(rng, 1:length(consumers))]
        j = sample_prey(i)
        if j ∉ prey[i]
            push!(prey[i], j)
            L += 1
        end
    end
    return prey, MODULE
end

function build_metaweb_heavytail(rng::AbstractRNG, C::Float64, basal_mask::BitVector)
    # Heavy-tailed consumer out-degree, overall Ltarget.
    prey = [Int[] for _ in 1:S]
    consumers = findall(!, basal_mask)
    nC = length(consumers)
    Ltarget = round(Int, C * S^2)

    # weights ~ Pareto-like
    w = zeros(Float64, nC)
    for k in 1:nC
        # Pareto with exponent gamma: inverse-CDF with u^( -1/(gamma-1) )
        u = rand(rng)
        w[k] = u^(-1/(HEAVYTAIL_GAMMA-1))
    end
    w ./= sum(w)

    deg = ones(Int, nC)
    remaining = max(0, Ltarget - nC)
    
    for _ in 1:remaining
        r = rand(rng)
        acc = 0.0
        idx = 1
        for k in 1:nC
            acc += w[k]
            if r <= acc
                idx = k
                break
            end
        end
        deg[idx] += 1
    end

    # cap degrees and reassign overflow to keep total ~ Ltarget
    kmax = max(2, round(Int, HEAVYTAIL_KMAX_FRAC * (S-1)))
    overflow = 0
    for k in 1:nC
        if deg[k] > kmax
            overflow += deg[k] - kmax
            deg[k] = kmax
        end
    end
    # redistribute overflow
    for _ in 1:overflow
        k = rand(rng, 1:nC)
        if deg[k] < kmax
            deg[k] += 1
        end
    end

    # choose prey sets
    for (kk, i) in enumerate(consumers)
        candidates = [j for j in 1:S if j != i]
        shuffle!(rng, candidates)
        d = min(deg[kk], length(candidates))
        prey[i] = candidates[1:d]
    end
    ensure_min1_prey!(rng, prey, basal_mask)
    return prey
end

function build_metaweb_cascade(rng::AbstractRNG, C::Float64, basal_mask::BitVector)
    # Cascade hierarchy: consumers feed on lower-ranked species.
    prey = [Int[] for _ in 1:S]
    consumers = findall(!, basal_mask)
    nC = length(consumers)
    Ltarget = round(Int, C * S^2)

    # ranks for all species; enforce every consumer has at least one lower-ranked prey candidate
    ranks = zeros(Float64, S)
    for attempt in 1:200
        for i in 1:S
            ranks[i] = rand(rng)
        end
        ok = true
        for i in consumers
            lower = findall(j -> ranks[j] < ranks[i] && j != i, 1:S)
            if isempty(lower)
                ok = false
                break
            end
        end
        ok && break
        attempt == 200 && error("Failed to sample cascade ranks with valid lower-prey candidates for all consumers")
    end

    # helper: sample prey among lower ranks with exponential bias toward nearby ranks
    function sample_lower_prey(i::Int)
        lower = Int[]
        w = Float64[]
        ri = ranks[i]
        for j in 1:S
            if j == i; continue; end
            if ranks[j] < ri
                push!(lower, j)
                # weight: exp(-λ * (ri-rj))
                push!(w, exp(-CASCADE_LAMBDA * (ri - ranks[j])))
            end
        end
        # normalize and sample
        sw = sum(w)
        r = rand(rng) * sw
        acc = 0.0
        for k in 1:length(lower)
            acc += w[k]
            if r <= acc
                return lower[k]
            end
        end
        return lower[end]
    end

    # Step 1: 1 prey each
    for i in consumers
        push!(prey[i], sample_lower_prey(i))
    end
    # Step 2: fill remaining edges
    L = nC
    while L < Ltarget
        i = consumers[rand(rng, 1:nC)]
        j = sample_lower_prey(i)
        if j ∉ prey[i]
            push!(prey[i], j)
            L += 1
        end
    end

    return prey, ranks
end

end