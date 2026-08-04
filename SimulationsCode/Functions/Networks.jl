module Networks

using Random
using ..Parameters: S, BASAL_FRAC, DEGREE_CLASSES, N_MODULES,
                    MODULAR_IN_BIAS, CASCADE_LAMBDA

export consumers_and_basal,
       assign_balanced_degrees,
       realized_degrees,
       realized_connectance,
       validate_metaweb,
       build_metaweb_random,
       build_metaweb_modular,
       build_metaweb_cascade

function consumers_and_basal()
    nb = round(Int, BASAL_FRAC * S)
    basal_mask = BitVector(falses(S))
    basal_mask[1:nb] .= true
    consumers = collect((nb + 1):S)
    return nb, basal_mask, consumers
end

"""Assign the requested degree classes as evenly as possible, then shuffle them."""
function assign_balanced_degrees(
    rng::AbstractRNG,
    basal_mask::BitVector;
    degree_classes::Vector{Int} = DEGREE_CLASSES
)
    isempty(degree_classes) && error("DEGREE_CLASSES must not be empty")
    any(d -> d < 1 || d >= S, degree_classes) &&
        error("Every degree class must be between 1 and S-1")

    consumers = findall(!, basal_mask)
    assigned = [degree_classes[mod1(i, length(degree_classes))] for i in eachindex(consumers)]
    shuffle!(rng, assigned)

    target_degree = zeros(Int, S)
    for (i, consumer) in enumerate(consumers)
        target_degree[consumer] = assigned[i]
    end
    return target_degree
end

realized_degrees(prey::Vector{Vector{Int}}) = length.(prey)

function realized_connectance(prey::Vector{Vector{Int}}, basal_mask::BitVector)
    links = sum(length(prey[i]) for i in eachindex(prey) if !basal_mask[i])
    return links / (S^2)
end

function validate_metaweb(
    prey::Vector{Vector{Int}},
    basal_mask::BitVector,
    target_degree::Vector{Int}
)
    length(prey) == S || error("Expected prey lists for S=$S species")
    length(target_degree) == S || error("Expected S=$S target degrees")

    for i in 1:S
        basal_mask[i] && !isempty(prey[i]) && error("Basal species $i has prey")
        !basal_mask[i] && length(prey[i]) != target_degree[i] &&
            error("Consumer $i has degree $(length(prey[i])); expected $(target_degree[i])")
        i in prey[i] && error("Self-link detected for species $i")
        length(unique(prey[i])) == length(prey[i]) ||
            error("Duplicate prey detected for consumer $i")
    end
    return true
end

function weighted_sample_without_replacement(
    rng::AbstractRNG,
    candidates::Vector{Int},
    weights::Vector{Float64},
    n::Int
)
    n <= length(candidates) || error("Requested degree $n exceeds available prey")
    pool = copy(candidates)
    w = copy(weights)
    selected = Int[]
    sizehint!(selected, n)

    for _ in 1:n
        total = sum(w)
        total > 0 || error("Prey-selection weights must have positive mass")
        draw = rand(rng) * total
        cumulative = 0.0
        chosen = length(pool)
        for idx in eachindex(pool)
            cumulative += w[idx]
            if draw <= cumulative
                chosen = idx
                break
            end
        end
        push!(selected, pool[chosen])
        deleteat!(pool, chosen)
        deleteat!(w, chosen)
    end
    return selected
end

function build_metaweb_random(
    rng::AbstractRNG,
    basal_mask::BitVector,
    target_degree::Vector{Int}
)
    prey = [Int[] for _ in 1:S]
    for i in findall(!, basal_mask)
        candidates = [j for j in 1:S if j != i]
        shuffle!(rng, candidates)
        prey[i] = candidates[1:target_degree[i]]
    end
    validate_metaweb(prey, basal_mask, target_degree)
    return prey
end

function build_metaweb_modular(
    rng::AbstractRNG,
    basal_mask::BitVector,
    target_degree::Vector{Int}
)
    prey = [Int[] for _ in 1:S]
    modules = [mod1(i, N_MODULES) for i in 1:S]

    for i in findall(!, basal_mask)
        candidates = [j for j in 1:S if j != i]
        weights = [modules[j] == modules[i] ? MODULAR_IN_BIAS : 1.0 for j in candidates]
        prey[i] = weighted_sample_without_replacement(
            rng, candidates, weights, target_degree[i]
        )
    end
    validate_metaweb(prey, basal_mask, target_degree)
    return prey, modules
end

function build_metaweb_cascade(
    rng::AbstractRNG,
    basal_mask::BitVector,
    target_degree::Vector{Int}
)
    prey = [Int[] for _ in 1:S]
    basal = findall(identity, basal_mask)
    consumers = findall(!, basal_mask)

    # All basal species are below consumers. Consumer order is randomized, so degree
    # remains independent of position while every consumer has at least |basal| prey.
    consumer_order = shuffle(rng, copy(consumers))
    ranks = zeros(Int, S)
    for (rank, species) in enumerate(basal)
        ranks[species] = rank
    end
    for (offset, species) in enumerate(consumer_order)
        ranks[species] = length(basal) + offset
    end

    for i in consumers
        candidates = [j for j in 1:S if ranks[j] < ranks[i]]
        weights = [exp(-CASCADE_LAMBDA * (ranks[i] - ranks[j]) / S) for j in candidates]
        prey[i] = weighted_sample_without_replacement(
            rng, candidates, weights, target_degree[i]
        )
        all(ranks[j] < ranks[i] for j in prey[i]) ||
            error("Cascade constraint failed for consumer $i")
    end
    validate_metaweb(prey, basal_mask, target_degree)
    return prey, ranks
end

end
