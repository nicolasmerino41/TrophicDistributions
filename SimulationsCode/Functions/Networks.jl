module Networks

using Random
using ..Parameters: S, BASAL_FRAC, DEGREE_CLASSES

export consumers_and_basal,
       assign_balanced_degrees,
       realized_degrees,
       validate_metaweb,
       build_trophic_community

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

"""
Build one neutral trophic community. Every consumer draws exactly its assigned
number of prey uniformly without replacement from the complete species pool,
excluding itself. Consumers can be resources and reciprocal or longer trophic
feedbacks can occur, but no architecture treatment is imposed.
"""
function build_trophic_community(
    rng::AbstractRNG,
    basal_mask::BitVector,
    target_degree::Vector{Int}
)
    prey = [Int[] for _ in 1:S]
    consumers = findall(!, basal_mask)

    for i in consumers
        candidates = [j for j in 1:S if j != i]
        shuffle!(rng, candidates)
        prey[i] = candidates[1:target_degree[i]]
    end
    validate_metaweb(prey, basal_mask, target_degree)
    return prey
end

end
