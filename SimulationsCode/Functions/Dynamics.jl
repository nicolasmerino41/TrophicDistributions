module Dynamics

using ..Parameters: S, NCELLS

export trophic_distributions

function constrained_distribution(
    abiotic::BitVector,
    resource_distributions::Vector{BitVector},
    resources::Vector{Int}
)
    available = BitVector(falses(NCELLS))
    @inbounds for resource in resources
        available .|= resource_distributions[resource]
    end
    return abiotic .& available
end

"""
Return direct and recursive trophically constrained distributions.

The direct version applies the trophic constraint once using each prey species'
abiotic distribution. The recursive version starts from the abiotic maps and
reapplies the constraint synchronously until the community reaches a fixed
point, allowing restrictions to propagate through consumer-resource feedbacks.
"""
function trophic_distributions(
    A::Vector{BitVector},
    prey::Vector{Vector{Int}},
    basal_mask::BitVector
)
    direct = [copy(A[i]) for i in 1:S]
    for i in 1:S
        basal_mask[i] && continue
        direct[i] = constrained_distribution(A[i], A, prey[i])
    end

    recursive = copy(direct)
    next_state = [copy(A[i]) for i in 1:S]
    for iteration in 1:S
        changed = false
        for i in 1:S
            if basal_mask[i]
                copyto!(next_state[i], A[i])
            else
                next_state[i] = constrained_distribution(A[i], recursive, prey[i])
            end
            changed |= next_state[i] != recursive[i]
        end
        recursive, next_state = next_state, recursive
        !changed && return direct, recursive
    end
    error("Trophic fixed point did not converge within S=$S iterations")
end

end
