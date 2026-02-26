module Environment

using Random
using ..Parameters: NCELLS, AUTOCORR_ITERS, AUTOCORR_ALPHA, E_MIN, E_MAX
using ..Grid: NEIGH_4

export rescale_to_range!, smooth_field_once!, make_environment

function rescale_to_range!(v::Vector{Float64}, lo::Float64, hi::Float64)
    mn = minimum(v); mx = maximum(v)
    if mx == mn
        fill!(v, 0.5*(lo+hi))
        return v
    end
    @inbounds for i in eachindex(v)
        v[i] = lo + (hi-lo) * (v[i] - mn) / (mx - mn)
    end
    return v
end

function smooth_field_once!(E::Vector{Float64}, tmp::Vector{Float64}, α::Float64)
    @inbounds for i in 1:NCELLS
        s = E[i]
        nb = NEIGH_4[i]
        m = s
        for j in nb
            m += E[j]
        end
        m /= (1 + length(nb))
        tmp[i] = (1-α)*E[i] + α*m
    end
    copyto!(E, tmp)
    return E
end

function make_environment(rng::AbstractRNG, kind::Symbol)
    E = randn(rng, NCELLS)
    if kind == :autocorr
        tmp = similar(E)
        for _ in 1:AUTOCORR_ITERS
            smooth_field_once!(E, tmp, AUTOCORR_ALPHA)
        end
    elseif kind != :random
        error("Unknown env kind: $kind")
    end
    rescale_to_range!(E, E_MIN, E_MAX)
    return E
end

end