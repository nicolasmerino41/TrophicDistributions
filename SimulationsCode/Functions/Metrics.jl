module Metrics

using ..Parameters: S
using Statistics

export gamma_richness_cons,
       mean_jaccard_mismatch,
       frac_affected,
       realized_overlap,
       jaccard_mismatch_vec,
       mismatch_q90,
       mismatch_frac_gt

function gamma_richness_cons(pres::Vector{BitVector}, basal_mask::BitVector)
    c = 0
    for i in 1:S
        basal_mask[i] && continue
        c += (count(pres[i]) > 0) ? 1 : 0
    end
    return c
end

function mean_jaccard_mismatch(A::Vector{BitVector}, AB::Vector{BitVector}, basal_mask::BitVector)
    vals = Float64[]
    for i in 1:S
        basal_mask[i] && continue
        Ai = A[i]
        count(Ai) == 0 && continue
        ABi = AB[i]
        inter = count(Ai .& ABi)
        uni   = count(Ai .| ABi)
        J = uni == 0 ? 1.0 : (inter / uni)
        push!(vals, 1 - J)
    end
    return isempty(vals) ? NaN : mean(vals)
end

function frac_affected(A::Vector{BitVector}, AB::Vector{BitVector}, basal_mask::BitVector)
    num = 0
    den = 0
    for i in 1:S
        basal_mask[i] && continue
        Ai = A[i]
        count(Ai) == 0 && continue
        den += 1
        num += (Ai != AB[i]) ? 1 : 0
    end
    return den == 0 ? NaN : num / den
end

function realized_overlap(A_raw::Vector{BitVector}, prey::Vector{Vector{Int}}, basal_mask::BitVector)
    vals = Float64[]
    for i in 1:S
        basal_mask[i] && continue
        Ai = A_raw[i]
        a = count(Ai)
        (a == 0 || isempty(prey[i])) && continue
        s = 0.0
        for j in prey[i]
            s += count(Ai .& A_raw[j]) / a
        end
        push!(vals, s / length(prey[i]))
    end
    return isempty(vals) ? NaN : mean(vals)
end

function jaccard_mismatch_vec(A::Vector{BitVector}, AB::Vector{BitVector}, basal_mask::BitVector)
    vals = Float64[]
    for i in 1:S
        basal_mask[i] && continue
        Ai = A[i]
        count(Ai) == 0 && continue
        ABi = AB[i]
        inter = count(Ai .& ABi)
        uni   = count(Ai .| ABi)
        J = (uni == 0) ? 1.0 : (inter / uni)
        push!(vals, 1.0 - J)
    end
    return vals
end

function q_from_sorted(v::Vector{Float64}, p::Float64)
    n = length(v)
    n == 0 && return NaN
    idx = clamp(round(Int, 1 + (n-1)*p), 1, n)
    return v[idx]
end

function mismatch_q90(vals::Vector{Float64})
    isempty(vals) && return NaN
    v = sort(vals)
    return q_from_sorted(v, 0.90)
end

function mismatch_frac_gt(vals::Vector{Float64}, thr::Float64)
    isempty(vals) && return NaN
    return mean(vals .> thr)
end

end