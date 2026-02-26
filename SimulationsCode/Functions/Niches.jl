module Niches

using Random
using ..Parameters: NCELLS, S, SUIT_THRESH, E_MIN, E_MAX
using Statistics

export suitability_mask_1d, BreadthRegime, regimes, draw_sigmas

@inline function suitability_mask_1d(E::Vector{Float64}, μ::Float64, σ::Float64, thresh::Float64)
    lim = sqrt(-2.0 * log(thresh))  # |(E-μ)/σ| <= lim
    invσ = 1.0 / max(σ, 1e-6)
    m = BitVector(undef, NCELLS)
    @inbounds for i in 1:NCELLS
        z = (E[i] - μ) * invσ
        m[i] = abs(z) <= lim
    end
    return m
end

struct BreadthRegime
    name::String
    meanσ::Float64
    logsd::Float64
end

const regimes = [
    BreadthRegime("Narrow + LowVar",  7.5, 0.20),
    BreadthRegime("Narrow + HighVar", 7.5, 0.55),
    BreadthRegime("Broad + LowVar",   16.0, 0.20),
    BreadthRegime("Broad + HighVar",  16.0, 0.55),
]

function draw_sigmas(rng::AbstractRNG, regime::BreadthRegime)
    σ = Vector{Float64}(undef, S)
    # lognormal around meanσ with logsd; clamp to reasonable bounds
    for i in 1:S
        val = exp(log(regime.meanσ) + regime.logsd * randn(rng))
        σ[i] = clamp(val, 1.5, 45.0)
    end
    return σ
end

end