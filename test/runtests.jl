using Test
using Random

include(joinpath(@__DIR__, "..", "SimulationsCode", "Functions.jl"))

using .Functions.Parameters: S, NCELLS, E_MIN, E_MAX, DEGREE_CLASSES,
                             TARGET_R_TOL, DEGREE_TARGET_R_TOL
using .Functions.Networks: consumers_and_basal, assign_balanced_degrees,
                           build_trophic_community, validate_metaweb
using .Functions.Dynamics: trophic_distributions
using .Functions.Metrics: jaccard_mismatch_by_species
using .Functions.MechanisticCorrelation: pearson_r, assign_mus_with_target_corr!,
                                         degree_corr_diagnostics,
                                         max_degree_target_error
using .Functions.Connectivity: make_workspaces
using .Functions.Niches: regimes
using .Functions.Simulation: simulate_one!

@testset "Balanced neutral trophic community" begin
    rng = MersenneTwister(41)
    _, basal_mask, consumers = consumers_and_basal()
    target_degree = assign_balanced_degrees(rng, basal_mask)
    prey = build_trophic_community(rng, basal_mask, target_degree)

    counts = [count(==(degree), target_degree[consumers]) for degree in DEGREE_CLASSES]
    @test maximum(counts) - minimum(counts) <= 1
    @test validate_metaweb(prey, basal_mask, target_degree)
    @test length.(prey)[consumers] == target_degree[consumers]
    @test all(i -> !(i in prey[i]), consumers)
    @test all(i -> length(prey[i]) == length(unique(prey[i])), consumers)
    @test any(i -> any(j -> j in consumers, prey[i]), consumers)
end

@testset "Degree-stratified niche correlation calibration" begin
    rng = MersenneTwister(31415)
    _, basal_mask, _ = consumers_and_basal()
    target_degree = assign_balanced_degrees(rng, basal_mask)
    prey = build_trophic_community(rng, basal_mask, target_degree)

    for target in (0.0, 0.5, 0.75, 0.95)
        mu = [rand(rng) * (E_MAX - E_MIN) + E_MIN for _ in 1:S]
        achieved = assign_mus_with_target_corr!(rng, mu, prey, basal_mask, target)
        diagnostics = degree_corr_diagnostics(mu, prey, basal_mask)
        @test abs(achieved - target) <= TARGET_R_TOL
        @test max_degree_target_error(diagnostics, target) <= DEGREE_TARGET_R_TOL
    end
end

@testset "Recursive constraints propagate through consumer resources" begin
    A = [BitVector(falses(NCELLS)) for _ in 1:S]
    prey = [Int[] for _ in 1:S]
    basal_mask = BitVector(falses(S))
    basal_mask[1] = true

    intermediate = 2
    top_consumer = 3
    prey[intermediate] = [1]
    prey[top_consumer] = [intermediate]
    A[1][1:5] .= true
    A[intermediate][1:10] .= true
    A[top_consumer][1:10] .= true

    direct, recursive = trophic_distributions(A, prey, basal_mask)
    @test count(direct[top_consumer]) == 10
    @test count(recursive[top_consumer]) == 5
    @test all(recursive[top_consumer] .<= direct[top_consumer])
end

@testset "Mismatch remains aligned with species identity" begin
    _, basal_mask, consumers = consumers_and_basal()
    A = [BitVector(falses(NCELLS)) for _ in 1:S]
    AB = [BitVector(falses(NCELLS)) for _ in 1:S]

    focal = consumers[1]
    A[focal][1:10] .= true
    mismatch, eligible = jaccard_mismatch_by_species(A, AB, basal_mask)

    @test eligible[focal]
    @test mismatch[focal] == 1.0
    @test !eligible[consumers[2]]
    @test isnan(mismatch[consumers[2]])
end

@testset "Pearson correlation uses the exact sample definition" begin
    values = collect(1.0:11.0)
    @test pearson_r(values, values) ≈ 1.0 atol=1e-12
    @test pearson_r(values, reverse(values)) ≈ -1.0 atol=1e-12
end

@testset "Reduced one-community smoke simulation" begin
    rng = MersenneTwister(2026)
    output = simulate_one!(
        rng, make_workspaces()[1], :random, regimes[1], 0.5, 1
    )
    _, _, consumers = consumers_and_basal()
    @test length(output.consumers) == length(consumers)
    @test all(row -> row.target_degree == row.realized_degree, output.consumers)
    @test all(row -> !row.eligible || row.indirect_increment >= -1e-12, output.consumers)
    @test output.community.community_id == 1
    @test output.community.correlation_ok
    @test Set(row.target_degree for row in output.consumers) == Set(DEGREE_CLASSES)
    @test Set(keys(output.degree_correlations)) == Set(DEGREE_CLASSES)
end
