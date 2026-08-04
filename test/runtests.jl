using Test
using Random

include(joinpath(@__DIR__, "..", "SimulationsCode", "Functions.jl"))

using .Functions.Parameters: S, NCELLS, E_MIN, E_MAX, DEGREE_CLASSES,
                             TARGET_R_TOL, DEGREE_TARGET_R_TOL
using .Functions.Networks: consumers_and_basal, assign_balanced_degrees,
                           build_metaweb_random, build_metaweb_modular,
                           build_metaweb_cascade, validate_metaweb
using .Functions.Metrics: jaccard_mismatch_by_species
using .Functions.MechanisticCorrelation: pearson_r, assign_mus_with_target_corr!,
                                         degree_corr_diagnostics,
                                         max_degree_target_error
using .Functions.Connectivity: make_workspaces
using .Functions.Niches: regimes
using .Functions.Simulation: simulate_one!

@testset "Balanced degree-controlled networks" begin
    rng = MersenneTwister(41)
    _, basal_mask, consumers = consumers_and_basal()
    target_degree = assign_balanced_degrees(rng, basal_mask)

    counts = [count(==(degree), target_degree[consumers]) for degree in DEGREE_CLASSES]
    @test maximum(counts) - minimum(counts) <= 1

    random_prey = build_metaweb_random(rng, basal_mask, target_degree)
    modular_prey, _ = build_metaweb_modular(rng, basal_mask, target_degree)
    cascade_prey, ranks = build_metaweb_cascade(rng, basal_mask, target_degree)

    for prey in (random_prey, modular_prey, cascade_prey)
        @test validate_metaweb(prey, basal_mask, target_degree)
        @test length.(prey)[consumers] == target_degree[consumers]
        @test all(i -> !(i in prey[i]), consumers)
        @test all(i -> length(prey[i]) == length(unique(prey[i])), consumers)
    end
    @test all(i -> all(ranks[j] < ranks[i] for j in cascade_prey[i]), consumers)
end

@testset "Degree-stratified niche correlation calibration" begin
    rng = MersenneTwister(31415)
    _, basal_mask, _ = consumers_and_basal()
    target_degree = assign_balanced_degrees(rng, basal_mask)
    random_prey = build_metaweb_random(rng, basal_mask, target_degree)
    modular_prey, _ = build_metaweb_modular(rng, basal_mask, target_degree)
    cascade_prey, _ = build_metaweb_cascade(rng, basal_mask, target_degree)

    for prey in (random_prey, modular_prey, cascade_prey),
        target in (0.0, 0.5, 0.75, 1.0)
        mu = [rand(rng) * (E_MAX - E_MIN) + E_MIN for _ in 1:S]
        achieved = assign_mus_with_target_corr!(rng, mu, prey, basal_mask, target)
        diagnostics = degree_corr_diagnostics(mu, prey, basal_mask)
        @test abs(achieved - target) <= TARGET_R_TOL
        @test max_degree_target_error(diagnostics, target) <= DEGREE_TARGET_R_TOL
    end
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
    workspace = make_workspaces()[1]
    output = simulate_one!(rng, workspace, :random, :random, regimes[1], 0.5, 1)
    _, _, consumers = consumers_and_basal()
    @test length(output.consumers) == length(consumers)
    @test all(row -> row.target_degree == row.realized_degree, output.consumers)
    @test output.community.community_id == 1
    @test output.community.overall_correlation_ok
    @test output.community.degree_correlation_ok
    @test output.community.correlation_ok
    @test Set(row.target_degree for row in output.consumers) == Set(DEGREE_CLASSES)
    @test Set(keys(output.degree_correlations)) == Set(DEGREE_CLASSES)
end
