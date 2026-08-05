using Test
using Random
using Statistics

include(joinpath(@__DIR__, "..", "SimulationsCode", "Functions.jl"))
include(joinpath(@__DIR__, "..", "sdm", "Parameters.jl"))
include(joinpath(@__DIR__, "..", "sdm", "Functions.jl"))

using .Functions.Parameters: NCELLS
using .Functions.Connectivity: make_workspaces
using .Functions.Niches: regimes
using .Functions.Simulation: simulate_world!, summarize_world
using .SDMParameters: DEGREES, LINK_COMPLETENESS_LEVELS
using .SDMFunctions: fit_ridge_logistic, predict_probability,
                     estimated_resource_availability, build_jobs, run_world_sdms

@testset "Ridge logistic SDM solver" begin
    predictor = collect(range(-3.0, 3.0, length=200))
    design = hcat(ones(length(predictor)), predictor)
    response = Float64.(predictor .> 0)
    model = fit_ridge_logistic(design, response)
    probability = predict_probability(model, design)
    @test model.converged
    @test all(0.0 .<= probability .<= 1.0)
    @test mean(probability[predictor .> 1]) > 0.95
    @test mean(probability[predictor .< -1]) < 0.05
end

@testset "Estimated resource combination" begin
    predictions = Dict(
        1 => fill(0.2, NCELLS),
        2 => fill(0.5, NCELLS)
    )
    combined = estimated_resource_availability(predictions, [1, 2])
    @test all(isapprox.(combined, 0.6; atol=1e-12))
    @test all(estimated_resource_availability(predictions, Int[]) .== 0.0)
end

@testset "Complete staged SDM handoff" begin
    jobs = build_jobs(
        environments=[:autocorr], networks=[:random], regime_indices=[3],
        correlations=[0.475], replicates=1
    )
    result = run_world_sdms(only(jobs), make_workspaces()[1])
    @test length(result.community_rows) == 1
    @test result.community_rows[1].correlation_ok
    @test !isempty(result.focal_rows)
    @test !isempty(result.comparison_rows)
    @test all(row -> isfinite(row.true_mismatch), result.comparison_rows)
    @test all(row -> row.abiotic_converged && row.biotic_converged,
              result.comparison_rows)
    expected_stages = Set(vcat("oracle", [
        "estimated_" * string(round(Int, 100 * level))
        for level in LINK_COMPLETENESS_LEVELS
    ]))
    @test Set(getproperty.(result.comparison_rows, :information_level)) == expected_stages
    @test all(row -> row.known_links <= row.true_resource_links,
              result.comparison_rows)
    @test all(row -> row.modeled_links <= row.known_links,
              result.comparison_rows)
end
