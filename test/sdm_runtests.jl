using Test
using Random
using Statistics

include(joinpath(@__DIR__, "..", "SimulationsCode", "Functions.jl"))
include(joinpath(@__DIR__, "..", "sdm", "Parameters.jl"))
include(joinpath(@__DIR__, "..", "sdm", "Functions.jl"))

using .Functions.Connectivity: make_workspaces
using .SDMParameters: DEGREES
using .SDMFunctions: fit_ridge_logistic, predict_probability,
                     build_jobs, run_world_sdms, summarize_communities,
                     summarize_results

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

@testset "Minimal framework-to-SDM handoff" begin
    job = only(build_jobs(
        environments=[:autocorr], regime_indices=[3],
        correlations=[0.475], replicates=1
    ))
    result = run_world_sdms(job, make_workspaces()[1])
    @test !isempty(result.rows)
    @test all(row -> row.degree in DEGREES, result.rows)
    @test all(row -> row.abiotic_converged && row.biotic_converged, result.rows)
    @test all(row -> isfinite(row.true_mismatch), result.rows)
    @test all(row -> isfinite(row.delta_auc), result.rows)
    @test all(row -> isfinite(row.delta_brier), result.rows)
    @test all(row -> 0.0 <= row.predictor_prevalence <= 1.0, result.rows)
    community_rows = summarize_communities(result.rows)
    @test length(community_rows) <= length(result.rows)
    @test all(row -> row.n_consumers >= 1, community_rows)
    summary = summarize_results(community_rows)
    @test !isempty(summary)
    @test all(row -> row.n >= 1, summary)
end

@testset "Correlation treatments use paired worlds" begin
    jobs = build_jobs(
        environments=[:random], regime_indices=[1],
        correlations=[0.0, 0.95], replicates=2
    )
    @test jobs[1].seed == jobs[3].seed
    @test jobs[2].seed == jobs[4].seed
    @test jobs[1].community_id != jobs[3].community_id
end
