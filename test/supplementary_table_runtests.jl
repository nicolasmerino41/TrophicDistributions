using Test
include(joinpath(@__DIR__, "..", "SimulationsCode", "Functions.jl"))
include(joinpath(@__DIR__, "..", "sdm", "Parameters.jl"))
include(joinpath(@__DIR__, "..", "sdm", "SupplementaryTable.jl"))

@testset "Full SDM grid and supplementary aggregation" begin
    @test length(SDMParameters.DEGREES) == 15
    @test length(SDMParameters.CORRELATIONS) == 15
    @test SDMParameters.FIGURE_DEGREES == [2, 6]
    @test SDMParameters.FIGURE_CORRELATIONS == [0.0, 0.95]
    input = [
        (degree=1, target_r=0.0, environment="random", regime="A", n_consumers=1,
         mean_delta_auc=0.1, mean_delta_brier=0.02, mean_true_mismatch=0.4),
        (degree=1, target_r=0.0, environment="random", regime="A", n_consumers=9,
         mean_delta_auc=0.3, mean_delta_brier=NaN, mean_true_mismatch=0.6)
    ]
    rows = supplementary_rows(input; degrees=[1, 2], correlations=[0.0], attempted=2)
    @test length(rows) == 2
    @test rows[1].mean_delta_auc ≈ 0.2 # Equal community weights, not consumer weights.
    @test rows[1].se_delta_auc ≈ 0.1
    @test rows[1].n_auc == 2
    @test rows[1].n_brier == 1
    @test isnan(rows[1].se_delta_brier)
    @test rows[1].n_consumers == 10
    @test rows[2].n_auc == 0
    @test isnan(rows[2].mean_delta_auc)
    mktempdir() do output
        full_rows = write_supplementary_table(input, output)
        @test length(full_rows) == 225
        @test length(readlines(joinpath(output, "supplementary_sdm_continuum.tsv"))) == 226
        @test occursin("NA", read(joinpath(output, "supplementary_sdm_continuum.md"), String))
    end
end
