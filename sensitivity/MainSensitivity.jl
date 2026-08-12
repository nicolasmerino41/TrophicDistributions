include(joinpath(@__DIR__, "Parameters.jl"))
using .SensitivityParameters

mkpath(SCENARIO_DIR)
repository_dir = normpath(joinpath(@__DIR__, ".."))
runner = joinpath(@__DIR__, "RunScenario.jl")
replicates = isempty(ARGS) ? N_REPLICATES : parse(Int, only(ARGS))
replicates > 0 || error("Replicates must be positive")

for (index, scenario) in enumerate(SCENARIOS)
    println("Scenario $(index) / $(length(SCENARIOS)): $(scenario.label)")
    environment = copy(ENV)
    for key in (
        "TD_NX", "TD_NY", "TD_S", "TD_EMIN_PATCH",
        "TD_SUIT_THRESH", "TD_AUTOCORR_ITERS"
    )
        pop!(environment, key, nothing)
    end
    merge!(environment, scenario.overrides)

    command = `$(Base.julia_cmd()) --project=$repository_dir $runner $(scenario.id) $(scenario.variable) $(scenario.value) $(scenario.label) $replicates`
    run(setenv(command, environment))
end

println("All sensitivity scenarios completed.")
println("Run sensitivity/PlotSensitivity.R to create the combined figure.")
