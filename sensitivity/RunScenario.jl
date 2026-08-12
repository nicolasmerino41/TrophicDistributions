length(ARGS) == 5 || error(
    "Usage: RunScenario.jl SCENARIO_ID VARIABLE VALUE LABEL REPLICATES"
)

scenario_id, variable, value_text, label, replicate_text = ARGS
replicates = parse(Int, replicate_text)
replicates > 0 || error("Replicates must be positive")

repository_dir = normpath(joinpath(@__DIR__, ".."))
include(joinpath(repository_dir, "SimulationsCode", "Functions.jl"))

using .Functions.Sweep: sweep_all
using .Functions.IO: save_table_tsv

results = sweep_all(replicates=replicates)
rows = [merge((
    scenario_id=scenario_id,
    sensitivity_variable=variable,
    sensitivity_value=parse(Float64, value_text),
    sensitivity_label=label
), row) for row in results.degree_summary]

output_dir = joinpath(@__DIR__, "Outputs", "scenarios")
mkpath(output_dir)
path = joinpath(output_dir, scenario_id * ".tsv")
save_table_tsv(path, rows)
println("Saved sensitivity scenario to: ", path)
