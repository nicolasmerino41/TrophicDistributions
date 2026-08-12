module SensitivityParameters

export N_REPLICATES, BASELINE, SCENARIOS, OUTPUT_DIR, SCENARIO_DIR

const N_REPLICATES = 50

const BASELINE = (
    id="baseline",
    variable="Baseline",
    value=0.0,
    label="Baseline",
    overrides=Dict{String,String}()
)

function grid_scenario(size::Int)
    # Diffusive smoothing distance scales with sqrt(iterations). Scaling the
    # iteration count by size^2 preserves its approximate landscape-relative
    # spatial scale when resolution changes.
    autocorrelation_iterations = max(1, round(Int, 18 * (size / 60)^2))
    return (
        id="grid_$(size)", variable="Grid dimension", value=Float64(size),
        label="$(size) x $(size)",
        overrides=Dict(
            "TD_NX" => string(size),
            "TD_NY" => string(size),
            "TD_AUTOCORR_ITERS" => string(autocorrelation_iterations),
            "TD_EMIN_PATCH" => string(round(Int, (50 / 3600) * size^2))
        )
    )
end

function patch_scenario(fraction::Float64)
    cells = round(Int, fraction * 60 * 60)
    return (
        id="patch_" * replace(string(fraction), "." => "p"),
        variable="Minimum patch", value=fraction,
        label=string(round(100 * fraction, digits=2), "%"),
        overrides=Dict("TD_EMIN_PATCH" => string(cells))
    )
end

function threshold_scenario(value::Float64)
    return (
        id="threshold_" * replace(string(value), "." => "p"),
        variable="Suitability threshold", value=value,
        label=string(value),
        overrides=Dict("TD_SUIT_THRESH" => string(value))
    )
end

function richness_scenario(value::Int)
    return (
        id="richness_$(value)", variable="Species richness", value=Float64(value),
        label=string(value), overrides=Dict("TD_S" => string(value))
    )
end

const SCENARIOS = vcat(
    [BASELINE],
    grid_scenario.([20, 40, 80, 100]),
    patch_scenario.([0.0, 0.005, 0.01, 0.015, 0.02, 0.03]),
    threshold_scenario.([0.15, 0.20, 0.30, 0.35]),
    richness_scenario.([150, 200, 300, 350, 400])
)

const OUTPUT_DIR = joinpath(@__DIR__, "Outputs")
const SCENARIO_DIR = joinpath(OUTPUT_DIR, "scenarios")

end
