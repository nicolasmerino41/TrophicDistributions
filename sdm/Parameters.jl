module SDMParameters

using ..Functions.Parameters: DEGREE_CLASSES, CORR_RANGE, N_CORR, NREP

export ENVIRONMENTS, NETWORKS, REGIME_INDICES, DEGREES, CORRELATIONS,
       N_REPLICATES, FOCAL_PRESENCE_SAMPLE_SIZE, BACKGROUND_SAMPLE_SIZE,
       MIN_FOCAL_PRESENCE_CELLS, SAMPLING_BIAS_STRENGTH,
       SAMPLING_BIAS_FLOOR, LOGISTIC_RIDGE, LOGISTIC_MAX_ITERS,
       LOGISTIC_TOL, BASE_SEED, PIPELINE_VERSION, OUTPUT_DIR,
       CHECKPOINT_DIR, FIGURE_DEGREES, FIGURE_CORRELATIONS

# The SDM experiment mirrors the framework grid. Environment, architecture,
# and niche-breadth regime are nuisance treatments over which Figure 5 averages.
const ENVIRONMENTS = [:random, :autocorr]
const NETWORKS = [:random, :modular, :cascade]
const REGIME_INDICES = collect(1:4)
const DEGREES = copy(DEGREE_CLASSES)
const CORRELATIONS = collect(range(CORR_RANGE[1], CORR_RANGE[2], length=N_CORR))
const N_REPLICATES = NREP

# Imperfect focal-species data.
const FOCAL_PRESENCE_SAMPLE_SIZE = 50
const BACKGROUND_SAMPLE_SIZE = 500
const MIN_FOCAL_PRESENCE_CELLS = 20
const SAMPLING_BIAS_STRENGTH = 7.0
const SAMPLING_BIAS_FLOOR = 0.05

# Ridge-stabilized presence-background logistic models.
const LOGISTIC_RIDGE = 0.1
const LOGISTIC_MAX_ITERS = 100
const LOGISTIC_TOL = 1e-8

# Four moderately contrasting framework cells highlighted in Figure 5.
# The correlation grid does not contain exactly 0.1 or 0.8, so use the nearest
# simulated treatments (approximately 0.068 and 0.814).
const FIGURE_DEGREES = [2, 6]
const FIGURE_CORRELATIONS = [
    CORRELATIONS[argmin(abs.(CORRELATIONS .- requested))]
    for requested in (0.1, 0.8)
]

const BASE_SEED = 20260805
const PIPELINE_VERSION = "figure5-simple-v1"
const OUTPUT_DIR = joinpath(@__DIR__, "Outputs")
const CHECKPOINT_DIR = joinpath(OUTPUT_DIR, "checkpoints", PIPELINE_VERSION)

end
