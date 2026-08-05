module SDMParameters

using ..Functions.Parameters: DEGREE_CLASSES, CORR_RANGE, N_CORR, NREP

export ENVIRONMENTS, NETWORKS, REGIME_INDICES, DEGREES, CORRELATIONS,
       N_REPLICATES, LINK_COMPLETENESS_LEVELS, FOCALS_PER_DEGREE,
       FOCAL_PRESENCE_SAMPLE_SIZE, PREY_PRESENCE_SAMPLE_SIZE,
       BACKGROUND_SAMPLE_SIZE, MIN_FOCAL_PRESENCE_CELLS,
       MIN_PREY_PRESENCE_CELLS, SAMPLING_BIAS_STRENGTH,
       SAMPLING_BIAS_FLOOR, BIASED_BACKGROUND,
       LOGISTIC_RIDGE, LOGISTIC_MAX_ITERS, LOGISTIC_TOL,
       BASE_SEED, PIPELINE_VERSION, OUTPUT_DIR, CHECKPOINT_DIR,
       MISMATCH_BIN_WIDTH

# Complete Figure 5 design. These treatments deliberately mirror the main
# simulation sweep, so every SDM result maps back to a degree-by-correlation
# heatmap cell.
const ENVIRONMENTS = [:random, :autocorr]
const NETWORKS = [:random, :modular, :cascade]
const REGIME_INDICES = collect(1:4)
const DEGREES = copy(DEGREE_CLASSES)
const CORRELATIONS = collect(range(CORR_RANGE[1], CORR_RANGE[2], length=N_CORR))
const N_REPLICATES = NREP
const FOCALS_PER_DEGREE = 1

# Biotic-information stages. A value of 1.0 means that every true trophic link
# is known; lower values retain a nested random subset of links. The oracle
# stage is always evaluated separately using the true prey distributions.
const LINK_COMPLETENESS_LEVELS = [1.0, 0.75, 0.50, 0.25]

# Imperfect occurrence data for focal consumers and their prey.
const FOCAL_PRESENCE_SAMPLE_SIZE = 50
const PREY_PRESENCE_SAMPLE_SIZE = 30
const BACKGROUND_SAMPLE_SIZE = 500
const MIN_FOCAL_PRESENCE_CELLS = 20
const MIN_PREY_PRESENCE_CELLS = 10

# Spatially biased sampling and target-group background.
const SAMPLING_BIAS_STRENGTH = 7.0
const SAMPLING_BIAS_FLOOR = 0.05
const BIASED_BACKGROUND = true

# Ridge-stabilized presence-background logistic models.
const LOGISTIC_RIDGE = 0.1
const LOGISTIC_MAX_ITERS = 100
const LOGISTIC_TOL = 1e-8

const BASE_SEED = 20260805
const PIPELINE_VERSION = "figure5-sdm-v1"
const OUTPUT_DIR = joinpath(@__DIR__, "Outputs")
const CHECKPOINT_DIR = joinpath(OUTPUT_DIR, "checkpoints", PIPELINE_VERSION)
const MISMATCH_BIN_WIDTH = 0.05

end
