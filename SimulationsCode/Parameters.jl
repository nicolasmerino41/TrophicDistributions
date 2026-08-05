module Parameters

export NX, NY, NCELLS,
       S, BASAL_FRAC,
       USE_CONNECTIVITY_FILTER, Emin_patch,
       E_MIN, E_MAX,
       SUIT_THRESH,
       AUTOCORR_ITERS, AUTOCORR_ALPHA,
       DEGREE_CLASSES, CORR_RANGE, N_CORR,
       NREP,
       N_MODULES, MODULAR_IN_BIAS,
       CASCADE_LAMBDA,
       CORR_CALIBRATION_ITERS, CORR_CALIBRATION_RESTARTS,
       CORR_CALIBRATION_DAMPING, CORR_WITHIN_DEGREE_SD,
       TARGET_R_TOL, DEGREE_TARGET_R_TOL,
       BASE_SEED,
       OUTDIR

# ============================================================
# 1) Parameters
# ============================================================
# Spatial grid
const NX = 60
const NY = 60
const NCELLS = NX * NY

# Species pool
const S = 250
const BASAL_FRAC = 0.30  # basal species fraction

# Spatial viability filter (movement/connectivity proxy)
const USE_CONNECTIVITY_FILTER = true
const Emin_patch = 50  # LCC threshold for persistence

# Environmental field domain
const E_MIN = 0.0
const E_MAX = 100.0

# Niche suitability: Gaussian with threshold
# suitability = exp(-0.5 * ((E-μ)/σ)^2 ) >= SUIT_THRESH
const SUIT_THRESH = 0.25

# Environmental autocorrelation
const AUTOCORR_ITERS = 18
const AUTOCORR_ALPHA = 0.55 # 0..1 (higher = smoother)

# Focal-consumer degree design and community-level correlation sweep.
# Degree classes are assigned as evenly as possible within every community.
const DEGREE_CLASSES = [1, 2, 3, 4, 5, 6, 8, 10, 13, 17, 22, 28, 35, 44, 55]
const CORR_RANGE       = (0.0, 0.95)
const N_CORR    = 15

# Replicates per heatmap cell
const NREP = 50

const TAIL_THRESH = 0.8 # threshold for tail detection

# Network-family knobs
const N_MODULES = 6
const MODULAR_IN_BIAS = 6.0 # >1 increases within-module links vs between
const CASCADE_LAMBDA = 2.5 # 0 = uniform among lower ranks; higher = interval-like

# Degree-stratified niche-correlation calibration
const CORR_CALIBRATION_ITERS = 120
const CORR_CALIBRATION_RESTARTS = 3
const CORR_CALIBRATION_DAMPING = 0.70
# Preserve the expected SD of the original Uniform(E_MIN, E_MAX) optima.
# The calibrator scales it down only when required to remain inside the domain.
const CORR_WITHIN_DEGREE_SD = (E_MAX - E_MIN) / sqrt(12)
const TARGET_R_TOL = 0.03
const DEGREE_TARGET_R_TOL = 0.05

# Thread-safe seeds
const BASE_SEED = 20260202

# Output directory
OUTDIR = joinpath(@__DIR__, "..", "Outputs")
isdir(OUTDIR) || mkpath(OUTDIR)

end # module Parameters
