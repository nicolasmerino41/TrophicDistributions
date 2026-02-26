module Parameters

export NX, NY, NCELLS,
       S, BASAL_FRAC,
       USE_CONNECTIVITY_FILTER, Emin_patch,
       E_MIN, E_MAX,
       SUIT_THRESH,
       AUTOCORR_ITERS, AUTOCORR_ALPHA,
       CONNECTANCE_RANGE, CORR_RANGE, N_CONNECT, N_CORR,
       NREP,
       N_MODULES, MODULAR_IN_BIAS,
       HEAVYTAIL_GAMMA, HEAVYTAIL_KMAX_FRAC,
       CASCADE_LAMBDA,
       RELAX_ITERS, MU_NOISE_SD, TARGET_R_TOL,
       BASE_SEED,
       OUTDIR

# ============================================================
# 1) Parameters (Have fun with these!)
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

# Sweep axes
const CONNECTANCE_RANGE = (0.005, 0.15)
const CORR_RANGE       = (0.0, 1.0)
const N_CONNECT = 2
const N_CORR    = 2

# Replicates per heatmap cell
const NREP = 2

const TAIL_THRESH = 0.8 # threshold for tail detection

# Network-family knobs
const N_MODULES = 6
const MODULAR_IN_BIAS = 6.0 # >1 increases within-module links vs between
const HEAVYTAIL_GAMMA = 2.2 # out-degree heaviness
const HEAVYTAIL_KMAX_FRAC = 0.35 # kmax = round(frac*(S-1))
const CASCADE_LAMBDA = 2.5 # 0 = uniform among lower ranks; higher = interval-like

# Mechanistic niche-correlation
const RELAX_ITERS = 30
const MU_NOISE_SD = 1.8            # consumer μ jitter during relaxation
const TARGET_R_TOL = 0.03

# Thread-safe seeds
const BASE_SEED = 20260202

# Output directory
OUTDIR = joinpath(@__DIR__, "..", "Outputs")
isdir(OUTDIR) || mkpath(OUTDIR)

end # module Parameters