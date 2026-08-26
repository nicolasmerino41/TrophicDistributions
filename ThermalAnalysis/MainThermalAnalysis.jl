using CSV, DataFrames, Statistics, Random

const SCRIPT_DIR = @__DIR__
const ROOT_DIR = normpath(joinpath(SCRIPT_DIR, ".."))
const DATA_DIR = joinpath(ROOT_DIR, "Data")
const OUTROOT = joinpath(ROOT_DIR, "Outputs", "thermal_metrics")
const MASTER_FILE = joinpath(DATA_DIR, "thermal_interactions_master.csv")
const COVERAGE_FILE = joinpath(DATA_DIR, "thermal_consumer_coverage.csv")

const TRAITS = ["ctmax", "lt50", "ctmin", "ltmax"]
const N_PERM = 1000
const N_SWEEPS = 10
const SEED = 1

isdir(OUTROOT) || mkpath(OUTROOT)
isfile(MASTER_FILE) || error(
    "Missing Data/thermal_interactions_master.csv. Run ThermalAnalysis/BuildThermalDataset.R first."
)

function edge_swap(predators::Vector{String}, rng::AbstractRNG; nsweeps::Int=10)
    shuffled = copy(predators)
    n = length(shuffled)
    for _ in 1:(nsweeps * n)
        i, j = rand(rng, 1:n, 2)
        i == j && continue
        shuffled[i], shuffled[j] = shuffled[j], shuffled[i]
    end
    shuffled
end

function null_mean_absdiff(
    predators::Vector{String}, resources::Vector{String},
    trait_lookup::Dict{String,Float64}; nperm::Int=1000, nsweeps::Int=10, seed::Int=1
)
    rng = MersenneTwister(seed)
    values = Vector{Float64}(undef, nperm)
    for permutation in 1:nperm
        shuffled = edge_swap(predators, rng; nsweeps=nsweeps)
        values[permutation] = mean(
            abs(trait_lookup[shuffled[i]] - trait_lookup[resources[i]])
            for i in eachindex(resources)
        )
    end
    values
end

function bool_column(df::DataFrame, name::Symbol)
    col = df[!, name]
    if eltype(col) <: Union{Missing,Bool}
        return coalesce.(col, false)
    end
    lowercase.(string.(coalesce.(col, "false"))) .== "true"
end

master = CSV.read(MASTER_FILE, DataFrame; missingstring="", ntasks=1, pool=false)
coverage = isfile(COVERAGE_FILE) ?
    CSV.read(COVERAGE_FILE, DataFrame; missingstring="", ntasks=1, pool=false) :
    DataFrame()

:consumer_taxon_status in propertynames(master) || error(
    "Master dataset lacks consumer_taxon_status. Run ThermalAnalysis/BuildThermalDataset.R first."
)
all(lowercase.(String.(coalesce.(master.consumer_taxon_status, ""))) .== "animal") || error(
    "Master dataset contains a consumer that is not a verified animal."
)

summary_all = DataFrame(
    metric=String[], edges=Int[], predators=Int[], prey=Int[],
    obs_mean_absdiff=Float64[], null_mean=Float64[], null_sd=Float64[],
    z=Float64[], p_one_sided=Float64[], median_recorded_prey_coverage=Float64[]
)

for metric in TRAITS
    outdir = joinpath(OUTROOT, metric)
    isdir(outdir) || mkpath(outdir)

    eligible_col = Symbol("eligible_", metric)
    consumer_value_col = Symbol("consumer_trait_value_", metric)
    resource_value_col = Symbol("resource_trait_value_", metric)
    consumer_source_col = Symbol("consumer_selected_trait_source_", metric)
    resource_source_col = Symbol("resource_selected_trait_source_", metric)

    required = [eligible_col, consumer_value_col, resource_value_col,
                consumer_source_col, resource_source_col]
    missing_cols = setdiff(required, propertynames(master))
    isempty(missing_cols) || error("Master dataset lacks columns: $(join(string.(missing_cols), ", "))")

    keep = bool_column(master, eligible_col)
    selected = master[keep, :]
    edges = DataFrame(
        pred=String.(selected.consumer),
        prey=String.(selected.resource),
        interaction_sources=String.(coalesce.(selected.interaction_sources, "")),
        interaction_types=String.(coalesce.(selected.interaction_types, "")),
        study_citations=String.(coalesce.(selected.study_citations, "")),
        pred_trait=Float64.(selected[!, consumer_value_col]),
        prey_trait=Float64.(selected[!, resource_value_col]),
        trait_source_pred=String.(selected[!, consumer_source_col]),
        trait_source_prey=String.(selected[!, resource_source_col])
    )
    edges.absdiff = abs.(edges.pred_trait .- edges.prey_trait)
    sort!(edges, [:pred, :prey])
    CSV.write(joinpath(outdir, "edge_table.csv"), edges)

    predator_level = if isempty(edges)
        DataFrame(pred=String[], pred_trait=Float64[], mean_prey_trait=Float64[],
                  trait_source_pred=String[], nprey=Int[])
    else
        combine(groupby(edges, :pred),
            :pred_trait => first => :pred_trait,
            :prey_trait => mean => :mean_prey_trait,
            :trait_source_pred => first => :trait_source_pred,
            nrow => :nprey
        )
    end
    CSV.write(joinpath(outdir, "predator_level.csv"), predator_level)

    observed = isempty(edges) ? NaN : mean(edges.absdiff)
    null_values = Float64[]
    null_mean = NaN
    null_sd = NaN
    z = NaN
    p = NaN

    if nrow(edges) >= 2
        trait_lookup = Dict{String,Float64}()
        for row in eachrow(edges)
            trait_lookup[row.pred] = row.pred_trait
            trait_lookup[row.prey] = row.prey_trait
        end
        null_values = null_mean_absdiff(
            collect(edges.pred), collect(edges.prey), trait_lookup;
            nperm=N_PERM, nsweeps=N_SWEEPS, seed=SEED
        )
        null_mean = mean(null_values)
        null_sd = std(null_values)
        z = null_sd > 0 ? (observed - null_mean) / null_sd : NaN
        p = mean(null_values .<= observed)
    end
    CSV.write(joinpath(outdir, "null_absdiff_values.csv"), DataFrame(null=null_values))

    metric_coverage = isempty(coverage) ? DataFrame() : coverage[coverage.metric .== metric, :]
    if !isempty(metric_coverage)
        CSV.write(joinpath(outdir, "consumer_coverage.csv"), metric_coverage)
    end
    included_coverage = if isempty(metric_coverage) || !(:included_in_analysis in propertynames(metric_coverage))
        metric_coverage
    else
        metric_coverage[bool_column(metric_coverage, :included_in_analysis), :]
    end
    median_coverage = isempty(included_coverage) ? NaN :
        median(skipmissing(included_coverage.prey_trait_coverage))

    summary_row = (
        metric=metric, edges=nrow(edges), predators=length(unique(edges.pred)),
        prey=length(unique(edges.prey)), obs_mean_absdiff=observed,
        null_mean=null_mean, null_sd=null_sd, z=z, p_one_sided=p,
        median_recorded_prey_coverage=median_coverage
    )
    summary_df = DataFrame([summary_row])
    CSV.write(joinpath(outdir, "summary.csv"), summary_df)
    push!(summary_all, summary_row)
end

CSV.write(joinpath(OUTROOT, "summary_all_metrics.csv"), summary_all)
println("Saved thermal analyses to: $OUTROOT")
