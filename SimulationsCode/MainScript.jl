include(joinpath(@__DIR__, "PackageLoading.jl"))
include(joinpath(@__DIR__, "Parameters.jl"))
include(joinpath(@__DIR__, "Functions.jl"))

using .Functions.Sweep: sweep_all
using .Functions.IO: save_cache, load_cache, save_all_tsv
using .Parameters: OUTDIR

results = sweep_all()

cache_path = joinpath(OUTDIR, "cache.jls")
save_cache(cache_path, results)

# Load a previous cache 
# cache_path = joinpath(OUTDIR, "cache.jls")
# data = load_cache(cache_path)
# results = load_cache(cache_path)

paths = save_all_tsv(results)
println("Saved outputs:")
foreach(path -> println("  ", path), values(paths))
println("Done.")