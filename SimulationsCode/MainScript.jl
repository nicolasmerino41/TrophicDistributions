include(joinpath(@__DIR__, "PackageLoading.jl"))
include(joinpath(@__DIR__, "Parameters.jl"))
include(joinpath(@__DIR__, "Functions.jl"))

using .Functions.Sweep: sweep_all
using .Functions.IO: save_cache, load_cache, save_all_tsv
using .Parameters: OUTDIR

store, Cvals, Rvals = sweep_all()

cache_path = joinpath(OUTDIR, "cache.jls")
save_cache(cache_path, (store=store, Cvals=Cvals, Rvals=Rvals))

# Load a previous cache 
# cache_path = joinpath(OUTDIR, "cache.jls")
# data = load_cache(cache_path)
# store = data.store
# Cvals = data.Cvals
# Rvals = data.Rvals

save_all_tsv(store)

println("Done.")