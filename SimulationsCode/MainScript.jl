include(joinpath(@__DIR__, "Functions.jl"))

using .Functions.Sweep: sweep_all
using .Functions.IO: save_all_tsv

results = sweep_all()

paths = save_all_tsv(results)
println("Saved outputs:")
foreach(path -> println("  ", path), values(paths))
println("Done.")