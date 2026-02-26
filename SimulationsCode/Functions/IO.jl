module IO

using ..Parameters: OUTDIR
using ..Sweep: ENVKINDS, NETFAMS
using ..Niches: regimes
using Serialization
using Printf

export save_matrix_tsv, save_all_tsv, save_cache, load_cache

function save_cache(cache_path::String, data)
    serialize(cache_path, data)
end

function load_cache(cache_path::String)
    return deserialize(cache_path)
end

function save_matrix_tsv(path::String, M::Matrix{Float64})
    open(path, "w") do io
        nr, nc = size(M)
        for r in 1:nr
            println(io, join([@sprintf("%.6f", M[r,c]) for c in 1:nc], '\t'))
        end
    end
end

function save_all_tsv(store)
    mean_dir = joinpath(OUTDIR, "meanMetrics")
    tail_dir = joinpath(OUTDIR, "tailMetrics")
    isdir(mean_dir) || mkpath(mean_dir)
    isdir(tail_dir) || mkpath(tail_dir)

    for env in ENVKINDS, net in NETFAMS, (ri, reg) in enumerate(regimes)
        for metric in [:dSrel, :mean_jaccard_mismatch, :frac_affected, :realized_overlap, :achieved_r, :Creal]
            M = store[(env, net, ri, metric)]
            fname = "mat_$(env)_$(net)_reg$(ri)_$(metric).tsv"
            save_matrix_tsv(joinpath(mean_dir, fname), M)
        end
        for metric in [:mismatch_q90, :mismatch_frac_gt]
            M = store[(env, net, ri, metric)]
            fname = "mat_$(env)_$(net)_reg$(ri)_$(metric).tsv"
            save_matrix_tsv(joinpath(tail_dir, fname), M)
        end
    end
end

end