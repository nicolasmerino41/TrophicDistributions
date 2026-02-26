module Sweep

using ..Parameters: CONNECTANCE_RANGE, CORR_RANGE, N_CONNECT, N_CORR, NREP, BASE_SEED, NX, NY, NCELLS, S, BASAL_FRAC, Emin_patch, USE_CONNECTIVITY_FILTER, OUTDIR
using ..Connectivity: make_workspaces
using ..Niches: regimes
using ..Simulation: simulate_one!
using Random
using Statistics

export sweep_all, ENVKINDS, NETFAMS, NETNAMES, regime_name

const ENVKINDS = [:random, :autocorr]
const NETFAMS  = [:random, :modular, :heavytail, :cascade]
const NETNAMES = Dict(
    :random=>"Random",
    :modular=>"Modular",
    :heavytail=>"Heavy-tail",
    :cascade=>"Cascade"
)

function regime_name(reg)
    return reg.name
end

function sweep_all()
    Cvals = collect(range(CONNECTANCE_RANGE[1], CONNECTANCE_RANGE[2], length=N_CONNECT))
    Rvals = collect(range(CORR_RANGE[1], CORR_RANGE[2], length=N_CORR))

    WSS = make_workspaces()

    metrics = [:dSrel, :mean_jaccard_mismatch, :frac_affected, :realized_overlap, :achieved_r, :Creal,
               :mismatch_q90, :mismatch_frac_gt]
    store = Dict{Tuple{Symbol,Symbol,Int,Symbol}, Matrix{Float64}}()

    for env in ENVKINDS, net in NETFAMS, (ri, reg) in enumerate(regimes), m in metrics
        store[(env, net, ri, m)] = fill(NaN, N_CORR, N_CONNECT)
    end

    println("Threads: ", Threads.nthreads())
    println("Grid: $(NX)x$(NY) cells=$(NCELLS), S=$(S), basal_frac=$(BASAL_FRAC), Emin=$(Emin_patch), connectivity_filter=$(USE_CONNECTIVITY_FILTER)")
    println("OUTDIR: ", OUTDIR)
    println("Sweep: env=$(length(ENVKINDS)) × net=$(length(NETFAMS)) × regimes=$(length(regimes)) × (C=$(N_CONNECT), r=$(N_CORR)) × reps=$(NREP)\n")

    total_cells = length(ENVKINDS)*length(NETFAMS)*length(regimes)*N_CONNECT*N_CORR
    cell_counter = 0

    for env in ENVKINDS
        for net in NETFAMS
            for (ri, reg) in enumerate(regimes)

                Threads.@threads for idx in 1:(N_CONNECT * N_CORR)
                    tid = Threads.threadid()
                    ws = WSS[tid]
                    rng = MersenneTwister(BASE_SEED + 100_000*tid + 17)

                    cind = (idx - 1) % N_CONNECT + 1
                    rind = (idx - 1) ÷ N_CONNECT + 1
                    C = Cvals[cind]
                    r = Rvals[rind]

                    dS = Float64[]
                    jm = Float64[]
                    fa = Float64[]
                    ov = Float64[]
                    rg = Float64[]
                    cr = Float64[]
                    q9 = Float64[]
                    fg = Float64[]

                    for rep in 1:NREP
                        seed = BASE_SEED +
                               10_000_000 * (findfirst(==(env), ENVKINDS)) +
                               1_000_000  * (findfirst(==(net), NETFAMS)) +
                               100_000    * ri +
                               10_000     * rind +
                               1_000      * cind +
                               rep
                        Random.seed!(rng, seed)
                        out = simulate_one!(rng, ws, env, net, reg, C, r)
                        push!(dS, out.dSrel)
                        push!(jm, out.mean_jaccard_mismatch)
                        push!(fa, out.frac_affected)
                        push!(ov, out.realized_overlap)
                        push!(rg, out.achieved_r)
                        push!(cr, out.Creal)
                        push!(q9, out.mismatch_q90)
                        push!(fg, out.mismatch_frac_gt)
                    end

                    store[(env, net, ri, :dSrel)][rind, cind] = mean(dS)
                    store[(env, net, ri, :mean_jaccard_mismatch)][rind, cind] = mean(jm)
                    store[(env, net, ri, :frac_affected)][rind, cind] = mean(fa)
                    store[(env, net, ri, :realized_overlap)][rind, cind] = mean(ov)
                    store[(env, net, ri, :achieved_r)][rind, cind] = mean(rg)
                    store[(env, net, ri, :Creal)][rind, cind] = mean(cr)
                    store[(env, net, ri, :mismatch_q90)][rind, cind] = mean(q9)
                    store[(env, net, ri, :mismatch_frac_gt)][rind, cind] = mean(fg)
                end

                cell_counter += N_CONNECT*N_CORR
                println("Done: env=$(env), net=$(net), regime=$(reg.name)  ($(cell_counter)/$(total_cells) cells)")
            end
        end
    end

    return store, Cvals, Rvals
end

end