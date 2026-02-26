module Connectivity

using ..Parameters: NCELLS, USE_CONNECTIVITY_FILTER, Emin_patch
using ..Grid: NEIGH_4

export CCWorkspace, make_workspaces,
       lcc_size, largest_component_mask,
       keep_components_ge_Emin, apply_connectivity_filter

mutable struct CCWorkspace
    seen::Vector{Int32}
    stamp::Int32
    queue::Vector{Int}
end

function make_workspaces()
    nt = Threads.maxthreadid()
    wss = Vector{CCWorkspace}(undef, nt)
    for t in 1:nt
        wss[t] = CCWorkspace(fill(Int32(0), NCELLS), Int32(0), Int[])
    end
    return wss
end

function lcc_size(ws::CCWorkspace, mask::BitVector)
    count(mask) == 0 && return 0
    ws.stamp += 1
    stamp = ws.stamp
    seen = ws.seen
    empty!(ws.queue)

    best = 0
    @inbounds for i in 1:NCELLS
        if mask[i] && seen[i] != stamp
            seen[i] = stamp
            push!(ws.queue, i)
            qpos = 1
            compsize = 0
            while qpos <= length(ws.queue)
                v = ws.queue[qpos]; qpos += 1
                compsize += 1
                for nb in NEIGH_4[v]
                    if mask[nb] && seen[nb] != stamp
                        seen[nb] = stamp
                        push!(ws.queue, nb)
                    end
                end
            end
            best = max(best, compsize)
            empty!(ws.queue)
        end
    end
    return best
end

function largest_component_mask(ws::CCWorkspace, mask::BitVector)
    count(mask) == 0 && return BitVector(falses(NCELLS))
    ws.stamp += 1
    stamp = ws.stamp
    seen = ws.seen
    empty!(ws.queue)

    best = 0
    best_nodes = Int[]
    tmp_nodes = Int[]

    @inbounds for i in 1:NCELLS
        if mask[i] && seen[i] != stamp
            seen[i] = stamp
            push!(ws.queue, i)
            qpos = 1
            empty!(tmp_nodes)
            while qpos <= length(ws.queue)
                v = ws.queue[qpos]; qpos += 1
                push!(tmp_nodes, v)
                for nb in NEIGH_4[v]
                    if mask[nb] && seen[nb] != stamp
                        seen[nb] = stamp
                        push!(ws.queue, nb)
                    end
                end
            end
            if length(tmp_nodes) > best
                best = length(tmp_nodes)
                best_nodes = copy(tmp_nodes)
            end
            empty!(ws.queue)
        end
    end

    out = BitVector(falses(NCELLS))
    @inbounds for v in best_nodes
        out[v] = true
    end
    return out
end

function keep_components_ge_Emin(ws::CCWorkspace, mask::BitVector, Emin::Int)
    count(mask) == 0 && return BitVector(falses(NCELLS))

    ws.stamp += 1
    stamp = ws.stamp
    seen = ws.seen
    empty!(ws.queue)

    out = BitVector(falses(NCELLS))
    comp_nodes = Int[]

    @inbounds for i in 1:NCELLS
        if mask[i] && seen[i] != stamp
            empty!(comp_nodes)
            seen[i] = stamp
            push!(ws.queue, i)
            qpos = 1

            while qpos <= length(ws.queue)
                v = ws.queue[qpos]; qpos += 1
                push!(comp_nodes, v)
                for nb in NEIGH_4[v]
                    if mask[nb] && seen[nb] != stamp
                        seen[nb] = stamp
                        push!(ws.queue, nb)
                    end
                end
            end

            if length(comp_nodes) >= Emin
                for v in comp_nodes
                    out[v] = true
                end
            end

            empty!(ws.queue)
        end
    end

    return out
end

function apply_connectivity_filter(ws::CCWorkspace, mask::BitVector, Emin::Int)
    if !USE_CONNECTIVITY_FILTER
        return mask
    end

    if count(mask) < Emin
        return BitVector(falses(NCELLS))
    end

    kept = keep_components_ge_Emin(ws, mask, Emin)
    return (count(kept) == 0) ? BitVector(falses(NCELLS)) : kept
end

end