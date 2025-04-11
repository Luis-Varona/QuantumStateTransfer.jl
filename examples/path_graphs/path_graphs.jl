# %%
using QuantumStateTransfer: OptimizedStateTransfer, optimized_state_transfer
using Graphs: adjacency_matrix, edges, path_graph

function ost_to_string(ost::OptimizedStateTransfer)
    round5(x::Float64) = round(x, digits=5)
    return "$(ost.source) -> $(ost.dest) with fidelity " *
        "$(round5(ost.maximum_fidelity)) at time $(round5(ost.optimal_time))"
end

N = 30
dest = joinpath(@__DIR__, "path_graph_data.txt")

open(dest, "w+") do f
    write(f, "Path Graph Data\n\n")
end

# TODO: Refactor this so that it works with the new `OptimizedStateTransfer` API
for n in 2:25
    g = path_graph(n)
    g_adj = BitMatrix(adjacency_matrix(g))
    pairs = [(u, v) for u in 1:(n - 1) for v in (u + 1):n]
    ost_vec = map(
        pair -> optimized_state_transfer(g_adj, pair..., max_time=20, initial_time=2),
        pairs
    )
    best_fidelity = maximum(map(ost -> ost.maximum_fidelity, ost_vec))
    
    function criteria(ost::OptimizedStateTransfer)
        return isapprox(ost.maximum_fidelity, best_fidelity, atol=1e-2, rtol=1e-5)
    end
    
    best_pairs = filter(criteria, ost_vec)
    io = IOBuffer()
    write(io, "Highest-fidelity qubit pairs in P_$n:\n")
    
    for (i, ost) in enumerate(best_pairs)
        write(io, "($i) $(ost_to_string(ost))\n")
    end
    
    n < N && write(io, "\n")
    out = String(take!(io))
    print(out)
    
    open(dest, "a") do f
        write(f, out)
    end
end