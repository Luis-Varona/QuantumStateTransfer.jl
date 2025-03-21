using QuantumStateTransfer: PSTPair, pst_pairs
using Graphs: complete_graph, cartesian_product


# %%
function hypercube_graph(n::Int)
    k2 = complete_graph(2)
    hc = k2
    
    for _ in 1:n - 1
        hc = cartesian_product(hc, k2)
    end
    
    return hc
end

for n in 1:5
    hc = hypercube_graph(n)
    n == 1 || println()
    @info "Finding all PST pairs for the $n-hypercube:"
    
    @time begin
        for (i, pair) in enumerate(collect(pst_pairs(hc).pst_pairs))
            println("($i) $pair")
        end
    end
end


# %%
using QuantumStateTransfer: optimized_state_transfer
using Graphs: adjacency_matrix, complete_multipartite_graph

G = complete_multipartite_graph([1, 4, 9])
G_adj = BitMatrix(adjacency_matrix(G))
n = size(G_adj, 1)
pairs = [(u, v) for u in 1:(n - 1) for v in (u + 1):n]
optim_vec = [optimized_state_transfer(G_adj, u, v, max_time=50) for (u, v) in pairs]

for (i, result) in enumerate(optim_vec)
    info = (
        result.source,
        result.dest,
        round(result.maximum_fidelity, digits=4),
        round(result.optimal_time, digits=4)
    )
    println("($i) $info")
end


# %%
G = hypercube_graph(5)
G_adj = BitMatrix(adjacency_matrix(G))
n = size(G_adj, 1)
pairs = [(u, v) for u in 1:(n - 1) for v in (u + 1):n]
optim_vec = [optimized_state_transfer(G_adj, u, v, max_time=50) for (u, v) in pairs]

for (i, result) in enumerate(optim_vec)
    info = (
        result.source,
        result.dest,
        round(result.maximum_fidelity, digits=4),
        round(result.optimal_time, digits=4)
    )
    println("($i) $info")
end