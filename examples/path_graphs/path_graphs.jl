# %%
using QuantumStateTransfer: QubitPairTransfer, optimized_state_transfer
using Graphs: adjacency_matrix, edges, path_graph

function qpt_to_string(qpt::QubitPairTransfer)
    round5(x::Float64) = round(x, digits=5)
    return "$(qpt.source) -> $(qpt.dest) with fidelity " *
        "$(round5(qpt.maximum_fidelity)) at time $(round5(qpt.optimal_time))"
end

N = 14
dest = joinpath(@__DIR__, "path_graph_data.txt")

open(dest, "w+") do f
    write(f, "Path Graph Data\n\n")
end

for n in 2:N
    g = path_graph(n)
    @time qpts = collect(optimized_state_transfer(g, max_time=1000, tol=0.05).qubit_pairs)
    best_fidelity = maximum(map(qpt -> qpt.maximum_fidelity, qpts))
    
    function criteria(qpt::QubitPairTransfer)
        return isapprox(qpt.maximum_fidelity, best_fidelity, atol=1e-2, rtol=1e-5)
    end
    
    best_pairs = filter(criteria, qpts)
    io = IOBuffer()
    write(io, "Highest-fidelity qubit pairs in P_$n:\n")
    
    for (i, qpt) in enumerate(best_pairs)
        write(io, "($i) $(qpt_to_string(qpt))\n")
    end
    
    n < N && write(io, "\n")
    out = String(take!(io))
    print(out)
    
    open(dest, "a") do f
        write(f, out)
    end
end
