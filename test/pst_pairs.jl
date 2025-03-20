module TestPSTPairs

using Test
using QuantumStateTransfer: pst_pairs
using Graphs

function hypercube_graph(n::Int)
    k2 = complete_graph(2)
    hc = k2
    
    for _ in 1:n - 1
        hc = cartesian_product(hc, k2)
    end
    
    return hc
end

function hypercube_pst_pairs(n::Int)
    hc = hypercube_graph(n)
    hc_adj = adjacency_matrix(hc)
    first_rows = map(k -> (hc_adj^k)[1, :], 1:n)
    nodes_exclude = vcat(map(v -> getproperty(v, :nzind), first_rows[1:(n - 1)])...)
    nodes = filter(u -> u != 1, setdiff(first_rows[n].nzind, nodes_exclude))
    return map(i -> (1, i), nodes)
end

for n in 1:5
    hc = hypercube_graph(n)
    theoretical = hypercube_pst_pairs(n)
    
    @time @testset "pst_hypercube_$n" begin
        output = map(pair -> (pair.source, pair.dest), collect(pst_pairs(hc).pst_pairs))
        @test isempty(setdiff(theoretical, output))
    end
end

end