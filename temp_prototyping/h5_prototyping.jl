using QuantumStateTransfer
using Graphs: Graph

function read_matrix(source::String, m::Int, n::Int)
    A = BitMatrix(undef, m, n)
    
    open(source, "r") do file
        for i in 1:m
            A[i, :] = parse.(Bool, split(readline(file)))
        end
    end
    
    return A
end

hypercube_5_source = "temp_prototyping/hypercube_5.txt"; n = 32
h5_adj_mat = read_matrix(hypercube_5_source, n, n)
h5_graph = Graph(h5_adj_mat)
h5_pst_pairs = pst_pairs(h5_graph)
h5_pair_gen = h5_pst_pairs.pst_pairs;
h5_pair_vec = collect(h5_pair_gen)

for pst_pair in h5_pair_vec
    println(pst_pair)
end

h5_unitary_evol = unitary_evolution(h5_graph)
println(h5_unitary_evol)