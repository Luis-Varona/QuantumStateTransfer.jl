module TestOptimizedStateTransfer

using Test
using QuantumStateTransfer: optimized_state_transfer, qubit_pair_transfer
using Graphs: Graph

C4_adj = BitMatrix([0 1 0 1; # Cycle graph on 4 vertices (adj. matrix format)
                    1 0 1 0;
                    0 1 0 1;
                    1 0 1 0]);
C4_graph = Graph(C4_adj); # Cycle graph on 4 vertices (graph format)

source = 1
dest1 = 2
dest2 = 3

maximum_fidelity1 = 1 / 4
is_pst2 = true
is_optimal_time1(time::Real) = time % 2π ≈ 3π / 4
is_optimal_time2(time::Real) = time % 2π ≈ π / 2

@time @testset "qubit_pair_transfer" begin
    qpt1 = qubit_pair_transfer(C4_adj, source, dest1)
    qpt2 = qubit_pair_transfer(C4_graph, source, dest2)
    
    @test qpt1.maximum_fidelity ≈ maximum_fidelity1
    @test qpt2.is_pst == is_pst2
    @test is_optimal_time1(qpt1.optimal_time)
    @test is_optimal_time2(qpt2.optimal_time)
end

# TODO: Make test for `optimized_state_transfer` function
@time @testset "optimized_state_transfer" begin
    @test 1 == 1
end

end