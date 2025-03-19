module TestOptimizedStateTransfer

using Test
using QuantumStateTransfer: optimized_state_transfer
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

@time @testset "optimized_state_transfer" begin
    results1 = optimized_state_transfer(C4_adj, source, dest1)
    results2 = optimized_state_transfer(C4_graph, source, dest2)
    
    @test results1.maximum_fidelity ≈ maximum_fidelity1
    @test results2.is_pst == is_pst2
    @test is_optimal_time1(results1.optimal_time)
    @test is_optimal_time2(results2.optimal_time)
end

end