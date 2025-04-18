@setup_workload begin
    using Graphs: cycle_graph
    g = cycle_graph(3)
    time_steps = 0:3
    min_time = 0.0
    max_time = 2.0
    tol = 0.1
    
    @compile_workload begin
        unitary_evolution(g, time_steps)
        ost = optimized_state_transfer(g; min_time=min_time, max_time=max_time, tol=tol)
        collect(ost.qubit_pairs) # Eager instantiation of `QubitPairTransfer` structs
    end
end
