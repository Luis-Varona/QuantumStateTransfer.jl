"""
    OptimizedStateTransfer

A struct to store data on optimized state transfer between two qubits in a quantum network.

# Fields
- `adjacency_matrix::AbstractMatrix{<:Real}`: The network in adjacency_matrix format.
- `source::Int`: The index of the source qubit.
- `dest::Int`: The index of the target qubit.
- `maximum_fidelity::Float64`: The (numerically approximated) maximum transfer fidelity from
the source qubit to the target qubit.
- `optimal_time::Float64`: The time at which the maximum transfer fidelity is achieved.
- `is_pst::Bool`: Whether there is perfect state transfer from the source qubit to the
target qubit (i.e., whether `maximum_fidelity` is sufficiently close to `1`).
"""
struct OptimizedStateTransfer
    adjacency_matrix::AbstractMatrix{<:Real}
    source::Int
    dest::Int
    maximum_fidelity::Float64
    optimal_time::Float64
    is_pst::Bool
end


"""
    optimized_state_transfer(
        adj_mat::AbstractMatrix{<:Real}, source::Int, dest::Int;
        min_time::Real=0, max_time::Real=4π, initial_time::Real=ℯ, tol::Real=1e-5,
    )
    optimized_state_transfer(
        graph::AbstractGraph{Int}, source::Int, dest::Int;
        min_time::Real=0, max_time::Real=4π, initial_time::Real=ℯ, tol::Real=1e-5,
    )

Find the maximized fidelity of state transfer between two qubits and the associated time.

# Arguments
- `adj_mat::AbstractMatrix{<:Real}`: The quantum network in adjacency_matrix format.
- `graph::AbstractGraph{Int}`: The network in graph format. (Only passed when `adj_mat` is
not provided.)
- `source::Int`: The index of the source qubit (assumed to have an initial state of the
`source`-th standard basis vector).
- `dest::Int`: The index of the target qubit (assumed to have an initial state of the
`dest`-th standard basis vector).

# Keywords
- `min_time::Real=0`: The minimum time at which state transfer is considered.
- `max_time::Real=4π`: The maximum time at which state transfer is considered.
- `initial_time::Real=ℯ`: The initial estimate passed to the optimization algorithm.
- `tol::Real=1e-5`: The margin of error for determining whether perfect state transfer
occurs between the source and target.

# Returns
- `OptimizedStateTransfer`: A struct containing data on the optimized state transfer.

# Examples
```jldoctest
julia> using Graphs

julia> C4_adj = BitMatrix([0 1 0 1; # Cycle graph on 4 vertices (adj. matrix format)
                           1 0 1 0;
                           0 1 0 1;
                           1 0 1 0]);

julia> C4_graph = Graph(C4_adj); # Cycle graph on 4 vertices (graph format)

julia> optimized_state_transfer(C4_adj, 1, 2) # There is no PST from node 1 to node 2
OptimizedStateTransfer(Bool[0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0], 1, 2, 0.25, 2.3561944921913542, false)

julia> optimized_state_transfer(C4_graph, 1, 3) # There is PST from node 1 to node 3 over time π/2
OptimizedStateTransfer(Bool[0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0], 1, 3, 1.0, 1.5707963268950111, true)


```
"""
function optimized_state_transfer(
    adj_mat::AbstractMatrix{<:Real}, source::Int, dest::Int;
    min_time::Real=0, max_time::Real=4π, initial_time::Real=ℯ, tol::Real=1e-5,
)
    if adj_mat != adj_mat'
        throw(ArgumentError("`adj_mat` must be symmetric"))
    end
    
    n = size(adj_mat, 1)
    source_state = falses(n); source_state[source] = true
    dest_state = falses(n); dest_state[dest] = true
    fidelity(time::Float64) = abs2(source_state' * exp(im * time * adj_mat) * dest_state)
    
    lower = [Float64(min_time)]
    upper = [Float64(max_time)]
    initial_x = [Float64(initial_time)]
    solver = IPNewton()
    
    res = optimize(x -> 1 - fidelity(x[1]), lower, upper, initial_x, solver)
    maximum_fidelity = 1 - res.minimum
    optimal_time = res.minimizer[1]
    is_pst = maximum_fidelity > 1 - tol
    
    return OptimizedStateTransfer(
        adj_mat, source, dest, maximum_fidelity, optimal_time, is_pst
    )
end

function optimized_state_transfer(
    graph::AbstractGraph{Int}, source::Int, dest::Int;
    min_time::Real=0, max_time::Real=4π, initial_time::Real=ℯ, tol::Real=1e-5,
)
    adj_mat = adjacency_matrix(graph)
    adj_mat = isa(graph, AbstractSimpleGraph) ? BitMatrix(adj_mat) : Matrix(adj_mat)
    return optimized_state_transfer(
        adj_mat, source, dest,
        min_time=min_time, max_time=max_time, initial_time=initial_time, tol=tol
    )
end