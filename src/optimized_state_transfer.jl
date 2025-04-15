const MIN_TIME = 0
const MAX_TIME = 4π
const DEFAULT_TOL = 1e-5


"""
    OptimizedStateTransfer

A struct to store data on optimized state transfer between all qubit pairs in a network.

# Fields
- `adjacency_matrix::AbstractMatrix{<:Real}`: The network in adjacency_matrix format.
- `exhibits_pst::Bool`: Whether the quantum network exhibits perfect state transfer
between any pair of qubits.
- `qubit_pairs::Base.Generator`: A lazily evaluated generator of structs containing
data on state transfer between every pair of qubits in the network.
"""
struct OptimizedStateTransfer
    adjacency_matrix::AbstractMatrix{<:Real}
    exhibits_pst::Bool
    qubit_pairs::Base.Generator
end

function Base.show(io::IO, ost::OptimizedStateTransfer)
    buffer = IOBuffer()
    show(IOContext(buffer), "text/plain", ost.adjacency_matrix)
    s1 = String(take!(buffer)); s1 = s1[1:(findfirst(isequal(':'), s1) - 1)]
    s2 = ost.exhibits_pst
    s3 = "QubitPairTransfer Base.Generator"
    out = "OptimizedStateTransfer($s1, $s2, $s3)"
    print(io, out)
end


"""
    QubitPairTransfer

A struct to store data on optimized state transfer between two qubits in a quantum network.

# Fields
- `source::Int`: The index of the source qubit.
- `dest::Int`: The index of the target qubit.
- `is_pst::Bool`: Whether there is perfect state transfer from the source qubit to the
target qubit (i.e., whether `maximum_fidelity` is sufficiently close to `1`).
- `maximum_fidelity::Float64`: The (numerically approximated) maximum transfer fidelity from
the source qubit to the target qubit.
- `optimal_time::Float64`: The time at which the maximum transfer fidelity is achieved.
"""
struct QubitPairTransfer
    source::Int
    dest::Int
    is_pst::Bool
    maximum_fidelity::Float64
    optimal_time::Float64
end

function Base.show(io::IO, qpt::QubitPairTransfer)
    (s1, s2, s3) = (qpt.source, qpt.dest, qpt.is_pst)
    (s4, s5) = round.([qpt.maximum_fidelity, qpt.optimal_time], digits=5)
    out = "QubitPairTransfer($s1, $s2, $s3, $s4, $s5)"
    print(io, out)
end


# TODO: Finalize docstring
"""
    optimized_state_transfer(adj_mat::AbstractMatrix{<:Real}; kwargs...)
    optimized_state_transfer(graph::AbstractGraph{Int}; kwargs...)

Add later

# Arguments
Add later

# Keywords
See `qubit_pair_transfer` documentation for more information. (The set of valid keyword
arguments is exactly the same for both functions, as `qubit_pair_transfer` is a helper
function for `optimized_state_transfer`.)

# Returns
Add later

# Examples
Add later
```
"""
function optimized_state_transfer(adj_mat::AbstractMatrix{<:Real}; kwargs...)
    n = size(adj_mat, 1)
    uvs = [(u, v) for u in 1:(n - 1) for v in (u + 1):n]
    qubit_pairs = Base.Generator(
        uv -> qubit_pair_transfer(adj_mat, uv[1], uv[2]; kwargs...), uvs
    )
    exhibits_pst = any(qpt -> qpt.is_pst, qubit_pairs)
    return OptimizedStateTransfer(adj_mat, exhibits_pst, qubit_pairs)
end

@inline function optimized_state_transfer(graph::AbstractGraph{Int}; kwargs...)
    MatType = hasproperty(graph, :weights) ? Matrix : BitMatrix
    adj_mat = MatType(adjacency_matrix(graph))
    return optimized_state_transfer(adj_mat; kwargs...)
end


"""
    qubit_pair_transfer(
        adj_mat::AbstractMatrix{<:Real}, source::Int, dest::Int;
        min_time::Real=$MIN_TIME,
        max_time::Real=$MAX_TIME,
        tol::Real=$DEFAULT_TOL,
    )
    qubit_pair_transfer(graph::AbstractGraph{Int}, source::Int, dest::Int; kwargs...)

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
- `min_time`::Real=$MIN_TIME: The minimum time at which state transfer is considered.
- `max_time::Real=$MAX_TIME`: The maximum time at which state transfer is considered.
- `tol::Real=$DEFAULT_TOL`: The margin of error for determining whether perfect state transfer
occurs between the source and target.

# Returns
- `QubitPairTransfer`: A struct containing data on the optimized state transfer. Also
indicates whether PST was found between the two qubits.

# Examples
```jldoctest
julia> using Graphs: Graph

julia> C4_adj = BitMatrix([0 1 0 1; # Cycle graph on 4 vertices (adj. matrix format)
                           1 0 1 0;
                           0 1 0 1;
                           1 0 1 0]);

julia> C4_graph = Graph(C4_adj); # Cycle graph on 4 vertices (graph format)

julia> qubit_pair_transfer(C4_adj, 1, 2) # There is no PST from node 1 to node 2
QubitPairTransfer(1, 2, false, 0.25, 2.35619)

julia> qubit_pair_transfer(C4_graph, 1, 3) # There is PST from node 1 to node 3 over time π/2
QubitPairTransfer(1, 3, true, 1.0, 1.5708)
"""
function qubit_pair_transfer(
    adj_mat::AbstractMatrix{<:Real}, source::Int, dest::Int;
    min_time::Real=MIN_TIME,
    max_time::Real=MAX_TIME,
    tol::Real=DEFAULT_TOL,
)
    (adj_mat == adj_mat') || throw(DomainError(adj_mat, ADJ_MAT_ERR))
    
    n = size(adj_mat, 1)
    source_state = falses(n); source_state[source] = true
    dest_state = falses(n); dest_state[dest] = true
    fidelity(time::Float64) = abs2(source_state' * exp(im * time * adj_mat) * dest_state)
    
    lower = Float64(min_time)
    upper = Float64(max_time)
    lipschitz = 2 * norm(adj_mat, 2)
    
    res = maximize_shubert(fidelity, lower, upper, lipschitz; tol=tol)
    maximum_fidelity = res.maximum
    optimal_time = res.maximizer
    is_pst = maximum_fidelity > 1 - tol
    
    return QubitPairTransfer(source, dest, is_pst, maximum_fidelity, optimal_time)
end

@inline function qubit_pair_transfer(
    graph::AbstractGraph{Int}, source::Int, dest::Int; kwargs...
)
    MatType = hasproperty(graph, :weights) ? Matrix : BitMatrix
    adj_mat = MatType(adjacency_matrix(graph))
    return qubit_pair_transfer(adj_mat, source, dest; kwargs...)
end
