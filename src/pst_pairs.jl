"""
    PSTPairs

A struct to store data on pairs of qubits in a network that exhibit perfect state transfer.

# Fields
- `adjacency_matrix::AbstractMatrix{<:Real}`: The network in adjacency_matrix format.
- `exhibits_pst::Bool`: Whether any pairs of qubits exhibit perfect state transfer.
- `pst_pairs::Base.Generator`: A generator of `PSTPair` structs representing (unordered)
pairs of qubits between which perfect state transfer occurs.
"""
struct PSTPairs
    adjacency_matrix::AbstractMatrix{<:Real}
    exhibits_pst::Bool
    pst_pairs::Base.Generator
end


"""
    PSTPair

A struct to store data on a pair of qubits that exhibit perfect state transfer.

# Fields
- `source::Int`: The index of the source qubit.
- `dest::Int`: The index of the target qubit.
- `maximum_fidelity::Float64`: The maximum transfer fidelity from the source to the target.
- `optimal_time::Float64`: The time at which the maximum transfer fidelity is achieved.
"""
struct PSTPair
    source::Int
    dest::Int
    maximum_fidelity::Float64
    optimal_time::Float64
end


"""
    pst_pairs(adj_mat::AbstractMatrix{<:Real}, tol::Real=1e-5)
    pst_pairs(graph::AbstractGraph{Int}, tol::Real=1e-5)

Generate all (unordered) pairs of qubits in a network that exhibit perfect state transfer.

# Arguments
- `adj_mat::AbstractMatrix{<:Real}`: The quantum network in adjacency_matrix format.
- `graph::AbstractGraph{Int}`: The network in graph format. (Only passed when `adj_mat` is
not provided.)

# Optional arguments
- `tol::Real=1e-5`: The margin of error for determining whether perfect state transfer occurs
between two qubits.

# Returns
- `PSTPairs`: A struct containing the adjacency matrix of the quantum network, a boolean
indicating whether any qubit pairs exhibit PST, and a generator of all PST pairs.
"""
function pst_pairs(adj_mat::AbstractMatrix{<:Real}, tol::Real=1e-5)
    function ost_result_to_pst_pair(result::OptimizedStateTransfer)
        return PSTPair(
            result.source, result.dest, result.maximum_fidelity, result.optimal_time
        )
    end
    
    n = size(adj_mat, 1)
    pairs = [(u, v) for u in 1:(n - 1) for v in (u + 1):n]
    optim_gen = (optimized_state_transfer(adj_mat, u, v, tol=tol) for (u, v) in pairs)
    pst_pairs = (ost_result_to_pst_pair(result) for result in optim_gen if result.is_pst)
    return PSTPairs(adj_mat, !isempty(pst_pairs), pst_pairs)
end

function pst_pairs(graph::AbstractGraph{Int}, tol::Real=1e-5)
    adj_mat = adjacency_matrix(graph)
    adj_mat = isa(graph, AbstractSimpleGraph) ? BitMatrix(adj_mat) : Matrix(adj_mat)
    return pst_pairs(adj_mat, tol)
end