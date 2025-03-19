"""
    TIME_STEPS::StepRangeLen{Float64}

The default times at which to track state amplitudes in a quantum network.
"""
const TIME_STEPS::StepRangeLen{Float64} = 0:(π / 200):5


"""
    UnitaryEvolution

A struct to store data on the unitary evolution of a quantum network.

# Fields
- `adjacency_matrix::AbstractMatrix{<:Real}`: The network in adjacency_matrix format.
- `time_steps::AbstractVector{<:Real}`: The times at which to track state amplitudes.
- `transfer_amplitudes::Array{ComplexF64, 3}`: The amplitudes of state transfer between
qubits. The first dimension corresponds to the time step, the second corresponds to the
target qubit, and the third corresponds to the source qubit.
- `transfer_fidelities::Array{Float64, 3}`: The magnitudes of the transfer amplitudes
squared, representing wave function probabilities rather than amplitudes.
"""
struct UnitaryEvolution
    adjacency_matrix::AbstractMatrix{<:Real}
    time_steps::AbstractVector{<:Real}
    transfer_amplitudes::Array{ComplexF64, 3}
    transfer_fidelities::Array{Float64, 3}
end


"""
    unitary_evolution(
        adj_mat::AbstractMatrix{<:Real}, time_steps::AbstractVector{<:Real}=1:(π / 200):5,
    )
    unitary_evolution(
        graph::AbstractGraph{Int}, time_steps::AbstractVector{<:Real}=1:(π / 200):5,
    )::UnitaryEvolution

[To do]

# Arguments
- `adj_mat::AbstractMatrix{<:Real}`: [To do]
- `graph::AbstractGraph{Int}`: [To do]

# Returns
- `UnitaryEvolution`: [To do]

# Examples
[To do]
"""
function unitary_evolution(
    adj_mat::AbstractMatrix{<:Real}, time_steps::AbstractVector{<:Real}=TIME_STEPS,
)
    amps_vec = map(u -> track_qubit_amplitude(adj_mat, u), 1:size(adj_mat, 1))
    transfer_amplitudes = stack(amps_vec)
    transfer_fidelities = abs2.(transfer_amplitudes)
    return UnitaryEvolution(adj_mat, time_steps, transfer_amplitudes, transfer_fidelities)
end

function unitary_evolution(
    graph::AbstractGraph{Int}, time_steps::AbstractVector{<:Real}=TIME_STEPS
)
    adj_mat = adjacency_matrix(graph)
    adj_mat = isa(graph, AbstractSimpleGraph) ? BitMatrix(adj_mat) : Matrix(adj_mat)
    return unitary_evolution(adj_mat, time_steps)
end


"""
    track_qubit_amplitude(
        adj_mat::AbstractMatrix{<:Real}, source::Int;
        dests::AbstractVector{Int}=1:size(adj_mat, 1),
        time_steps::AbstractVector{<:Real}=0:(π / 200):5,
    )
    track_qubit_amplitude(
        graph::AbstractGraph{Int}, source::Int;
        dests::AbstractVector{Int}=1:nv(graph),
        time_steps::AbstractVector{<:Real}=0:(π / 200):5,
    )

[To do]

# Arguments
- `adj_mat::AbstractMatrix{<:Real}`: [To do]
- `graph::AbstractGraph{Int}`: [To do]
- `source::Int`: [To do]

# Keywords
[To do]

# Returns
- `Matrix{ComplexF64}`: [To do]

# Examples
[To do]
"""
function track_qubit_amplitude(
    adj_mat::AbstractMatrix{<:Real}, source::Int;
    dests::AbstractVector{Int}=1:size(adj_mat, 1),
    time_steps::AbstractVector{<:Real}=TIME_STEPS,
)
    if adj_mat != adj_mat'
        throw(ArgumentError("`adj_mat` must be symmetric"))
    end
    
    identity_mat = I(size(adj_mat, 1))
    source_state = @view identity_mat[:, source]
    dest_states = @view identity_mat[:, dests]
    
    amps = map(time -> source_state' * exp(im * time * adj_mat) * dest_states, time_steps)
    return vcat(amps...)
end

function track_qubit_amplitude(
    graph::AbstractGraph{Int}, source::Int;
    dests::AbstractVector{Int}=1:nv(graph),
    time_steps::AbstractVector{<:Real}=TIME_STEPS,
)
    adj_mat = adjacency_matrix(graph)
    adj_mat = isa(graph, AbstractSimpleGraph) ? BitMatrix(adj_mat) : Matrix(adj_mat)
    return track_qubit_amplitude(adj_mat, source; dests=dests, time_steps=time_steps)
end