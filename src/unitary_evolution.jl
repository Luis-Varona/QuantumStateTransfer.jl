const TIME_STEPS::StepRangeLen{Float64} = 0:(π / 200):5
const TIME_STEPS_STR::String = "0:(π / 200):5"


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
        adj_mat::AbstractMatrix{<:Real}, time_steps::AbstractVector{<:Real}=$TIME_STEPS_STR,
    )
    unitary_evolution(
        graph::AbstractGraph{Int}, time_steps::AbstractVector{<:Real}=$TIME_STEPS_STR,
    )

Simulates unitary evolution on a quantum network in the idealized absence of noise.

# Arguments
- `adj_mat::AbstractMatrix{<:Real}`: The quantum network in adjacency_matrix format.
- `graph::AbstractGraph{Int}`: The network in graph format. (Only passed when `adj_mat` is
not provided.)

# Keywords
- `time_steps::AbstractVector{<:Real}=$TIME_STEPS_STR`: The times at which to track state
transfer between qubits.

# Returns
- `UnitaryEvolution`: A struct storing data on the unitary evolution of the quantum network.
Most notably, contains the amplitude and fidelity of state transfer between each pair of
qubits at every time step.


# Examples
```jldoctest
julia> using Graphs: complete_bipartite_graph

julia> g = complete_bipartite_graph(1, 2); time_steps = LinRange(0, π, 7);

julia> evolution = unitary_evolution(g, time_steps); # Simulate state transfer on K₁₂

julia> evolution.adjacency_matrix # Adjacency matrix of the complete bipartite graph K₁₂
3×3 BitMatrix:
 0  1  1
 1  0  0
 1  0  0

julia> evolution.time_steps # Time steps at which to measure state transfer
7-element LinRange{Float64, Int64}:
 0.0, 0.523599, 1.0472, 1.5708, 2.0944, 2.61799, 3.14159

julia> evolution.transfer_amplitudes # State transfer amplitudes between qubits
7×3×3 Array{ComplexF64, 3}:
[:, :, 1] =
       1.0-0.0im  0.0-0.0im       0.0-0.0im
  0.738144-0.0im  0.0+0.477044im  0.0+0.477044im
 0.0897146-0.0im  0.0+0.704255im  0.0+0.704255im
   -0.6057-0.0im  0.0+0.56264im   0.0+0.56264im
 -0.983903-0.0im  0.0+0.126364im  0.0+0.126364im
 -0.846825-0.0im  0.0-0.37609im   0.0-0.37609im
 -0.266255-0.0im  0.0-0.681582im  0.0-0.681582im

[:, :, 2] =
 0.0-0.0im             1.0-0.0im        0.0-0.0im
 0.0+0.477044im   0.869072-0.0im  -0.130928-0.0im
 0.0+0.704255im   0.544857-0.0im  -0.455143-0.0im
 0.0+0.56264im     0.19715-0.0im   -0.80285-0.0im
 0.0+0.126364im  0.0080487-0.0im  -0.991951-0.0im
 0.0-0.37609im   0.0765877-0.0im  -0.923412-0.0im
 0.0-0.681582im   0.366872-0.0im  -0.633128-0.0im

[:, :, 3] =
 0.0-0.0im             0.0-0.0im        1.0-0.0im
 0.0+0.477044im  -0.130928-0.0im   0.869072-0.0im
 0.0+0.704255im  -0.455143-0.0im   0.544857-0.0im
 0.0+0.56264im    -0.80285-0.0im    0.19715-0.0im
 0.0+0.126364im  -0.991951-0.0im  0.0080487-0.0im
 0.0-0.37609im   -0.923412-0.0im  0.0765877-0.0im
 0.0-0.681582im  -0.633128-0.0im   0.366872-0.0im

julia> evolution.transfer_fidelities # State transfer fidelities between qubits
7×3×3 Array{Float64, 3}:
[:, :, 1] =
 1.0        0.0        0.0
 0.544857   0.227571   0.227571
 0.0080487  0.495976   0.495976
 0.366872   0.316564   0.316564
 0.968064   0.0159678  0.0159678
 0.717112   0.141444   0.141444
 0.0708919  0.464554   0.464554

[:, :, 2] =
 0.0        1.0         0.0
 0.227571   0.755287    0.0171421
 0.495976   0.296869    0.207155
 0.316564   0.0388681   0.644568
 0.0159678  6.47816e-5  0.983967
 0.141444   0.00586567  0.85269
 0.464554   0.134595    0.400851

[:, :, 3] =
 0.0        0.0        1.0
 0.227571   0.0171421  0.755287
 0.495976   0.207155   0.296869
 0.316564   0.644568   0.0388681
 0.0159678  0.983967   6.47816e-5
 0.141444   0.85269    0.00586567
 0.464554   0.400851   0.134595
```
"""
function unitary_evolution(
    adj_mat::AbstractMatrix{<:Real}, time_steps::AbstractVector{<:Real}=TIME_STEPS,
)
    amps_vec = map(
        u -> track_qubit_amplitude(adj_mat, u, time_steps=time_steps),
        1:size(adj_mat, 1)
    )
    transfer_amplitudes = stack(amps_vec)
    transfer_fidelities = abs2.(transfer_amplitudes)
    return UnitaryEvolution(adj_mat, time_steps, transfer_amplitudes, transfer_fidelities)
end

@inline function unitary_evolution(
    graph::AbstractGraph{Int}, time_steps::AbstractVector{<:Real}=TIME_STEPS,
)
    MatType = hasproperty(graph, :weights) ? Matrix : BitMatrix
    adj_mat = MatType(adjacency_matrix(graph))
    return unitary_evolution(adj_mat, time_steps)
end


"""
    track_qubit_amplitude(
        adj_mat::AbstractMatrix{<:Real}, source::Int;
        dests::AbstractVector{Int}=1:size(adj_mat, 1),
        time_steps::AbstractVector{<:Real}=$TIME_STEPS_STR,
    )
    track_qubit_amplitude(graph::AbstractGraph{Int}, source::Int; kwargs...)

Tracks the wave function amplitude of a qubit in a quantum network over time.

# Arguments
- `adj_mat::AbstractMatrix{<:Real}`: The quantum network in adjacency_matrix format.
- `graph::AbstractGraph{Int}`: The network in graph format. (Only passed when `adj_mat` is
not provided.)
- `source::Int`: The qubit whose wave function amplitude to track.

# Keywords
- `dests::AbstractVector{Int}=1:size(adj_mat, 1)`: The qubits at which to track the wave
function amplitude of `source`.
- `time_steps::AbstractVector{<:Real}=$TIME_STEPS_STR`: The times at which to track the wave
function amplitude of `source`.

# Returns
- `Matrix{ComplexF64}`: The transfer amplitudes of `source` to each qubit in `dests` at each
time step. The first dimension corresponds to the time step, and the second corresponds to
the target qubit.

# Examples
```jldoctest
julia> using Graphs: wheel_graph

julia> g = wheel_graph(8); source = 3; dests = [2, 4, 7]; time_steps = 0:1:5;

julia> amps_3_to_247 = track_qubit_amplitude(g, source, dests=dests, time_steps=time_steps)
6×3 Matrix{ComplexF64}:
        0.0-0.0im               0.0-0.0im               0.0-0.0im
 -0.0282559+0.340869im   -0.0282559+0.340869im   0.00691979-0.364799im
  0.0270747+0.159635im    0.0270747+0.159635im     0.353068-0.20453im
  -0.248756-0.288904im    -0.248756-0.288904im     0.291226-0.129163im
  -0.158849+0.0938213im   -0.158849+0.0938213im   -0.210858+0.121758im
   0.500422+0.160662im     0.500422+0.160662im    -0.259143-0.0882654im
```
"""
function track_qubit_amplitude(
    adj_mat::AbstractMatrix{<:Real}, source::Int;
    dests::AbstractVector{Int}=1:size(adj_mat, 1),
    time_steps::AbstractVector{<:Real}=TIME_STEPS,
)
    (adj_mat == adj_mat') || throw(ArgumentError("`adj_mat` must be symmetric"))
    
    identity_mat = I(size(adj_mat, 1))
    source_state = @view identity_mat[:, source]
    dest_states = @view identity_mat[:, dests]
    
    amps = map(time -> source_state' * exp(im * time * adj_mat) * dest_states, time_steps)
    return vcat(amps...)
end

@inline function track_qubit_amplitude(graph::AbstractGraph{Int}, source::Int; kwargs...)
    MatType = hasproperty(graph, :weights) ? Matrix : BitMatrix
    adj_mat = MatType(adjacency_matrix(graph))
    return track_qubit_amplitude(adj_mat, source; kwargs...)
end