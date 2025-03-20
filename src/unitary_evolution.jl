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

Simulates unitary evolution on a quantum network in the idealized absence of noise.

# Arguments
- `adj_mat::AbstractMatrix{<:Real}`: The quantum network in adjacency_matrix format.
- `graph::AbstractGraph{Int}`: The network in graph format. (Only passed when `adj_mat` is
not provided.)

# Keywords
- `time_steps::AbstractVector{<:Real}=1:(π / 200):5`: The times at which to track state
transfer between qubits.

# Returns
- `UnitaryEvolution`: A struct storing data on the unitary evolution of the quantum network.
Most notably, contains the amplitude and fidelity of state transfer between each pair of
qubits at every time step.


# Examples
```jldoctest
julia> using Graphs: cycle_graph

julia> c4_evolution = unitary_evolution(cycle_graph(4));

julia> c4_evolution.adjacency_matrix
4×4 BitMatrix:
 0  1  0  1
 1  0  1  0
 0  1  0  1
 1  0  1  0

julia> c4_evolution.time_steps
0.0:0.015707963267948967:4.995132319207771 # Equal to 0:(π / 200):5

julia> c4_evolution.transfer_amplitudes
319×4×4 Array{ComplexF64, 3}:
[:, :, 1] =
       1.0-0.0im  0.0-0.0im                 0.0-0.0im  0.0-0.0im
  0.999753-0.0im  0.0+0.0157054im   -0.00024672-0.0im  0.0+0.0157054im
  0.999013-0.0im  0.0+0.0313953im  -0.000986636-0.0im  0.0+0.0313953im
  0.997781-0.0im  0.0+0.0470542im   -0.00221902-0.0im  0.0+0.0470542im
  0.996057-0.0im  0.0+0.0626666im   -0.00394265-0.0im  0.0+0.0626666im
          ⋮
 0.0544967-0.0im  0.0-0.226995im      -0.945503-0.0im  0.0-0.226995im
 0.0618467-0.0im  0.0-0.240877im      -0.938153-0.0im  0.0-0.240877im
  0.069629-0.0im  0.0-0.254521im      -0.930371-0.0im  0.0-0.254521im
  0.077836-0.0im  0.0-0.267913im      -0.922164-0.0im  0.0-0.267913im

[:, :, 2] =
 0.0-0.0im              1.0-0.0im  0.0-0.0im                 0.0-0.0im
 0.0+0.0157054im   0.999753-0.0im  0.0+0.0157054im   -0.00024672-0.0im
 0.0+0.0313953im   0.999013-0.0im  0.0+0.0313953im  -0.000986636-0.0im
 0.0+0.0470542im   0.997781-0.0im  0.0+0.0470542im   -0.00221902-0.0im
 0.0+0.0626666im   0.996057-0.0im  0.0+0.0626666im   -0.00394265-0.0im
    ⋮
 0.0-0.226995im   0.0544967-0.0im  0.0-0.226995im      -0.945503-0.0im
 0.0-0.240877im   0.0618467-0.0im  0.0-0.240877im      -0.938153-0.0im
 0.0-0.254521im    0.069629-0.0im  0.0-0.254521im      -0.930371-0.0im
 0.0-0.267913im    0.077836-0.0im  0.0-0.267913im      -0.922164-0.0im

[:, :, 3] =
          0.0-0.0im  0.0-0.0im              1.0-0.0im  0.0-0.0im
  -0.00024672-0.0im  0.0+0.0157054im   0.999753-0.0im  0.0+0.0157054im
 -0.000986636-0.0im  0.0+0.0313953im   0.999013-0.0im  0.0+0.0313953im
  -0.00221902-0.0im  0.0+0.0470542im   0.997781-0.0im  0.0+0.0470542im
  -0.00394265-0.0im  0.0+0.0626666im   0.996057-0.0im  0.0+0.0626666im
             ⋮
    -0.945503-0.0im  0.0-0.226995im   0.0544967-0.0im  0.0-0.226995im
    -0.938153-0.0im  0.0-0.240877im   0.0618467-0.0im  0.0-0.240877im
    -0.930371-0.0im  0.0-0.254521im    0.069629-0.0im  0.0-0.254521im
    -0.922164-0.0im  0.0-0.267913im    0.077836-0.0im  0.0-0.267913im

[:, :, 4] =
 0.0-0.0im                 0.0-0.0im  0.0-0.0im              1.0-0.0im
 0.0+0.0157054im   -0.00024672-0.0im  0.0+0.0157054im   0.999753-0.0im
 0.0+0.0313953im  -0.000986636-0.0im  0.0+0.0313953im   0.999013-0.0im
 0.0+0.0470542im   -0.00221902-0.0im  0.0+0.0470542im   0.997781-0.0im
 0.0+0.0626666im   -0.00394265-0.0im  0.0+0.0626666im   0.996057-0.0im
    ⋮
 0.0-0.226995im      -0.945503-0.0im  0.0-0.226995im   0.0544967-0.0im
 0.0-0.240877im      -0.938153-0.0im  0.0-0.240877im   0.0618467-0.0im
 0.0-0.254521im      -0.930371-0.0im  0.0-0.254521im    0.069629-0.0im
 0.0-0.267913im      -0.922164-0.0im  0.0-0.267913im    0.077836-0.0im

julia> c4_evolution.transfer_fidelities
319×4×4 Array{Float64, 3}:
[:, :, 1] =
 1.0         0.0          0.0         0.0
 0.999507    0.000246659  6.08707e-8  0.000246659
 0.998028    0.000985662  9.7345e-7   0.000985662
 0.995567    0.00221409   4.92404e-6  0.00221409
 0.99213     0.0039271    1.55445e-5  0.0039271
 ⋮
 0.00296989  0.0515268    0.893976    0.0515268
 0.00382501  0.0580217    0.880132    0.0580217
 0.0048482   0.0647808    0.86559     0.0647808
 0.00605845  0.0717776    0.850386    0.0717776

[:, :, 2] =
 0.0          1.0         0.0          0.0
 0.000246659  0.999507    0.000246659  6.08707e-8
 0.000985662  0.998028    0.000985662  9.7345e-7
 0.00221409   0.995567    0.00221409   4.92404e-6
 0.0039271    0.99213     0.0039271    1.55445e-5
 ⋮
 0.0515268    0.00296989  0.0515268    0.893976
 0.0580217    0.00382501  0.0580217    0.880132
 0.0647808    0.0048482   0.0647808    0.86559
 0.0717776    0.00605845  0.0717776    0.850386

[:, :, 3] =
 0.0         0.0          1.0         0.0
 6.08707e-8  0.000246659  0.999507    0.000246659
 9.7345e-7   0.000985662  0.998028    0.000985662
 4.92404e-6  0.00221409   0.995567    0.00221409
 1.55445e-5  0.0039271    0.99213     0.0039271
 ⋮
 0.893976    0.0515268    0.00296989  0.0515268
 0.880132    0.0580217    0.00382501  0.0580217
 0.86559     0.0647808    0.0048482   0.0647808
 0.850386    0.0717776    0.00605845  0.0717776

[:, :, 4] =
 0.0          0.0         0.0          1.0
 0.000246659  6.08707e-8  0.000246659  0.999507
 0.000985662  9.7345e-7   0.000985662  0.998028
 0.00221409   4.92404e-6  0.00221409   0.995567
 0.0039271    1.55445e-5  0.0039271    0.99213
 ⋮
 0.0515268    0.893976    0.0515268    0.00296989
 0.0580217    0.880132    0.0580217    0.00382501
 0.0647808    0.86559     0.0647808    0.0048482
 0.0717776    0.850386    0.0717776    0.00605845


```
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

Tracks the wave function amplitude of a qubit in a quantum network over time.

# Arguments
- `adj_mat::AbstractMatrix{<:Real}`: The quantum network in adjacency_matrix format.
- `graph::AbstractGraph{Int}`: The network in graph format. (Only passed when `adj_mat` is
not provided.)
- `source::Int`: The qubit whose wave function amplitude to track.

# Keywords
- `dests::AbstractVector{Int}=1:size(adj_mat, 1)`: The qubits at which to track the wave
function amplitude of `source`.
- `time_steps::AbstractVector{<:Real}=0:(π / 200):5`: The times at which to track the wave
function amplitude of `source`.

# Returns
- `Matrix{ComplexF64}`: The transfer amplitudes of `source` to each qubit in `dests` at each
time step. The first dimension corresponds to the time step, and the second corresponds to
the target qubit.

# Examples
```jldoctest
julia> using Graphs: complete_multipartite_graph

julia> k122_node2_amps = track_qubit_amplitude(
    complete_multipartite_graph([1, 2, 2]), 2, time_steps=0:0.5:5
)
11×5 Matrix{ComplexF64}:
       0.0-0.0im              1.0-0.0im              0.0-0.0im               0.0-0.0im               0.0-0.0im
 -0.192803+0.352923im     0.68285-0.0697057im   -0.31715-0.0697057im  -0.0873013+0.35103im    -0.0873013+0.35103im
 -0.296067+0.190103im    0.238568-0.309653im   -0.761432-0.309653im   -0.0533586+0.144996im   -0.0533586+0.144996im
 0.0940871-0.00667217im  0.258737-0.280712im   -0.741263-0.280712im     0.253733-0.210152im     0.253733-0.210152im
   0.39497+0.180761im    0.460087+0.180299im   -0.539913+0.180299im     0.286909-0.198102im     0.286909-0.198102im
  0.170988+0.228892im    0.459578+0.412059im   -0.540422+0.412059im    -0.182253-0.0674034im   -0.182253-0.0674034im
 -0.026023-0.182558im    0.508058+0.0563556im  -0.491942+0.0563556im   -0.472027-0.0833521im   -0.472027-0.0833521im
  0.156815-0.418635im    0.721185-0.271338im   -0.278815-0.271338im    -0.155766+0.0571552im   -0.155766+0.0571552im
  0.156442-0.135118im    0.647642-0.113343im   -0.352358-0.113343im     0.220392+0.381336im     0.220392+0.381336im
 -0.260205+0.056111im    0.249513+0.107376im   -0.750487+0.107376im     0.205079+0.313436im     0.205079+0.313436im
 -0.421545-0.124698im    0.197874+0.0607882im  -0.802126+0.0607882im     0.11741-0.211222im      0.11741-0.211222im
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