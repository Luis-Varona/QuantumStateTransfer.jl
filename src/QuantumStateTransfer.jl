module QuantumStateTransfer

using Graphs: AbstractGraph, adjacency_matrix, nv
using LinearAlgebra: I, norm

include("ErrorMessages.jl")
include("ShubertPiyavskii.jl")
using .ErrorMessages
using .ShubertPiyavskii: maximize_shubert

include("unitary_evolution.jl")
include("optimized_state_transfer.jl")

export UnitaryEvolution, unitary_evolution, track_qubit_amplitude
export OptimizedStateTransfer, QubitPairTransfer,
    optimized_state_transfer, qubit_pair_transfer

end
