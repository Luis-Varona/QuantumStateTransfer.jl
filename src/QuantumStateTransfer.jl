module QuantumStateTransfer

using Graphs: AbstractGraph, adjacency_matrix, nv
using LinearAlgebra: I, norm
using PrecompileTools: @setup_workload, @compile_workload

include("ErrorMessages.jl")
include("ShubertPiyavskii.jl")
using .ErrorMessages.QuantumStateTransfer
using .ShubertPiyavskii: maximize_shubert

include("QuantumStateTransfer/unitary_evolution.jl")
include("QuantumStateTransfer/optimized_state_transfer.jl")
include("QuantumStateTransfer/precompile_workload.jl")

export UnitaryEvolution, unitary_evolution, track_qubit_amplitude
export OptimizedStateTransfer, QubitPairTransfer,
    optimized_state_transfer, qubit_pair_transfer

end
