using QuantumStateTransfer

@info "Testing optimized_state_transfer"
include("optimized_state_transfer.jl")

println()

@info "Testing for PST on hypercube graphs"
include("pst_pairs.jl")