module ErrorMessages

module QuantumStateTransfer

export ADJ_MAT_ERR

const ADJ_MAT_ERR = "Adjacency matrix must be symmetric"
# TODO: I'm sure there will be more...

end # module QuantumStateTransfer


module ShubertPiyavskii

export LIPSCHITZ_ERR, OPTIM_RANGE_ERR, TOL_ERR

const LIPSCHITZ_ERR = "Lipschitz constant must be non-negative"
const OPTIM_RANGE_ERR = "Optimization range must be non-empty"
const TOL_ERR = "Tolerance must be positive"

end # module ShubertPiyavskii

end # module ErrorMessages
