module ErrorMessages

export ADJ_MAT_ERR, INTERVAL_ERR, TOL_ERR, LIPSCHITZ_ERR

const ADJ_MAT_ERR = "Adjacency matrix must be symmetric"
const INTERVAL_ERR = "Optimization range must be non-empty"
const LIPSCHITZ_ERR = "Lipschitz constant must be positive"
const TOL_ERR = "Tolerance must be positive"

end
