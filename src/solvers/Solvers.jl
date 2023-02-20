module Solvers

include("acceleratedmethods/acceleratedmethods.jl")
include("sweepmethods/sweepmethods.jl")
export ScalingMethod, DecompositionMethod
export BoundaryPointsSM, BoundaryPointsDM
export evaluate_points, construct_matrices
export solve, solve_vect
export solve_wavenumber, solve_spectrum

end