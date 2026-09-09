module QuantumBilliards
using Bessels
using SpecialFunctions
using QuadGK
using CoordinateTransformations, Rotations
using LinearAlgebra, StaticArrays, CircularArrays
using Optim
using FFTW
using Logging, TimerOutputs
using Random, Distributions
using BilliardGeometry
using KrylovKit

#abstract types
include("abstracttypes.jl")
export AbsBasis, AbsSolver, AbsBasisSolver, AbsBIMSolver

include("utils/coordinatesystems.jl")
include("utils/geometryutils.jl")
include("utils/typeutils.jl")
include("utils/macros.jl")
include("utils/billiardutils.jl")
export make_triangle_and_basis, adapt_basis

include("basis/planewaves/realplanewaves.jl")
export RealPlaneWaves
include("basis/fourierbessel/corneradapted.jl")
export CornerAdaptedFourierBessel
export resize_basis, basis_fun, dk_fun, gradient, basis_and_gradient 

include("solvers/boundarypoints.jl")
export BoundaryPoints
export boundary_s, component_offsets, points_in_billiard
export estimate_rmin_rmax
include("solvers/boundarygeomcache.jl")
export BoundaryPanelArrays, BoundaryGeomCache, boundary_geom_cache, component_normals
export flatten_boundary_components, flatten_boundary_ds
include("solvers/decompositions.jl")
include("solvers/matrixconstructors.jl")
export basis_matrix, basis_and_gradient_matrices, dk_matrix

include("solvers/acceleratedmethods/acceleratedmethods.jl")
include("solvers/sweepmethods/sweepmethods.jl")
export SweepBasisSolver, AcceleratedBasisSolver
export SweepBIMSolver, AcceleratedBIMSolver
export VerginiSaracenoSolver, print_benchmark_info
export DecompositionMethodSolver
export ParticularSolutionsMethod
export BoundaryGrading, SmoothPeriodicGrading, CornerGrading, GlobalCornerGrading
export DoubleLayerPotentialSolver, CombinedFieldIntegralEquationSolver, CompositeBIMSolver
export ExpandedBIMSolver, BeynSolver
export BoundaryPointsSM, BoundaryPointsDM
export evaluate_points, construct_matrices, construct_matrices_benchmark
export solve, solve_vect, solve_vectors, solve_state
export solve_wavenumber, solve_spectrum
export k_sweep
export boundary_matrix_size
export weyl_window_width, plan_weyl_windows, beyn_disks_from_windows, beyn_buffer_matrices

include("spectra/spectralutils.jl")
export SpectralData, compute_spectrum, merge_spectra, overlap_and_merge!
include("spectra/unfolding.jl")
export weyl_law, area, fundamental_area

include("states/eigenstates.jl")
include("states/basisstates.jl")

export BasisEigenstate, BasisState
export compute_eigenstate
export BIMEigenstate
include("states/symmetry/reflections.jl")
include("states/wavefunctions.jl")
include("states/boundaryfunctions.jl")
include("states/husimifunctions.jl")

export wavefunction, compute_psi, boundary_limits #wavefunction_norm 
export get_boundary_curves_with_ignored, boundary_function, momentum_function, husimi_function


end