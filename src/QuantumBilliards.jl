module QuantumBilliards

using FFTW
using Bessels
using SpecialFunctions
using LinearAlgebra
using SparseArrays
using StaticArrays
using Arpack
using Random
using ForwardDiff
using QuadGK
using FastGaussQuadrature
using Optim
using ProgressMeter
using KrylovKit
using LinearMaps
using CoordinateTransformations
using Rotations
using CircularArrays
using Logging,TimerOutputs
using BilliardGeometry
import BilliardGeometry: boundary_coords

# abstract types
include("abstracttypes.jl")
export CoordinateSystem,AbsBasis,AbsSolver,AbsPoints,SweepSolver,AcceleratedSolver,AbsState,StationaryState,CFIE
include("utils/special_function_calls.jl")
include("utils/other.jl")

# utilities
include("utils/macros.jl")
export try_MKL!

include("utils/coordinatesystems.jl")
export CartesianCS,PolarCS,polar_to_cartesian,cartesian_to_polar

include("utils/kress_grading_single_corner.jl")
include("utils/kress_grading_multi_corner.jl")

# basis
include("basis/planewaves/realplanewaves.jl")
export RealPlaneWaves

include("basis/fourierbessel/corneradapted.jl")
include("basis/fourierbessel/multi_cafb.jl")
export CornerAdaptedFourierBessel,MultiCornerAdaptedFourierBessel
export resize_basis,basis_fun,dk_fun,gradient,basis_and_gradient

# boundary geometry
include("solvers/boundary_points.jl")
export BoundaryPoints,boundary_matrix_size,boundary_coords,boundary_s,component_offsets
export points_in_billiard
export kress_R!,kress_R_even!,kress_R_odd!
export BoundaryPanelArrays,BoundaryGeomCache,boundary_geom_cache
export component_lengths,component_normals,flatten_boundary_components,flatten_boundary_ds
export print_component_junctions

# symmetry
include("states/symmetry/symmetry.jl")
export DiagonalReflection,AntiDiagonalReflection
include("states/symmetry/symmetry_orbits.jl")
export SymmetryOrbitMap,symmetry_index_orbits,symmetry_orbit

# spectral geometry utilities
include("spectra/unfolding.jl")
export corner_correction,weyl_law,k_at_state,area

include("states/symmetry/reflections.jl")
export apply_symmetries_to_wavefunction,apply_symmetries_to_boundary_function,apply_symmetries_to_boundary_points

# matrix helpers
include("solvers/decompositions.jl")
export generalized_eigen,generalized_eigvals,generalized_eigen_all,adjust_scaling_and_samplers

include("solvers/matrixconstructors.jl")
export filter_matrix!,basis_matrix,gradient_matrices,basis_and_gradient_matrices,dk_matrix

# basis sweep methods
include("solvers/sweepmethods/basis_sweep/particular_solutions_method.jl")
export ParticularSolutionsMethod
export evaluate_points,construct_matrices,construct_matrices_benchmark
export solve_full,solve_with_rank_reduction,solve,solve_INFO,solve_vect

include("solvers/sweepmethods/basis_sweep/decomposition_method.jl")
export DecompositionMethodSolver

include("solvers/sweepmethods/basis_sweep/basis_sweep_methods.jl")
export solve_wavenumber,k_sweep

# direct boundary-integral method
include("solvers/sweepmethods/dlp/dlp.jl")
export BoundaryIntegralMethod,AbstractHankelBasis
export default_helmholtz_kernel_matrix,default_helmholtz_kernel_derivative_matrix,default_helmholtz_kernel_second_derivative_matrix
export compute_kernel_matrix,compute_kernel_matrix!,compute_kernel_matrix_with_derivatives!
export fredholm_matrix!,fredholm_matrix_with_derivatives!,fredholm_matrix,fredholm_matrix_with_derivatives
export adjoint_fredholm_matrix!,smallest_nullvec_krylov!,construct_matrices!

# DLP Kress
include("solvers/sweepmethods/dlp/dlp_kress.jl")
export DLP_kress,DLP_kress_global_corners
export DLPKressWorkspace,DLPKressReducedWorkspace
export build_Rmat_dlp_kress,build_dlp_kress_workspace_full,build_dlp_kress_reduced_workspace,build_dlp_kress_workspace
export construct_dlp_matrix!,construct_dlp_split!,construct_fredholm_matrix!
export construct_dlp_matrix_derivatives!,construct_fredholm_matrix_derivatives!

# CFIE Kress
include("solvers/sweepmethods/cfie/cfie_kress.jl")
export CFIE_kress,CFIE_kress_corners,CFIE_kress_global_corners,CFIE_kress_composite_solver
export CFIEKressWorkspace,build_cfie_kress_workspace,build_Rmat_kress,cfie_reduced_orbit_size
export construct_matrices_reduced!,construct_matrices_reduced_deriv!

include("solvers/sweepmethods/sweepmethods.jl")

# Chebyshev core
include("chebyshev/chebyshev_core.jl")
export _cheb_clenshaw,_chebfit!,_breaks_uniform

include("chebyshev/chebyshev_bessels.jl")
export ChebHankelTableH,ChebJTable,ChebHankelPlanH,ChebJPlan
export plan_h,plan_j,panel_indices,precompute_geom,panel_and_geom
export eval_h!,eval_j!,eval_h,eval_j,eval_h_multi_ks!,eval_j_multi_ks!
export h0_h1_j0_j1_multi_ks_at_r!,h0_h1_multi_ks_at_r!,h1_j1_multi_ks_at_r!
export h0_h1_h2_at_r,h0_h1_h2_multi_ks_at_r!,h1_multi_ks_at_r!,h1_at_r
export SLPWavefunctionChebPlan,CFIEWavefunctionChebPlan

# Chebyshev DLP
include("chebyshev/chebyshev_dlp.jl")
export DLPDerivChebWorkspace,DLPDerivativeChebyshevWorkspace
export compute_kernel_matrix_complex_k!,fredholm_matrix_complex_k!
export compute_kernel_matrices_DLP_chebyshev!,compute_kernel_matrices_DLP_chebyshev_derivatives!
export assemble_fredholm_matrices!,assemble_fredholm_matrices_with_derivatives!
export build_derivative_chebyshev_workspace
export construct_matrices_chebyshev!,construct_matrices_chebyshev_with_derivatives!
export construct_matrix_chebyshev_with_derivatives_at!
export adjoint_fredholm_matrix_from_bim_chebyshev!

# Chebyshev DLP Kress
include("chebyshev/chebyshev_dlp_kress.jl")
export DLPKressBlockCache,DLPKressSystemCache
export DLPKressH1J1BesselWorkspace,DLPKressH0H1J0J1BesselWorkspace
export DLPKressH1J1ChebWorkspace,DLPKressH0H1J0J1ChebWorkspace
export DLPKressReducedH1J1ChebWorkspace,DLPKressReducedH0H1J0J1ChebWorkspace
export DLPKressValueChebWorkspace,DLPKressDerivativeChebWorkspace
export build_dlp_kress_block_cache,build_dlp_kress_plans_h1_j1,build_dlp_kress_plans_h0_h1_j0_j1
export build_dlp_kress_h1_j1_cheb_workspace,build_dlp_kress_h0_h1_j0_j1_cheb_workspace
export construct_dlp_kress_matrices_chebyshev!,construct_dlp_kress_matrices_derivatives_chebyshev!
export adjoint_fredholm_matrix_from_dlp_chebyshev!

# Chebyshev CFIE Kress
include("chebyshev/chebyshev_cfie_kress.jl")
export CFIEKressBlockCache,CFIEKressSystemCache,CFIEKressReducedWorkspace
export CFIEKressH0H1J0J1BesselWorkspace,CFIEKressChebWorkspace
export build_cfie_kress_plans,build_cfie_kress_block_caches,build_cfie_kress_reduced_workspace,build_cfie_kress_cheb_workspace
export compute_kernel_matrices_CFIE_kress_chebyshev!
export construct_matrix_chebyshev_at!,construct_matrix_chebyshev_with_derivatives_at!

# optimal Chebyshev panelization
include("chebyshev/chebyshev_optimal_panelization.jl")
export chebyshev_params,construct_boundary_matrices!

# Vergini-Saraceno
include("solvers/acceleratedmethods/vergini_saraceno.jl")
export AbsScalingMethod,VerginiSaracenoSolver
export sm_results,solve_vectors,match_wavenumbers,match_wavenumbers_with_X
export overlap_and_merge!,overlap_and_merge_state!
export SpectralData,StateData,solve_state_data_bundle,compute_spectrum

# EBIM
include("solvers/acceleratedmethods/ebim.jl")
export EBIMSolver,EBIMChebBatchCache,build_ebim_cheb_cache
export construct_ebim_cheb_matrices!,construct_ebim_cheb_matrices,construct_ebim_cheb_matrix_at!
export solve_full!,solve_krylov!,solve!
export overlap_and_merge_ebim!,solve_spectrum_ebim,ebim_inv_diff

# Beyn
include("solvers/acceleratedmethods/beyn.jl")
export BeynSolver
export weyl_window_width,plan_weyl_windows,beyn_disks_from_windows,beyn_buffer_matrices
export construct_B_matrix,residual_and_norm_select,imag_k_check,solve_spectrum_beyn

# common accelerated-method interface
include("solvers/acceleratedmethods/accelerated_methods.jl")
export solve_wavenumber_beyn,solve_wavenumber_ebim,solve_spectrum

# states
include("states/eigenstates.jl")
export Eigenstate,BasisState,compute_eigenstate

include("states/boundary_and_layer_density_functions.jl")
export regularize!,boundary_function,momentum_function,symmetrize_layer_density

include("states/wavefunctions.jl")
export boundary_limits,wavefunction,wavefunctions,compute_psi
export CFIEWavefunctionCache

include("states/husimifunctions.jl")
export husimi_function

end