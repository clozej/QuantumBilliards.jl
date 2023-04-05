module QuantumBilliards
#using Reexport

#abstract types
include("abstracttypes.jl")
#utils must be included here so modules work
#export AbsBasis
include("utils/coordinatesystems.jl")
include("utils/geometryutils.jl")
include("utils/billiardutils.jl")
include("utils/typeutils.jl")
include("utils/gridutils.jl")
include("utils/symmetry.jl")
export Reflection, XReflection, YReflection, XYReflection
export real_length, is_inside

include("billiards/boundarypoints.jl")
export BoundaryPoints
export boundary_coords, dilated_boundary_points
#basis
#include("basis/Basis.jl")
#@reexport using .Basis
include("basis/planewaves/realplanewaves.jl")
export RealPlaneWaves
include("basis/fourierbessel/corneradapted.jl")
export CornerAdaptedFourierBessel
include("basis/fundamental/fundamentalbessels.jl")
export FundamentalBessel
export resize_basis, basis_fun, dk_fun, gradient, basis_and_gradient 

#billiards
#include("billiards/Billiards.jl")
#@reexport using .Billiards

include("billiards/geometry/geometry.jl")
export LineSegment, VirtualLineSegment
export CircleSegment, VirtualCircleSegment
export DispersingCircleSegment, VirtualDispersingCircleSegment
export curve, arc_length, tangent, tangent_vec, normal_vec
include("billiards/stadium.jl")
include("billiards/lemon.jl")
include("billiards/sinai.jl")
include("billiards/triangle.jl")
export adapt_basis
#include("limacon.jl")
#include("rectangle.jl")
export Stadium, Lemon, Triangle, Sinai
export curve, tangent, normal, arc_length
export tangent_vec, normal_vec
#convenience functions may be moved somewhere else
#export make_stadium_and_basis, make_triangle_and_basis 



#solvers
#include("solvers/Solvers.jl")
#@reexport using .Solvers
include("solvers/samplers.jl")
export linear_nodes, gauss_legendre_nodes, fourier_nodes
include("solvers/decompositions.jl")
include("solvers/matrixconstructors.jl")
export basis_matrix, basis_and_gradient_matrices, dk_matrix

include("solvers/acceleratedmethods/acceleratedmethods.jl")
include("solvers/sweepmethods/sweepmethods.jl")
export ScalingMethodA, ScalingMethodB 
export DecompositionMethod, ParticularSolutionsMethod
export BoundaryPointsSM, BoundaryPointsDM
export evaluate_points, construct_matrices, construct_matrices_benchmark
export solve, solve_vect
export solve_wavenumber, solve_spectrum
export k_sweep

#spectra
#include("spectra/Spectra.jl")
#@reexport using .Spectra

include("spectra/spectralutils.jl")
include("spectra/unfolding.jl")
export weyl_law
export compute_spectrum
#states
#include("states/States.jl")
#@reexport using .States
include("states/eigenstates.jl")
include("states/basisstates.jl")
include("states/randomstates.jl")

export Eigenstate, EigenstateBundle, BasisState, GaussianRandomState
export compute_eigenstate, compute_eigenstate_bundle

include("states/wavefunctions.jl")
include("states/boundaryfunctions.jl")
include("states/husimifunctions.jl")

export wavefunction #wavefunction_norm 
export boundary_function, momentum_function, husimi_function
export boundary_coords
#plotting functions in Makie
#include("plotting/Plotting.jl")
#@reexport using .Plotting
include("plotting/plottingmakie.jl")
export plot_curve!, plot_boundary!
export plot_domain_fun!, plot_domain!
export plot_lattice!
export plot_wavefunction!, plot_wavefunction_gradient!, plot_probability!
export plot_boundary_function!, plot_momentum_function!, plot_husimi_function!
export plot_heatmap!, plot_heatmap_balaced!
include("plotting/testplotting.jl")
export  plot_geometry_test!, plot_basis_test!, plot_solver_test!, plot_state_test!, plot_matrix!

include("utils/benchmarkutils.jl")
export BenchmarkInfo
export benchmark_solver, compute_benchmarks

include("plotting/benchmarkplotting.jl")
export plot_benchmarks!

#convenience functions
function make_stadium_and_basis(half_width;radius=1.0,x0=zero(half_width),y0=zero(half_width), rot_angle=zero(half_width))
    billiard = Stadium(half_width; radius=radius,x0=x0,y0=y0)
    basis = CornerAdaptedFourierBessel(1, pi/2.0, SVector(x0,y0),rot_angle) 
    return billiard, basis 
end
export make_stadium_and_basis

function make_lemon_and_basis(half_separation;full_domain=false,radius=1.0,x0=zero(half_separation),y0=zero(half_separation), rot_angle=zero(half_separation))
    billiard = Lemon(half_separation; radius=radius,x0=x0,y0=y0)
    basis = CornerAdaptedFourierBessel(1, pi/2.0, SVector(x0,y0),rot_angle) 
    return billiard, basis 
end
export make_lemon_and_basis

function make_triangle_and_basis(gamma,chi; edge_i=1)
    cor = Triangle(gamma,chi).corners
    x0,y0 = cor[mod1(edge_i+2,3)]
    re = [:Virtua, :Virtual, :Virtual]
    re[edge_i] = :Real 
    tr = Triangle(gamma,chi; curve_types = re, x0 = x0, y0 =y0)
    #dim = round(Int, tr.boundary[edge_i].length*k*solver.dim_scaling_factor/(2*pi))
    basis = CornerAdaptedFourierBessel(1, adapt_basis(tr,edge_i+2)...) 
    return tr, basis 
end
export make_triangle_and_basis

end
