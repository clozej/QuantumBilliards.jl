#using  CoordinateTransformations, Rotations, StaticArrays
using GLMakie
using Revise, BenchmarkTools
#include("../src/billiards/coordinatesystems.jl")
#include("../src/billiards/geometry.jl")
include("../src/billiards/triangle.jl")
include("../src/basis/fourierbessel/corneradapted.jl")
include("../src/solvers/acceleratedmethods/acceleratedmethods.jl")
include("../src/solvers/acceleratedmethods/scalingmethod.jl")
include("../src/solvers/sweepmethods/sweepmethods.jl")
include("../src/solvers/sweepmethods/decompositionmethod.jl")
include("../src/solvers/matrixconstructors.jl")
include("../src/states/eigenstates.jl")
include("../src/states/wavefunctions.jl")
include("../src/plotting/plottingmakie.jl")

gamma = sqrt(2)/2 * pi
chi  = 2.0

billiard, basis = make_triangle_and_basis(gamma, chi; edge_i =3)
basis = rescale_basis(basis, 24)
#basis.cs
d = 5.0
b = 5.0
acc_solver = ScalingMethod(d,b)
sw_solver = DecompositionMethod(d,b)

xlim,ylim = boundary_limits(billiard.boundary; grd=1000, type=typeof(d))
nx = 512
ny = 512
x_grid = range(xlim... , nx)
y_grid = range(ylim... , ny)
pts = [SVector(x,y) for y in y_grid for x in x_grid]
pts_mask = is_inside(billiard,x_grid,y_grid)
pts = pts[pts_mask]
pts
pts_pol = collect(cartesian_to_polar(pt) for pt in pts)

function test_basis(i::Int,k,x_grid,y_grid)
    let i = i, k = k, x_grid=x_grid, y_grid=y_grid
    return collect(sin(i*k*x)*cos(i*k*y) for y in y_grid for x in x_grid)
    end
end

#@code_warntype basis_matrix(basis,10.0, pts, 1:5)
#=
function test_basis(i::Int,k,x_grid,y_grid)
    [sin(k*x)*cos(k*y) for y in y_grid for x in x_grid]
end
=#
println("Basis function")
@btime basis_fun(basis,1,10.0,pts)
@btime basis_fun(basis,1,10.0,x_grid,y_grid)
@btime basis_fun(basis,1:24,10.0,pts)
@btime basis_fun(basis,1:24,10.0,x_grid,y_grid)
println("Basis matrix")
@btime basis_matrix(basis,10.0, pts, 1:8)
@btime basis_matrix(basis,10.0, pts)
@btime basis_matrix(basis,10.0, x_grid,y_grid, 1:8)
@btime basis_matrix(basis,10.0, x_grid,y_grid)
println("dB_dk")
@btime dk_fun(basis,1,10.0,pts)
@btime dk_fun(basis,1:24,10.0,pts)
println("dB_dk matrix")
@btime dk_matrix(basis,10.0,pts, 1:8)
@btime dk_matrix(basis,10.0,pts) 
println("Gradient")
@btime gradient(basis,1,10.0,pts)
@btime gradient(basis,1:24,10.0,pts)
println("Gradient matrices")
@btime gradient_matrices(basis,10.0,pts, 1:8)
@btime gradient_matrices(basis,10.0,pts)
println("Basis and gradient")
@btime basis_and_gradient(basis,1,10.0,pts)
@btime basis_and_gradient(basis,1:24,10.0,pts) 
println("Basis and gradient matrices")
@btime basis_and_gradient_matrices(basis,10.0,pts, 1:8)
@btime basis_and_gradient_matrices(basis,10.0,pts)
println("plane waves test")
@btime test_basis(1,10.0,x_grid,y_grid);


B = basis_fun(basis,1:10,10.0,x_grid,y_grid)
Z = reshape(B[:,2], nx,ny)
f = Figure(resolution = (1000,1000));
hmap, ax = plot_heatmap_balaced!(f,x_grid,y_grid,Z)
plot_boundary!(ax,billiard; dens= 25.0)
display(f)

#@code_warntype DecompositionMethod(d,b)
#@code_warntype ScalingMethod(d,b)
