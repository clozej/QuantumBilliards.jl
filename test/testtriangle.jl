#using  CoordinateTransformations, Rotations, StaticArrays
using GLMakie
#using Revise
#include("../src/billiards/coordinatesystems.jl")
#include("../src/billiards/geometry.jl")
include("../src/billiards/triangle.jl")
include("../src/solvers/acceleratedmethods/acceleratedmethods.jl")
include("../src/solvers/acceleratedmethods/scalingmethod.jl")
include("../src/solvers/sweepmethods/sweepmethods.jl")
include("../src/solvers/sweepmethods/decompositionmethod.jl")
include("../src/states/eigenstates.jl")
include("../src/states/basisstates.jl")
include("../src/states/wavefunctions.jl")
include("../src/states/gradients.jl")
include("../src/states/boundaryfunctions.jl")
include("../src/plotting/plottingmakie.jl")
include("../src/utils/benchmarkutils.jl")
#include("../src/utils/gridutils.jl")
#testing all curve types from geometry.jl


gamma = 2/3*pi #sqrt(2)/2 * pi
chi  = 2.0

billiard, basis = make_triangle_and_basis(gamma, chi)
basis32 = Float32(basis)
d = 5.0
b = 10.0
acc_solver = ScalingMethod(d,b)
sw_solver = DecompositionMethod(d,b)

ks = collect(LinRange(3.0, 10.0, 1000))
ts = k_sweep(sw_solver,basis,billiard,ks)
f = Figure(resolution = (1000,500));
#axis = [Axis(f[1,i]) for i in 1:3]
ax = Axis(f[1,1],xlabel=L"k", ylabel=L"\log(t)")
lines!(ax, ks, log10.(ts))
#lines!(ax, ks, log10.(ten2))
#lines!(ax, ks, ten)
display(f)

k0 = 5.5
dk = 0.5
k, ten = solve_wavenumber(sw_solver,basis,billiard,k0,dk)
state = compute_eigenstate(sw_solver, basis, billiard, k)
state32 = Eigenstate(Float32(state.k), Float32.(state.vec))


k0 = 1005.8
dk = 0.05
k, ten = solve_wavenumber(acc_solver,basis,billiard,k0,dk)
state = compute_eigenstate(sw_solver, basis, billiard, k)
state32 = Eigenstate(Float32(state.k), Float32.(state.vec))


f = Figure(resolution = (1000,1000));
plot_probability!(f,state32, basis32, billiard; b=10.0, log =(true,-5), plot_normal=false, inside_only=true) 
display(f)

f = Figure(resolution = (1000,500));
ax = Axis(f[1,1],xlabel=L"s", ylabel=L"u")
plot_boundary_function!(ax,state, basis, billiard; b=50.0) 
display(f)


f = Figure(resolution = (1000,1000));
plot_wavefunction!(f,state32, basis32, billiard, b=10.0, plot_normal=false, inside_only=true) 
display(f)

f = Figure(resolution = (1000,1000));
plot_wavefunction_gradient!(f,state32, basis32, billiard; b=100.0, lengthscale= 0.0025, plot_normal=false, inside_only=true) 
display(f)

basisstate = BasisState(10.0, 1, 10)
dX, dY, x_grid, y_grid = wavefunction_gradient(basisstate, basis, b=200.0,  inside_only=false)
f = Figure(resolution = (2000,1000));
hmap, ax = plot_heatmap_balaced!(f[1,1],x_grid, y_grid,dX )
hmap, ax = plot_heatmap_balaced!(f[1,2],x_grid, y_grid,dY )
display(f)



f = Figure(resolution = (1000,1000));
plot_wavefunction!(f,basisstate, basis, billiard, b=10.0, plot_normal=false, inside_only=false) 
display(f)

f = Figure(resolution = (1000,500));
ax = Axis(f[1,1],xlabel=L"s", ylabel=L"u")
plot_boundary_function!(ax,basisstate, basis, billiard; b=10.0) 
display(f)


#=

pts = boundary_coords(billiard, 20)
pts.normal

new_basis = rescale_basis(basis, 10)
dX, dY = gradient_matrices(new_basis, 10.0, pts.xy)
dX1, dY1 = gradient(new_basis, 1:10, 10.0, pts.xy)
dx, dy = gradient(new_basis, 1, 10.0, pts.xy)
nx = getindex.(pts.normal,1)
ny = getindex.(pts.normal,2)

dX[:,1] - dx
dX[:,1] - dX1[:,1]

dY[:,1] - dy
dY[:,1] - dY1[:,1]

dY[:,1] - dy
dY


=#