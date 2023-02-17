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
include("../src/states/wavefunctions.jl")
include("../src/plotting/plottingmakie.jl")
include("../src/utils/benchmarkutils.jl")
#include("../src/utils/gridutils.jl")
#testing all curve types from geometry.jl


gamma = sqrt(2)/2 * pi
chi  = 2.0

billiard, basis = make_triangle_and_basis(gamma, chi)
basis32 = Float32(basis)
d = 5.0
b = 10.0
acc_solver = ScalingMethod(d,b)
sw_solver = DecompositionMethod(d,b)
#=
ks = collect(LinRange(3.0, 10.0, 1000))
ts = k_sweep(sw_solver,basis,billiard,ks)
f = Figure(resolution = (1000,500));
#axis = [Axis(f[1,i]) for i in 1:3]
ax = Axis(f[1,1],xlabel=L"k", ylabel=L"\log(t)")
lines!(ax, ks, log10.(ts))
#lines!(ax, ks, log10.(ten2))
#lines!(ax, ks, ten)
display(f)
=#

k0 = 600.8
dk = 0.05
@time k, ten = solve_wavenumber(acc_solver,basis,billiard,k0,dk)
state = compute_eigenstate(sw_solver, basis, billiard, k)
state32 = Eigenstate(Float32(state.k), Float32.(state.vec))
#=
type = eltype(state.vec)
L = real_length(billiard)
xlim,ylim = boundary_limits(billiard.boundary; grd=round(Int, k*L*b/(2*pi)), type=type)
dy = ylim[2] - ylim[1]
nx = max(round(Int, k*dx*b/(2*pi)), 512)
ny = max(round(Int, k*dy*b/(2*pi)), 512)
x_grid = collect(type,range(xlim... , nx))
y_grid = collect(type,range(ylim... , ny))
new_basis = rescale_basis(basis, state.dim)
#B = basis_matrix(new_basis, k, x_grid, y_grid)
memory_size(B)
=#

@code_warntype compute_psi(state, basis, billiard,x_grid,y_grid)
@code_warntype wavefunction(state,basis,billiard;b=10.0);


@time wavefunction(state,basis,billiard;b=10.0);
@time wavefunction(state32,basis32,billiard;b=10.0);
@time wavefunction(state32,basis32,billiard;b=10.0,inside_only=false);
Psi, x, y = wavefunction(state,basis,billiard;b=10.0)
memory_size(Psi)
Psi
Psi32, x, y = wavefunction(state32,basis32,billiard;b=10.0);
memory_size(Psi32)

P = abs2.(Psi)
sum(P)


inside = reshape((is_inside(billiard, x, y)),length(x),length(y))
outside = .!inside


f = Figure(resolution = (1000,1000));
plot_probability!(f,state, basis, billiard, b=10.0, log =(true,-5), plot_normal=false, inside_only=true) 
display(f)

f = Figure(resolution = (1000,1000));
plot_probability!(f,state32, basis32, billiard, b=10.0, log =(true,-5), plot_normal=false, inside_only=true) 
display(f)
