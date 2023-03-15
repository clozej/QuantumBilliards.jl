using MKL
using Revise
include("../src/QuantumBilliards.jl")

using .QuantumBilliards
#using Revise 
using CairoMakie
#using Latexify


gamma = pi/2 #sqrt(2)/2 * pi
chi  = 4.5#sqrt(5)

billiard1, basis1 = make_triangle_and_basis(gamma, chi; edge_i =1)
billiard2, basis2 = make_triangle_and_basis(gamma, chi; edge_i =2)
billiard3, basis3 = make_triangle_and_basis(gamma, chi; edge_i =3)

f = Figure(resolution = (1000,500))
ax = Axis(f[1,1])
plot_geometry_test!(ax, billiard2)
display(f)

f = Figure(resolution = (1000,500))
plot_basis_test!(f, basis2, billiard2)
display(f)

d = 3.0
b = 10.0
sw_solver = DecompositionMethod(d,b)
acc_solverA = ScalingMethodA(d,b)
acc_solverB = ScalingMethodB(d,b)
acc_solver = acc_solverA

k_range = collect(range(5.0,10.0,step=0.01))
ten1 = k_sweep(sw_solver,basis1,billiard1,k_range)
ten2 = k_sweep(sw_solver,basis2,billiard2,k_range)
ten3 = k_sweep(sw_solver,basis3,billiard3,k_range)

f = Figure(resolution = (1000,500))
ax = Axis(f[1,1])
lines!(ax,k_range,log10.(ten1))
lines!(ax,k_range,log10.(ten2))
lines!(ax,k_range,log10.(ten3))
display(f)


k0 = 500.00
dk = 0.05
acc_infoA = benchmark_solver(acc_solverA, basis1, billiard1, gauss_legendre_nodes, k0, dk; plot_matrix=true);
acc_infoB = benchmark_solver(acc_solverB, basis1, billiard1, gauss_legendre_nodes, k0, dk; plot_matrix=true);

acc_infoA = benchmark_solver(acc_solverA, basis2, billiard2, gauss_legendre_nodes, k0, dk; plot_matrix=true);
acc_infoB = benchmark_solver(acc_solverB, basis2, billiard2, gauss_legendre_nodes, k0, dk; plot_matrix=true);

acc_infoA = benchmark_solver(acc_solverA, basis3, billiard3, gauss_legendre_nodes, k0, dk; plot_matrix=true);
acc_infoB = benchmark_solver(acc_solverB, basis3, billiard3, gauss_legendre_nodes, k0, dk; plot_matrix=true);


k0 = 1000.001
dk = 0.001
acc_info = benchmark_solver(acc_solver, basis, billiard, gauss_legendre_nodes, k0, dk; plot_matrix=true);
sw_info = benchmark_solver(sw_solver, basis, billiard, gauss_legendre_nodes, k0, dk; plot_matrix=true, log=true);

k1a, ten1a = solve_spectrum(acc_solverA,basis1, billiard1,k0,dk)
k1b, ten1b = solve_spectrum(acc_solverB,basis1, billiard1,k0,dk)
k2a, ten2a = solve_spectrum(acc_solverA,basis2, billiard2,k0,dk)
k2b, ten2b = solve_spectrum(acc_solverB,basis2, billiard2,k0,dk)
k3a, ten3a = solve_spectrum(acc_solverA,basis3, billiard3,k0,dk)
k3b, ten3b = solve_spectrum(acc_solverB,basis3, billiard3,k0,dk)

f = Figure(resolution = (1000,500))
ax = Axis(f[1,1])
scatter!(ax,k1a,log10.(ten1a))
scatter!(ax,k1b,log10.(ten1b))
scatter!(ax,k2a,log10.(ten2a))
scatter!(ax,k2b,log10.(ten2b))
scatter!(ax,k3a,log10.(ten3a))
scatter!(ax,k3b,log10.(ten3b))
display(f)

f = Figure(resolution = (1000,500))
plot_solver_test!(f,sw_solver,basis,billiard,5.0,10.0,0.01)
display(f)

f = Figure(resolution = (1000,500))
plot_solver_test!(f,acc_solverA,basis,billiard,100.0,101.0,0.1)
display(f)

f = Figure(resolution = (1000,500))
plot_solver_test!(f,acc_solverB,basis,billiard,100.0,101.0,0.1)
display(f)

k0 = 500.10
dk = 0.01
k, ten = solve_wavenumber(acc_solver,basis, billiard,k0,dk)
k, ten = solve_spectrum(acc_solver,basis, billiard,k0,dk)
k, ten = solve_wavenumber(sw_solver,basis, billiard,k0,dk)
state = compute_eigenstate(sw_solver, basis, billiard, k)
state1 = compute_eigenstate(acc_solver, basis, billiard, k0)

states = compute_eigenstate_bundle(acc_solver, basis, billiard, k0;dk =0.1, tol=0.0005)
states.ks
states.tens


f = Figure(resolution = (1500,1500))
plot_husimi_function!(f,states,basis,billiard)
display(f)

f = Figure(resolution = (1500,1500))
plot_boundary_function!(f,states,basis,billiard)
display(f)

f = Figure(resolution = (1500,1500))
plot_momentum_function!(f,states,basis,billiard;log=true)
display(f)

f = Figure(resolution = (1500,1500))
plot_probability!(f,states,basis,billiard)
display(f)

f = Figure(resolution = (1500,1500))
plot_state_test!(f,state,basis,billiard)
display(f)

f = Figure(resolution = (1500,1500))
plot_state_test!(f,state1,basis,billiard)
display(f)

b_range =collect(range(2.0,6.0,step=0.5))
f = Figure(resolution = (1000,500))
plot_benchmarks!(f, sw_solver, basis, billiard, gauss_legendre_nodes, k0, dk, 3.5, b_range)
display(f)

