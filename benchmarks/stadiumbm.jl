using MKL
include("../src/QuantumBilliards.jl")
#using Revise
using .QuantumBilliards
#using Revise 
#using GLMakie
using CairoMakie
#using Latexify



eps = 0.1 #sqrt(2)/2 * pi
billiard, fb_basis = make_stadium_and_basis(eps)
odd_x = [XReflection(-1)]
even_x = [XReflection(1)]
odd_y = [YReflection(-1)]
even_y = [YReflection(1)]
odd_x_odd_y = [XReflection(-1), YReflection(-1) , XYReflection(-1,-1)]
odd_x_even_y = [XReflection(-1), YReflection(1) , XYReflection(-1,1)]
even_x_odd_y = [XReflection(1), YReflection(-1) , XYReflection(1,-1)]
even_x_even_y = [XReflection(1), YReflection(1) , XYReflection(1,1)]
sym_sectors = [odd_x, even_x, odd_y, even_y, odd_x_odd_y,odd_x_even_y,even_x_odd_y, even_x_even_y]
sym_idx = 6
basis = RealPlaneWaves(10,sym_sectors[sym_idx];angle_arc=pi/2.0)

f = Figure(resolution = (1000,1000))

plot_geometry_test!(f, billiard)
display(f)

f = Figure(resolution = (1000,1000))
plot_basis_test!(f[1,1], basis, billiard; i=1)
plot_basis_test!(f[1,2], basis, billiard; i=2)
plot_basis_test!(f[2,1], basis, billiard; i=3)
plot_basis_test!(f[2,2], basis, billiard; i=4)
display(f)


d = 1.5
b = [5.0,10.0]
sw_solver = DecompositionMethod(d,b)
acc_solverA = ScalingMethodA(d,b)
acc_solverB = ScalingMethodB(d,b)
acc_solver = acc_solverA

k0 = 500.00
dk = 0.1
acc_infoA = benchmark_solver(acc_solverA, basis, billiard, k0, dk; plot_matrix=true);
acc_infoB = benchmark_solver(acc_solverB, basis, billiard, k0, dk; plot_matrix=true);

sw_info = benchmark_solver(sw_solver, basis, billiard, k0, dk; plot_matrix=true, log=false);

acc_infoA = benchmark_solver(acc_solverA, basis, billiard, k0, dk; plot_matrix=true);
acc_infoA = benchmark_solver(acc_solverA, fb_basis, billiard, k0, dk; plot_matrix=true);


f = Figure(resolution = (1000,500))
plot_solver_test!(f,sw_solver,basis,billiard,5.0,10.0,0.01)
display(f)

f = Figure(resolution = (1000,500))
plot_solver_test!(f,acc_solverA,basis,billiard,100.0,101.0,0.2)
display(f)

f = Figure(resolution = (1000,500))
plot_solver_test!(f,acc_solverB,basis,billiard,100.0,101.0,0.2)
display(f)

f = Figure(resolution = (1000,500))
plot_solver_test!(f,acc_solver,basis,billiard,500.0,501.0,0.05, tol = 1e-3)
display(f)

k0 = 302.0
dk = 0.1
k_rpw, ten = solve_wavenumber(acc_solver, basis, billiard,k0,dk)
k_fb, ten = solve_wavenumber(acc_solver, fb_basis, billiard,k0,dk)

state_rpw = compute_eigenstate(sw_solver, basis, billiard, k_rpw)
state_fb = compute_eigenstate(sw_solver, fb_basis, billiard, k_fb)

f = Figure(resolution = (1500,1500))
plot_state_test!(f,state_rpw; b_u= 20.0)
display(f)

f = Figure(resolution = (1500,1500))
plot_wavefunction!(f,state_rpw; b= 5.0, dens = 100.0, fundamental_domain=false)
display(f)


f = Figure(resolution = (1500,1500))
plot_state_test!(f,state_fb; b_u= 10.0)
display(f)


ks, ten = solve_spectrum(acc_solver,basis, billiard,k0,dk)
k, ten = solve_wavenumber(sw_solver,basis, billiard,k0,dk)
state = compute_eigenstate(sw_solver, basis, billiard, k)


states = compute_eigenstate_bundle(acc_solver, basis, billiard, k0;dk =0.1, tol=0.02)
states.X
states.ks
states.tens

f = Figure(resolution = (800,2500))
plot_probability!(f,states,b=10.0, log=(true,-5))
display(f)

f = Figure(resolution = (1500,1500))
plot_momentum_function!(f,states;log=true)
display(f)

f = Figure(resolution = (1500,1500))
plot_momentum_function!(f,states;log=false)
display(f)



f = Figure(resolution = (1500,1500))
plot_state_test!(f,state; b_u= 10.0)
display(f)

k0 = 100.0
dk = 0.1
k, ten = solve_wavenumber(acc_solver,basis, billiard,k0,dk)
state = compute_eigenstate(sw_solver, basis, billiard, k)
state1 = compute_eigenstate(acc_solver, basis, billiard, k0)

f = Figure(resolution = (1500,1500))
plot_state_test!(f,state; b_u= 10.0)
display(f)

f = Figure(resolution = (1500,1500))
plot_state_test!(f,state1; b_u= 10.0)
display(f)


f = Figure(resolution = (1000,1000))
plot_wavefunction!(f,state,b=10.0)
display(f)

