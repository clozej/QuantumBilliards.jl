using MKL
using Revise
include("../src/QuantumBilliardsTest.jl")
using CairoMakie

gamma = sqrt(2)/2 * pi #2/3*pi 
chi  = 2.0

billiard, basis = make_triangle_and_basis(gamma, chi)

d = 1.6
b = 5.0
sw_solver = DecompositionMethod(d,b)
acc_solver = ScalingMethodA(d,b)

k0 = 1000.00
dk = 0.1
acc_infoA = benchmark_solver(acc_solver, basis, billiard, k0, dk; plot_matrix=true);

state = 1000
n = 30
k0 = k_at_state(state, billiard.area, billiard.length, billiard.angles)
k1 = k_at_state(state+n, billiard.area, billiard.length, billiard.angles)
n_expected = 2.0
dk = n_expected * 2.0*pi / (billiard.area * k1)

f = Figure(resolution = (1000,500))
plot_solver_test!(f,acc_solver,basis,billiard,k0,k1,dk)
display(f)
