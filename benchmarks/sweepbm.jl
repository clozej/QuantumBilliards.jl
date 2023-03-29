using MKL
#using Revise
include("../src/QuantumBilliards.jl")

using .QuantumBilliards
#using Revise 
using CairoMakie
#using Latexify

gamma = 2/3*pi #sqrt(2)/2 * pi
chi  = 2.0

billiard, basis = make_triangle_and_basis(gamma, chi)

f = Figure(resolution = (1000,500))
plot_basis_test!(f, basis, billiard)
display(f)

d = 5.0
b = 10.0
b_int = b
dm_solver = DecompositionMethod(d,b)
psm_solver = ParticularSolutionsMethod(d,b,b)

k0 = 7.2
dk = 0.2
k_dm, ten_dm = solve_wavenumber(dm_solver,basis, billiard, k0, dk)
k_psm, ten_psm = solve_wavenumber(psm_solver,basis, billiard, k0, dk)
#k, ten = solve_wavenumber(sw_solver,basis, billiard,k0,dk)
k_range = collect(range(5.0,10.0,step=0.005))
tens_dm = k_sweep(dm_solver,basis,billiard,k_range)
tens_psm = k_sweep(psm_solver,basis,billiard,k_range)
#tens_sym = k_sweep(sw_solver,basis_sym,billiard,k_range)

f = Figure(resolution = (1000,500))
ax = Axis(f[1,1])
lines!(ax,k_range,log10.(tens_dm))
lines!(ax,k_range,log10.(tens_psm))
scatter!(ax,k_dm,log10.(ten_dm))
scatter!(ax,k_psm,log10.(ten_psm))
vlines!(ax, [k0]; color=:black, linewidth=0.5)
display(f)

state_dm = compute_eigenstate(dm_solver, basis, billiard, k_dm)
state_psm = compute_eigenstate(psm_solver, basis, billiard, k_psm)
k_dm
k_psm
state_psm.vec
state_dm.vec

f = Figure(resolution = (1500,1500))
plot_probability!(f[1,1],state_dm;log=false,inside_only=false)
plot_probability!(f[2,1],state_dm;log=(true,-5),inside_only=false)
display(f)

f = Figure(resolution = (1500,1500))
plot_probability!(f[1,1],state_psm;log=false,inside_only=false)
plot_probability!(f[2,1],state_psm;log=(true,-5),inside_only=false)
display(f)

sw_info = benchmark_solver(dm_solver, basis, billiard, gauss_legendre_nodes, k0, dk; plot_matrix=true, log=false);
sw_info = benchmark_solver(psm_solver, basis, billiard, gauss_legendre_nodes, k0, dk; plot_matrix=true, log=false);


using LinearAlgebra
B = rand(10,10)
B_int = rand(5,10)
F = svd(B, B_int)
#H = F.R*F.Q'
idx = 1:F.k + F.l #inidices containing the singular values we need
sv = F.a[idx] ./ F.b[idx] #generalized singular values
U = F.U[:,idx]
U*
i_min = argmin(sv)
F = svd(B, B_int)
#H = F.R*F.Q'
idx = 1:F.k + F.l #inidices containing the singular values we need
sv = F.a[idx] ./ F.b[idx] #generalized singular values
U = F.U[:,idx]
i_min = argmin(sv)
H = F.R*F.Q'
F.Q'*F.Q
X = F.U*F.D1
X