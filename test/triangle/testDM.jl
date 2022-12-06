using MKL
using GLMakie
using Revise, BenchmarkTools
include("../../src/QuantumBilliards.jl")
import .QuantumBilliards as qb

#Strong-Mixing
gamma = pi/2*(1+sqrt(5))/2
chi  = (1+sqrt(2))/2
edge_i=1

#Weak-Mixing
gamma = sqrt(2)/2 * pi
chi  = 2.0
edge_i=1

#Veech
kind =7
n = 4 + kind
gamma = pi/2
chi = (n-2)/2
edge_i=3

#solver parameters
d = 4.0
b = 4.0
b_int = 4.0
acc_solver = qb.ScalingMethod(d,b) #accelerated solver
sw_solver = qb.DecompositionMethod(d,b) #sweep solver

#typeof(sw_solver)<:qb.SweepSolver

billiard, basis = qb.make_triangle_and_basis(gamma,chi; edge_i=edge_i)

ks = collect(LinRange(4.0, 10.0, 1000))
ts = qb.k_sweep(sw_solver,basis,billiard,ks)

f = Figure(resolution = (1000,500));
#axis = [Axis(f[1,i]) for i in 1:3]
ax = Axis(f[1,1],xlabel=L"k", ylabel=L"\log(t)")
lines!(ax, ks, log10.(ts))
#lines!(ax, ks, log10.(ten2))
#lines!(ax, ks, ten)
display(f)

k0 = 6.6
dk = 0.1
k, ten = qb.solve_wavenumber(sw_solver,basis,billiard,k0,dk)
state = qb.compute_eigenstate(sw_solver, basis, billiard, k)
state32 = qb.Eigenstate(Float32(state.k), Float32.(state.vec))  

#eltype(k) == Float64
Psi, x, y = qb.wavefunction(state,basis,billiard, inside_only=true)

k0 = 300.8
dk = 0.05
k, ten = qb.solve_wavenumber(acc_solver,basis,billiard,k0,dk)
ks_sm, ts_sm = qb.solve_spectrum(acc_solver,basis,billiard,k0,dk)



f = Figure(resolution = (1000,1000));
qb.plot_wavefunction!(f,state, basis, billiard) 
display(f)

f = Figure(resolution = (1000,1000));
qb.plot_probability!(f,state32, basis, billiard, b=25.0, log =(true,-5), plot_normal=false, inside_only=true) 
display(f)

#work in progress

s, u, U = qb.boundary_function(state, basis, billiard;include_virtual=true, sampler=qb.gauss_legendre_nodes)
f = Figure(resolution = (1000,1000));
#axis = [Axis(f[1,i]) for i in 1:3]
ax = Axis(f[1,1],xlabel=L"s", ylabel=L"u")
lines!(ax,s, u)
display(f)
s

bs = Basisstate(10.0,3,100)
f = Figure(resolution = (1000,1000));
plot_wavefunction!(f,bs, basis, billiard) 
display(f)


s, u = boundary_function(bs, basis, billiard)
f = Figure(resolution = (1000,1000));
#axis = [Axis(f[1,i]) for i in 1:3]
ax = Axis(f[1,1],xlabel=L"s", ylabel=L"u")
lines!(ax,s, u)
display(f)

