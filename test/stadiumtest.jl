using MKL
using GLMakie
using Revise, BenchmarkTools
include("../src/QuantumBilliards.jl")
import .QuantumBilliards as qb


#solver parameters
d = 10.0
b = 20.0
b_int = 10.0
acc_solver = qb.ScalingMethod(d,b) #accelerated solver
sw_solver = qb.DecompositionMethod(d,b) #sweep solver
#sw_solver = qb.ParticularSolutionsMethod(d,b,b_int)

eps = 0.5
billiard, basis = qb.make_stadium_and_basis(eps)

f = Figure(resolution = (1000,1000));
#axis = [Axis(f[1,i]) for i in 1:3]
ax = Axis(f[1,1],xlabel=L"x", ylabel=L"y")
qb.plot_billiard!(ax,billiard)
#lines!(ax, ks, ten)
display(f)

ks = collect(LinRange(3.0, 10.0, 1000))
ts = qb.k_sweep(sw_solver,basis,billiard,ks)
ts2 = qb.k_sweep(sw_solver,basis,billiard,ks;sampler=qb.linear_nodes)

f = Figure(resolution = (1000,500));
#axis = [Axis(f[1,i]) for i in 1:3]
ax = Axis(f[1,1],xlabel=L"k", ylabel=L"\log(t)")
lines!(ax, ks, log10.(ts))
lines!(ax, ks, log10.(ts2))
#lines!(ax, ks, log10.(ten2))
#lines!(ax, ks, ten)
display(f)

k0 = 4.3
dk = 1.0
k, ten = qb.solve_wavenumber(sw_solver,basis,billiard,k0,dk)
k, ten = qb.solve_wavenumber(sw_solver,basis,billiard,k0,dk;sampler=qb.linear_nodes)

k0 = 302.6
dk = 0.05
k, ten = qb.solve_wavenumber(acc_solver,basis,billiard,k0,dk)
ks_sm, ts_sm = qb.solve_spectrum(acc_solver,basis,billiard,k0,dk)


state = qb.compute_eigenstate(sw_solver, basis, billiard, k)
state = qb.compute_eigenstate(sw_solver, basis, billiard, k;sampler=qb.linear_nodes)

state32 = qb.Eigenstate(Float32(state.k), Float32.(state.vec))  


f = Figure(resolution = (1000,1000));
qb.plot_probability!(f,state32, basis, billiard, b=25.0, log =(true,-5), plot_normal=false, inside_only=true) 
display(f)

state.vec

f = Figure(resolution = (1000,500));
#axis = [Axis(f[1,i]) for i in 1:3]
ax = Axis(f[1,1],xlabel=L"i", ylabel=L"\log(c)")
lines!(ax, state.vec)
#lines!(ax, ks, log10.(ten2))
#lines!(ax, ks, ten)
display(f)