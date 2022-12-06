using MKL
using GLMakie
using Revise, BenchmarkTools
include("../src/QuantumBilliards.jl")
import .QuantumBilliards as qb

#solver parameters
d = 10.0
b = 10.0
b_int = 10.0
acc_solver = qb.ScalingMethod(d,b) #accelerated solver
sw_solver = qb.DecompositionMethod(d,b) #sweep solver
sw_solver = qb.ParticularSolutionsMethod(d,b,b_int)

chi  = (1+sqrt(2))/2
billiard, basis = qb.make_rectangle_and_basis(chi)

ks = collect(LinRange(0.1, 10.0, 1000))
ts = qb.k_sweep(sw_solver,basis,billiard,ks)


f = Figure(resolution = (1000,500));
#axis = [Axis(f[1,i]) for i in 1:3]
ax = Axis(f[1,1],xlabel=L"k", ylabel=L"\log(t)")
lines!(ax, ks, log10.(ts))
#lines!(ax, ks, log10.(ten2))
#lines!(ax, ks, ten)
display(f)


k0 = 4.6
dk = 0.5
k, ten = qb.solve_wavenumber(sw_solver,basis,billiard,k0,dk)
state = qb.compute_eigenstate(sw_solver, basis, billiard, k)
state32 = qb.Eigenstate(Float32(state.k), Float32.(state.vec))  


k0 = 300.8
dk = 0.05
k, ten = qb.solve_wavenumber(acc_solver,basis,billiard,k0,dk)
ks_sm, ts_sm = qb.solve_spectrum(acc_solver,basis,billiard,k0,dk)



f = Figure(resolution = (1000,1000));
qb.plot_probability!(f,state32, basis, billiard, b=25.0, log =(false,-5), plot_normal=false, inside_only=true) 
display(f)
