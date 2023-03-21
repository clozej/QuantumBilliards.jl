#using MKL
include("../src/QuantumBilliards.jl")
#using Revise
using .QuantumBilliards
#using Revise 
using CairoMakie
using StaticArrays
#using Latexify

angle1 = 0.7
angle2 = 0.4

billiard = Sinai(angle1,angle2;full_domain=false)

f = Figure(resolution = (1000,500))
ax = Axis(f[1,1])
plot_geometry_test!(ax, billiard)
display(f)
odd_x_odd_y = [XReflection(-1), YReflection(-1) , XYReflection(-1,-1)]
basis = FundamentalBessel(odd_x_odd_y)
k = 50.0
basis = resize_basis(basis,billiard,100,k)


f = Figure(resolution = (1000,500))
plot_basis_test!(f, basis, billiard, i=50)
display(f)

d = 5.0
b = 10.0
sw_solver = DecompositionMethod(d,b)

k0 = 100.15
dk = 0.05
k, ten = solve_wavenumber(sw_solver,basis, billiard,k0,dk)
#k, ten = solve_wavenumber(sw_solver,basis, billiard,k0,dk)
k_range = collect(range(2.0,10.0,step=0.01))
tens = k_sweep(sw_solver,basis,billiard,k_range)

f = Figure(resolution = (1000,500))
ax = Axis(f[1,1])
lines!(ax,k_range,log10.(tens))
#lines!(ax,k_range,log10.(tens2))
#vlines!(ax, [k]; color=:black, linewidth=0.5)
display(f)
#println((k,ten))


k0 = 100.15
dk = 0.05
k, ten = solve_wavenumber(sw_solver,basis, billiard,k0,dk)
#k, ten = solve_wavenumber(sw_solver,basis, billiard,k0,dk)
k_range = collect(range(100.0,100.5,step=0.005))
tens = k_sweep(sw_solver,basis,billiard,k_range)

f = Figure(resolution = (1000,500))
ax = Axis(f[1,1])
lines!(ax,k_range,log10.(tens))
#lines!(ax,k_range,log10.(tens2))
vlines!(ax, [k]; color=:black, linewidth=0.5)
display(f)
println((k,ten))

state = compute_eigenstate(sw_solver, basis, billiard, k)
f = Figure(resolution = (1500,1500))
plot_probability!(f,state;log=(true,-5),inside_only=false)
display(f)
state.basis.symmetries
f = Figure(resolution = (1500,1500))
plot_probability!(f,state;log=false,inside_only=true)
display(f)

k0 = 500.15
dk = 0.05
k, ten = solve_wavenumber(sw_solver,basis, billiard,k0,dk)
#k, ten = solve_wavenumber(sw_solver,basis, billiard,k0,dk)
k_range = collect(range(500.0,500.2,step=0.001))
tens = k_sweep(sw_solver,basis,billiard,k_range)

f = Figure(resolution = (1000,500))
ax = Axis(f[1,1])
lines!(ax,k_range,log10.(tens))
#lines!(ax,k_range,log10.(tens2))
vlines!(ax, [k]; color=:black, linewidth=0.5)
display(f)
println((k,ten))

state = compute_eigenstate(sw_solver, basis, billiard, k)
f = Figure(resolution = (1500,1500))
plot_probability!(f,state;log=(true,-5),inside_only=false)
display(f)

f = Figure(resolution = (1500,1500))
plot_probability!(f,state;log=false,inside_only=true)
display(f)