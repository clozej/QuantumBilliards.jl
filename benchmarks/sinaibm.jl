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

billiard_full = Sinai(angle1,angle2;full_domain=true)
billiard = Sinai(angle1,angle2;full_domain=false)

f = Figure(resolution = (1000,500))
ax = Axis(f[1,1])
plot_geometry_test!(ax, billiard)
display(f)

basis = FundamentalBessel()
odd_x_odd_y = [XReflection(-1), YReflection(-1) , XYReflection(-1,-1)]
basis_sym = FundamentalBessel(odd_x_odd_y)

k = 25.0
basis = resize_basis(basis,billiard,100,k)


f = Figure(resolution = (1000,500))
plot_basis_test!(f, basis, billiard_full, i=25, k = k)
display(f)

d = 5.0
b = 10.0
dm_solver = DecompositionMethod(d,b)
psm_solver = ParticularSolutionsMethod(d,b,b)

k0 = 3.2
dk = 0.2
k, ten = solve_wavenumber(dm_solver,basis, billiard_full, k0, dk)
k, ten = solve_wavenumber(psm_solver,basis, billiard_full, k0, dk)
#k, ten = solve_wavenumber(sw_solver,basis, billiard,k0,dk)
k_range = collect(range(2.0,10.0,step=0.005))
tens_dm = k_sweep(dm_solver,basis,billiard_full,k_range)
tens_psm = k_sweep(psm_solver,basis,billiard_full,k_range)
#tens_sym = k_sweep(sw_solver,basis_sym,billiard,k_range)

f = Figure(resolution = (1000,500))
ax = Axis(f[1,1])
lines!(ax,k_range,log10.(tens_dm))
lines!(ax,k_range,log10.(tens_psm))
#lines!(ax,k_range,log10.(tens2))
vlines!(ax, [k]; color=:black, linewidth=0.5)
display(f)
#println((k,ten))

state = compute_eigenstate(dm_solver, basis, billiard, k)
state = compute_eigenstate(psm_solver, basis, billiard, k)
f = Figure(resolution = (1500,1500))
plot_probability!(f,state;log=(true,-5),inside_only=false)
display(f)
state.basis.symmetries
f = Figure(resolution = (1500,1500))
plot_probability!(f,state;log=false,inside_only=true)
display(f)


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

rand(SVector{2, Float64}, 5)