using MKL
include("../src/QuantumBilliards.jl")
#using Revise
using .QuantumBilliards
#using Revise 
using CairoMakie
#using Latexify

B = 0.42 #sqrt(2)/2 * pi
billiard, basis = make_lemon_and_basis(B;full_domain=false)

f = Figure(resolution = (1000,500))
ax = Axis(f[1,1])
plot_geometry_test!(ax, billiard)
display(f)

f = Figure(resolution = (1000,500))
plot_basis_test!(f, basis, billiard)
display(f)

d = 3.0
b = 5.0
sw_solver = DecompositionMethod(d,b)
acc_solverA = ScalingMethodA(d,b)
acc_solverB = ScalingMethodB(d,b)
acc_solver = acc_solverA

k0 = 500.00
dk = 0.1
acc_infoA = benchmark_solver(acc_solverA, basis, billiard, gauss_legendre_nodes, k0, dk; plot_matrix=true);
acc_infoB = benchmark_solver(acc_solverB, basis, billiard, gauss_legendre_nodes, k0, dk; plot_matrix=true);

sw_info = benchmark_solver(sw_solver, basis, billiard, gauss_legendre_nodes, k0, dk; plot_matrix=true, log=false);


f = Figure(resolution = (1000,500))
plot_solver_test!(f,sw_solver,basis,billiard,5.0,10.0,0.01)
display(f)

f = Figure(resolution = (1000,500))
plot_solver_test!(f,acc_solverA,basis,billiard,100.0,101.0,0.25)
display(f)

f = Figure(resolution = (1000,500))
plot_solver_test!(f,acc_solverB,basis,billiard,100.0,101.0,0.25)
display(f)

f = Figure(resolution = (1000,500))
plot_solver_test!(f,acc_solver,basis,billiard,500.0,501.0,0.125, tol = 1e-3)
display(f)

k0 = 1000.0
dk = 0.1
k, ten = solve_wavenumber(acc_solver,basis, billiard,k0,dk)
ks, ten = solve_spectrum(acc_solver,basis, billiard,k0,dk)
k, ten = solve_wavenumber(sw_solver,basis, billiard,k0,dk)
state = compute_eigenstate(sw_solver, basis, billiard, k)

states = compute_eigenstate_bundle(acc_solver, basis, billiard, k0;dk =0.1, tol=0.0025)
states.X
states.ks
states.tens

f = Figure(resolution = (800,2500))
plot_probability!(f,states,b=10.0, log=(true,-5))
display(f)

f = Figure(resolution = (1500,1500))
plot_momentum_function!(f,states;log=false)
display(f)

f = Figure(resolution = (1500,1500))
plot_husimi_function!(f,states;log=false,include_virtual=false)
display(f)

f = Figure(resolution = (1500,1500))
plot_state_test!(f,state; b_u= 10.0, include_virtual=false)
display(f)

#=
using LinearAlgebra, StaticArrays
using TimerOutputs

function construct_matrices_benchmark1(solver::ScalingMethod, basis::Ba, pts::BoundaryPointsSM, k) where {Ba<:AbsBasis}
    to = TimerOutput()
    #type = eltype(pts.w)
    xy, w = pts.xy, pts.w
    #M =  length(xy)
    #basis and gradient matrices
    @timeit to "basis_and_gradient_matrices" B, dX, dY = basis_and_gradient_matrices(basis, k, xy)
    N = basis.dim
    type = eltype(B)
    F = zeros(type,(N,N))
    Fk = similar(F)
    @timeit to "F construction" begin 
        @timeit to "weights" T = (w.*B) #reused later
        #@timeit to "copy" Bt = copy(B')
        @timeit to "product" mul!(F,B',T) #boundary norm matrix
    end
    #reuse B
    @timeit to "Fk construction" begin 
        r = xy# normalize.(xy)
        @timeit to "dilation derivative" x = getindex.(r,1)
        @timeit to "dilation derivative" y = getindex.(r,2)
        #inplace modifications
        @timeit to "dilation derivative" dX = x .* dX 
        @timeit to "dilation derivative" dY = y .* dY
        #reuse B
        @timeit to "dilation derivative" B = dX .+ dY
        @timeit to "product" mul!(Fk,B',T) #B is now derivative matrix
        #symmetrize matrix
        @timeit to "addition" Fk = Fk+Fk'
        Fk = Fk ./ k
    end
    print_timer(to)    
    return F, Fk        
end

function construct_matrices_benchmark2(solver::ScalingMethod, basis::Ba, pts::BoundaryPointsSM, k) where {Ba<:AbsBasis}
    to = TimerOutput()
    #type = eltype(pts.w)
    xy, w = pts.xy, pts.w
    #M =  length(xy)
    N = basis.dim
    #basis matrix
    @timeit to "basis_matrix" B = basis_matrix(basis, k, xy)
    type = eltype(B)
    F = zeros(type,(N,N))
    Fk = similar(F)
    @timeit to "F construction" begin 
        @timeit to "weights" T = (w.*B) #reused later
        #@timeit to "copy" Bt = copy(B')
        @timeit to "product" mul!(F,B',T) #boundary norm matrix
    end
    #reuse B
    @timeit to "dk_matrix" B = dk_matrix(basis,k, xy)
    @timeit to "Fk construction" begin 
        @timeit to "product" mul!(Fk,B',T) #B is now derivative matrix
        #symmetrize matrix
        @timeit to "addition" Fk = Fk+Fk'
    end
    print_timer(to)    
    return F, Fk        
end
=#

#=
k0 = 500.10
dk = 0.01
new_basis = resize_basis(basis,200)
pts = evaluate_points(acc_solver, billiard, gauss_legendre_nodes,k0)

F1, Fk1 = construct_matrices_benchmark1(acc_solver, new_basis, pts, k0)
F2, Fk2 = construct_matrices_benchmark2(acc_solver, new_basis, pts, k0)



f = Figure(resolution = (1000,2000))
plot_matrix!(f[1,1], F1)
plot_matrix!(f[1,2], F2)
display(f)

f = Figure(resolution = (1000,2000))
plot_matrix!(f[1,1], Fk1)
plot_matrix!(f[1,2], Fk2)
display(f)

=#