include("../abstracttypes.jl")
include("../billiards/billiard.jl")
#include("decompositions.jl")
include("samplers.jl")
include("particularsolutionsmethod.jl")
include("decompositionmethod.jl")
using LinearAlgebra, Optim

function solve_wavenumber(solver::SweepSolver,basis::AbsBasis, billiard::AbsBilliard, k, dk; sampler=gauss_legendre_nodes)
    dim = round(Int, real_length(billiard)*k*solver.dim_scaling_factor/(2*pi))
    new_basis = rescale_basis(basis,dim)
    pts = evaluate_points(solver, billiard, sampler, k)
    function f(k)
        return solve(solver,new_basis,pts,k)
    end
    res =  optimize(f, k-0.5*dk, k+0.5*dk)
    k0,t0 = res.minimizer, res.minimum
    return k0, t0
end

function k_sweep(solver::SweepSolver, basis::AbsBasis, billiard::AbsBilliard, ks; sampler=gauss_legendre_nodes)
    k = maximum(ks)
    dim = round(Int, real_length(billiard)*k*solver.dim_scaling_factor/(2*pi))
    new_basis = rescale_basis(basis,dim)
    pts = evaluate_points(solver, billiard, sampler, k)
    res = similar(ks)
    for (i,k) in enumerate(ks)
        res[i] = solve(solver,new_basis,pts,k)
    end
    return res
end
