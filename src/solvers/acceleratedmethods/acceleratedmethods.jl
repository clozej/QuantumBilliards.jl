#include("../../abstracttypes.jl")
#include("../../utils/billiardutils.jl")
#include("decompositions.jl")
#include("../samplers.jl")
include("scalingmethod.jl")

function solve_wavenumber(solver::AcceleratedSolver,basis::AbsBasis, billiard::AbsBilliard, k, dk; sampler=gauss_legendre_nodes)
    dim = round(Int, real_length(billiard)*k*solver.dim_scaling_factor/(2*pi))
    new_basis = resize_basis(basis,dim)
    pts = evaluate_points(solver, billiard, sampler, k)
    ks, ts = solve(solver,new_basis,pts,k,dk)
    idx = findmin(abs.(ks.-k))[2]
    return ks[idx], ts[idx]
end


function solve_spectrum(solver::AcceleratedSolver,basis::AbsBasis, billiard::AbsBilliard, k, dk; sampler=gauss_legendre_nodes)
    dim = round(Int, real_length(billiard)*k*solver.dim_scaling_factor/(2*pi))
    new_basis = resize_basis(basis,dim)
    pts = evaluate_points(solver, billiard, sampler, k)
    ks, ts = solve(solver,new_basis,pts,k,dk)
    return ks, ts
end