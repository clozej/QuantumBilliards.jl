#include("../abstracttypes.jl")
#include("../utils/billiardutils.jl")
#include("../utils/typeutils.jl")

struct Eigenstate{K,T,Ba,Bi} <: StationaryState
    k::K
    k_basis::K
    vec::Vector{T}
    ten::T
    dim::Int64
    eps::T
    basis::Ba
    billiard::Bi
end

function Eigenstate(k, vec, ten, basis, billiard)  
    eps = set_precision(vec[1])
    if eltype(vec) <: Real
        filtered_vec = eltype(vec).([abs(v)>eps ? v : zero(vec[1]) for v in vec])
    else 
        filtered_vec = vec
    end
    return Eigenstate(k, k, filtered_vec,ten, length(vec), eps, basis, billiard)
end

function Eigenstate(k, k_basis, vec, ten, basis, billiard)  
    eps = set_precision(vec[1])
    if eltype(vec) <: Real
        filtered_vec = eltype(vec).([abs(v)>eps ? v : zero(vec[1]) for v in vec])
    else 
        filtered_vec = vec
    end
    return Eigenstate(k, k_basis, filtered_vec, ten, length(vec), eps, basis, billiard)
end

function compute_eigenstate(solver::SweepSolver, basis::AbsBasis, billiard::AbsBilliard,k)
    L = billiard.length
    dim = max(solver.min_dim,round(Int, L*k*solver.dim_scaling_factor/(2*pi)))
    basis_new = resize_basis(basis,billiard, dim, k)
    pts = evaluate_points(solver, billiard, k)
    ten, vec = solve_vect(solver, basis_new, pts, k)
    return Eigenstate(k, vec, ten, basis_new, billiard)
end

function compute_eigenstate(solver::AcceleratedSolver, basis::AbsBasis, billiard::AbsBilliard, k; dk = 0.1)
    L = billiard.length
    dim = max(solver.min_dim,round(Int, L*k*solver.dim_scaling_factor/(2*pi)))
    basis_new = resize_basis(basis,billiard,dim,k)
    pts = evaluate_points(solver, billiard, k)
    ks, tens, X = solve_vectors(solver,basis_new, pts, k, dk)
    idx = findmin(abs.(ks.-k))[2]
    k_state = ks[idx]
    ten = tens[idx]
    vec = X[:,idx]
    return Eigenstate(k_state, k, vec, ten, basis_new, billiard)
end

struct EigenstateBundle{K,T,Ba,Bi} <: AbsState 
    ks::Vector{K}
    k_basis::K
    X::Matrix{T}
    tens::Vector{T}
    dim::Int64
    eps::T
    basis::Ba
    billiard::Bi
end

function EigenstateBundle(ks, k_basis, X, tens, basis, billiard)  
    eps = set_precision(X[1,1])
    type = eltype(X)
    if  type <: Real
        filtered_array = type.([abs(x)>eps ? x : zero(type) for x in X])
    else 
        filtered_array = X
    end
    return EigenstateBundle(ks, k_basis, filtered_array, tens, length(X[:,1]), eps, basis, billiard)
end

function compute_eigenstate_bundle(solver::AcceleratedSolver, basis::AbsBasis, billiard::AbsBilliard, k; dk = 0.1, tol=1e-5)
    L = billiard.length
    dim = max(solver.min_dim,round(Int, L*k*solver.dim_scaling_factor/(2*pi)))
    basis_new = resize_basis(basis,billiard, dim,k)
    pts = evaluate_points(solver, billiard, k)
    ks, tens, X = solve_vectors(solver,basis_new, pts, k, dk)
    idx = abs.(tens) .< tol
    ks = ks[idx]
    tens = tens[idx]
    X = X[:,idx]
    return EigenstateBundle(ks, k, X, tens, basis_new, billiard)
end