#include("../abstracttypes.jl")
#include("../utils/billiardutils.jl")
#include("../utils/typeutils.jl")

"""
The struct defines the parameters an eigenstate must have and also certain types designed for streamline computations.
Comment: This struct is low-level and is not instantied by it, rather a function (high level wrapper) exists for it (`compute_eigenstate`)

# Fields
- `k`: The wave number.
- `vec`: The eigenstate vector.
- `ten`: The tension at that wavenumber `k`.
- `basis`: The basis used to generate the eigenstate.
- `billiard`: The billiard geometry in question.
- `eps`: The precision of the given float type
- `dim`: The dimension of the basis
"""
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

"""
Constructor for the `Eigenstate` struct.

# Arguments:
- `k::T`: The wave number.
- `vec::Vector{T}`: The eigenstate vector.
- `ten::T`: The tension at that wavenumber `k`.
- `basis<:AbsBasis`: The basis used to generate the eigenstate.
- `billiard<:AbsBilliard{T}`: The billiard geometry in question.

# Returns:
- An instance of `Eigenstate`.
"""
function Eigenstate(k, vec, ten, basis, billiard)  
    eps = set_precision(vec[1])
    if eltype(vec) <: Real
        filtered_vec = eltype(vec).([abs(v)>eps ? v : zero(vec[1]) for v in vec])
    else 
        filtered_vec = vec
    end
    return Eigenstate(k, k, filtered_vec,ten, length(vec), eps, basis, billiard)
end

#TODO WTF is k_basis, no idea from struct definition
function Eigenstate(k, k_basis, vec, ten, basis, billiard)  
    eps = set_precision(vec[1])
    if eltype(vec) <: Real
        filtered_vec = eltype(vec).([abs(v)>eps ? v : zero(vec[1]) for v in vec])
    else 
        filtered_vec = vec
    end
    return Eigenstate(k, k_basis, filtered_vec, ten, length(vec), eps, basis, billiard)
end

"""
This function computes a quantum eigenstate for a given wavenumber `k` using a `SweepSolver`.

# Logic
- First it calucaltes the dimension of the basis and the lenth of the boundary from the solver and the boundary params.
- It then resizes the basis to match the new dimension and the given wavenumber `k (if below a minimum defined by the basis rules then use that)
- Evaluates the points for the billiard at that specific k wavenumber.
- Solves a SVD problem and gathers the minimum singular value (the "tension") and also the associated right eigenvector.
- From these constructs the Eigenstate struct

# Arguments
- `solver::SweepSolver`: The solver used to compute the eigenstate.
- `basis::AbsBasis`: The basis functions used.
- `billiard::AbsBilliard`: The billiard geometry on which the problem is defined.
- `k`: Initial wavenumber guess.

# Returns
- `Eigenstate{K,T,Ba,Bi}`: The computed quantum eigenstate.
"""
function compute_eigenstate(solver::SweepSolver, basis::AbsBasis, billiard::AbsBilliard,k)
    L = billiard.length
    dim = max(solver.min_dim,round(Int, L*k*solver.dim_scaling_factor/(2*pi)))
    basis_new = resize_basis(basis,billiard, dim, k)
    pts = evaluate_points(solver, billiard, k)
    ten, vec = solve_vect(solver, basis_new, pts, k)
    return Eigenstate(k, vec, ten, basis_new, billiard)
end

"""
Compute a quantum eigenstate for a given wavenumber `k` using the logic of Vergini-Saraceno method with an `AcceleratedSolver`.

# Logic
- The method calculates the dimension of the basis and the length of the boundary using the solver and billiard parameters.
- The basis is resized to match the new dimension and the given wavenumber `k`. If the dimension is below a minimum defined by the basis rules, it uses that minimum.
- The function evaluates the points on the billiard at the specific wavenumber `k`.
- It then solves a generalized eigenvalue problem across a small range of wavenumbers around `k`.
- Among the obtained eigenstates, the one with the smallest tension is selected.
- The selected eigenstate is then used to construct and return the `Eigenstate` struct.

# Arguments
- `solver::AcceleratedSolver`: The solver used to compute the eigenstate.
- `basis::AbsBasis`: The basis functions used.
- `billiard::AbsBilliard`: The billiard geometry on which the problem is defined.
- `k`: Initial wavenumber guess.
- `dk`: Range around `k` within which to search for the eigenstate (default is 0.1).

# Returns
- `Eigenstate{K,T,Ba,Bi}`: The computed quantum eigenstate with the smallest tension.
"""
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

"""
A struct representing a bundle of quantum eigenstates.

# Fields
- `ks::Vector{K}`: A vector of wavenumbers associated with the eigenstates.
- `X::Matrix{T}`: A matrix where each column represents an eigenvector associated with a particular wavenumber in `ks`.
- `tens::Vector{T}`: A vector of tension values corresponding to each eigenstate.
- `dim::Int64`: The dimension of the basis used to generate the eigenstates.
- `eps::T`: The precision used for filtering small values in the eigenstate vectors.
- `basis::Ba`: The basis used to represent the eigenstates.
- `billiard::Bi`: The billiard geometry associated with the eigenstates.
"""
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

"""
Construct an `EigenstateBundle` containing a set of quantum eigenstates.

# Logic
- The constructor initializes an `EigenstateBundle` struct, given a set of wavenumbers, eigenvectors, tension values, basis, and billiard geometry.
- It sets a precision value `eps` based on the type of the eigenvector matrix `X` to filter out small values in the eigenvectors.
- If the type of `X` is real, the constructor filters out elements in `X` that are smaller than `eps`, replacing them with zeros.
- The filtered eigenvector matrix, along with the other provided parameters, is then used to construct and return an `EigenstateBundle` instance.

# Arguments
- `ks::Vector{K}`: A vector of wavenumbers associated with the eigenstates.
- `X::Matrix{T}`: A matrix where each column represents an eigenvector associated with a particular wavenumber in `ks`.
- `tens::Vector{T}`: A vector of tension values corresponding to each eigenstate.
- `basis::Ba`: The basis used to represent the eigenstates.
- `billiard::Bi`: The billiard geometry associated with the eigenstates.

# Returns
- `EigenstateBundle{K,T,Ba,Bi}`: An `EigenstateBundle` instance representing a bundle of quantum eigenstates.
"""
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

"""
Compute a bundle of quantum eigenstates around a given wavenumber `k` using the logic of the Vergini-Saraceno method with an `AcceleratedSolver`.

# Logic
- The method calculates the dimension of the basis and the length of the boundary using the solver and billiard parameters.
- The basis is resized to match the new dimension and the given wavenumber `k`. If the dimension is below a minimum defined by the basis rules, it uses that minimum.
- The function evaluates the points on the billiard at the specific wavenumber `k`.
- It then solves a generalized eigenvalue problem across a small range of wavenumbers around `k`.
- Among the obtained eigenstates, those with tensions below a specified tolerance `tol` are selected.
- The selected eigenstates are then used to construct and return an `EigenstateBundle` struct.

# Arguments
- `solver::AcceleratedSolver`: The solver used to compute the eigenstates.
- `basis::AbsBasis`: The basis functions used.
- `billiard::AbsBilliard`: The billiard geometry on which the problem is defined.
- `k`: Initial wavenumber guess.
- `dk`: Range around `k` within which to search for eigenstates (default is 0.1).
- `tol`: Tolerance for selecting eigenstates based on their tension (default is 1e-5).

# Returns
- `EigenstateBundle{K,T,Ba,Bi}`: A bundle of computed quantum eigenstates with tensions below the specified tolerance.
"""
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