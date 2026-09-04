#include("../abstracttypes.jl")
#include("../utils/billiardutils.jl")
#include("../utils/typeutils.jl")

"""
    Eigenstate{K,T,S,Bi,Ba} <: StationaryState

`Eigenstate` is a concrete type representing a numerically computed eigenstate
of a quantum billiard at a given wavenumber.

## Description
Eigenstates are produced by [`compute_eigenstate`](@ref), which combines a
sweep or accelerated solver, a basis and a billiard geometry to obtain the
expansion coefficients `vec` in `basis` and an estimate of the tension `ten`,
quantifying how well the boundary condition is satisfied. Coefficients
smaller in magnitude than the numerical precision `eps` (given by
`set_precision`) are set to zero when the state is constructed.

## Attributes
* `k`: The wavenumber of the eigenstate, as refined by the solver.
* `k_basis`: The wavenumber at which `basis` was evaluated to obtain `vec` (may differ slightly from `k` for accelerated solvers).
* `vec`: Expansion coefficients of the eigenstate in `basis`.
* `ten`: Tension of the solution, measuring the residual boundary condition violation.
* `dim`: Dimension of `vec` (and of `basis`).
* `eps`: Numerical precision threshold below which coefficients of `vec` are treated as zero.
* `solver`: The solver (`S<:AbsBasisSolver`) used to compute the eigenstate.
* `basis`: The basis (`Ba<:AbsBasis`), resized/evaluated at `k_basis`, in which `vec` is expressed.
* `billiard`: The billiard (`Bi<:AbsBilliard`) the eigenstate is defined on.

## API
The following functions can be evaluated for this type:
- [`compute_eigenstate`](@ref)
- [`boundary_function`](@ref)
- [`momentum_function`](@ref)
- [`wavefunction`](@ref)
- [`husimi_function`](@ref)
"""
struct Eigenstate{K,T,S,Bi,Ba} <: StationaryState
    k::K
    k_basis::K
    vec::Vector{K}
    ten::T
    dim::Int64
    eps::T
    solver::S
    basis::Ba
    billiard::Bi
end

"""
    Eigenstate(k, vec, ten, solver, basis, billiard) → state::Eigenstate

Construct an [`Eigenstate`](@ref) with `k_basis` set equal to `k`, filtering
out negligible coefficients of `vec`.

## Description
The numerical precision `eps` is obtained from `set_precision` applied to
`vec[1]`. If `vec` has real entries, any entry with absolute value not
exceeding `eps` is replaced by zero; a complex-valued `vec` is left
unfiltered.

## Arguments
* `k`: The wavenumber of the eigenstate.
* `vec`: Expansion coefficients of the eigenstate in `basis`.
* `ten`: Tension of the solution.
* `solver`: The solver used to compute the eigenstate.
* `basis`: The basis in which `vec` is expressed.
* `billiard`: The billiard the eigenstate is defined on.

## Returns
*  `state` : A new [`Eigenstate`](@ref) with `k_basis = k` and filtered coefficients.
"""
function Eigenstate(k, vec, ten, solver, basis, billiard)  
    eps = set_precision(vec[1])
    if eltype(vec) <: Real
        filtered_vec = eltype(vec).([abs(v)>eps ? v : zero(vec[1]) for v in vec])
    else 
        filtered_vec = vec
    end
    return Eigenstate(k, k, filtered_vec,ten, length(vec), eps, solver, basis, billiard)
end

"""
    Eigenstate(k, k_basis, vec, ten, solver, basis, billiard) → state::Eigenstate

Construct an [`Eigenstate`](@ref) allowing the refined wavenumber `k` and the
basis-evaluation wavenumber `k_basis` to differ, filtering out negligible
coefficients of `vec`.

## Description
Behaves as [`Eigenstate(k, vec, ten, solver, basis, billiard)`](@ref Eigenstate),
except that `k_basis` is taken as given instead of being set equal to `k`.
This is used by accelerated solvers, where `basis` is evaluated at a fixed
scaling wavenumber `k_basis` while the eigenstate itself is refined to a
nearby wavenumber `k`.

## Arguments
* `k`: The (refined) wavenumber of the eigenstate.
* `k_basis`: The wavenumber at which `basis` was evaluated to obtain `vec`.
* `vec`: Expansion coefficients of the eigenstate in `basis`.
* `ten`: Tension of the solution.
* `solver`: The solver used to compute the eigenstate.
* `basis`: The basis in which `vec` is expressed.
* `billiard`: The billiard the eigenstate is defined on.

## Returns
*  `state` : A new [`Eigenstate`](@ref) with filtered coefficients.
"""
function Eigenstate(k, k_basis, vec, ten, solver, basis, billiard)  
    eps = set_precision(vec[1])
    if eltype(vec) <: Real
        filtered_vec = eltype(vec).([abs(v)>eps ? v : zero(vec[1]) for v in vec])
    else 
        filtered_vec = vec
    end
    return Eigenstate(k, k_basis, filtered_vec, ten, length(vec), eps, solver, basis, billiard)
end

"""
    compute_eigenstate(solver::SweepBasisSolver, basis::AbsBasis, billiard::AbsBilliard, k; multithreaded::Bool = true) → state::Eigenstate

Computes the [`Eigenstate`](@ref) of `billiard` at wavenumber `k` using a
sweep-method `solver` (e.g. `DecompositionMethodSolver`).

## Description
The basis dimension is set to
`dim = max(solver.min_dim, round(Int, L*k*solver.dim_scaling_factor/(2*pi)))`,
with `L` the total boundary length, and `basis` is resized to this dimension
with `resize_basis`. Boundary points are sampled with `evaluate_points`, and
the generalized eigenvalue problem is solved at `k` with `solve_vect` to
obtain the tension `ten` and coefficient vector `vec`.

## Arguments
* `solver`: The `SweepBasisSolver` used to solve the eigenvalue problem.
* `basis`: The basis used to approximate the eigenstate.
* `billiard`: The billiard the eigenstate is computed on.
* `k`: The wavenumber at which the eigenstate is computed.

## Keyword arguments
*  `multithreaded::Bool = true` : Whether the matrix construction is multithreaded.

## Returns
*  `state` : The computed [`Eigenstate`](@ref) at wavenumber `k`.
"""
function compute_eigenstate(solver::SweepBasisSolver, basis::AbsBasis, billiard::AbsBilliard,k; multithreaded = true)
    L = CompositeCurve(get_boundary_curves(billiard)).length
    dim = max(solver.min_dim,round(Int, L*k*solver.dim_scaling_factor/(2*pi)))
    basis_new = resize_basis(basis,billiard, dim, k)
    pts = evaluate_points(solver, billiard, k)
    ten, vec = solve_vect(solver, basis_new, pts, k; multithreaded)
    return Eigenstate(k, vec, ten, solver, basis_new, billiard)
end

"""
    compute_eigenstate(solver::AcceleratedBasisSolver, basis::AbsBasis, billiard::AbsBilliard, k; dk::Real = 0.1, multithreaded::Bool = true) → state::Eigenstate

Computes the [`Eigenstate`](@ref) of `billiard` closest to wavenumber `k`
using an accelerated `solver` (e.g. `VerginiSaracenoSolver`).

## Description
The basis dimension is set to
`dim = max(solver.min_dim, round(Int, L*k*solver.dim_scaling_factor/(2*pi)))`,
with `L` the total boundary length, and `basis` is resized to this dimension
with `resize_basis`. Boundary points are sampled with `evaluate_points`, and
`solve_vectors` is used to find all candidate wavenumbers `ks`, tensions
`tens` and eigenvectors `X` within `dk` of `k`. The candidate `k_state`
closest to `k` is selected and used to build the resulting [`Eigenstate`](@ref),
whose `k_basis` is set to the requested `k` (the wavenumber at which `basis`
was evaluated).

## Arguments
* `solver`: The `AcceleratedBasisSolver` used to solve the eigenvalue problem.
* `basis`: The basis used to approximate the eigenstate.
* `billiard`: The billiard the eigenstate is computed on.
* `k`: The target wavenumber around which the eigenstate is searched for.

## Keyword arguments
*  `dk::Real = 0.1` : Half-width of the wavenumber window around `k` within which candidate eigenstates are searched.
*  `multithreaded::Bool = true` : Whether the matrix construction is multithreaded.

## Returns
*  `state` : The computed [`Eigenstate`](@ref) closest to wavenumber `k`.
"""
function compute_eigenstate(solver::AcceleratedBasisSolver, basis::AbsBasis, billiard::AbsBilliard, k; dk = 0.1, multithreaded = true)
    L = CompositeCurve(get_boundary_curves(billiard)).length
    dim = max(solver.min_dim,round(Int, L*k*solver.dim_scaling_factor/(2*pi)))
    basis_new = resize_basis(basis,billiard,dim,k)
    pts = evaluate_points(solver, billiard, k)
    ks, tens, X = solve_vectors(solver,basis_new, pts, k, dk; multithreaded)
    idx = findmin(abs.(ks.-k))[2]
    k_state = ks[idx]
    ten = tens[idx]
    vec = X[:,idx]
    return Eigenstate(k_state, k, vec, ten, solver, basis_new, billiard)
end
