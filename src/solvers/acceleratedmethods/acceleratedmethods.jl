#include("../../abstracttypes.jl")
#include("../../utils/billiardutils.jl")
#include("decompositions.jl")
#include("../samplers.jl")
include("verginisaraceno.jl")

"""
    solve_wavenumber(solver::AcceleratedBasisSolver, basis::AbsBasis, billiard::AbsBilliard, k, dk; multithreaded::Bool = true) → (k0, t0)

Finds the eigenvalue candidate `k0` closest to the target wavenumber `k`, together
with its associated tension `t0`, among the eigenvalues found by an accelerated
`solver` in a single diagonalization near `k`.

## Description
The basis dimension is scaled with the boundary length and wavenumber via
`solver.dim_scaling_factor` (bounded below by `solver.min_dim`), the basis is
resized to this dimension with [`resize_basis`](@ref), boundary points are
generated with [`evaluate_points`](@ref), and the accelerated eigenvalue problem is
solved with [`solve`](@ref) to obtain candidate wavenumbers `ks` and tensions `ts`
within the window `dk` of `k`. The candidate closest to `k` is returned.

## Arguments
* `solver`: The [`AcceleratedBasisSolver`](@ref) used to solve the eigenvalue problem.
* `basis`: The basis used to approximate the eigenstates.
* `billiard`: The billiard whose boundary is discretized.
* `k`: The target wavenumber around which the search is performed.
* `dk`: Half-width of the wavenumber window around `k`.

## Keyword arguments
* `multithreaded::Bool = true`: Whether the matrix construction is multithreaded.

## Returns
* `k0`: The candidate wavenumber closest to `k`.
* `t0`: The tension associated with `k0`.
"""
function solve_wavenumber(solver::AcceleratedBasisSolver,basis::AbsBasis, billiard::AbsBilliard, k, dk; multithreaded = true)
    L = CompositeCurve(get_boundary_curves(billiard)).length
    dim = max(solver.min_dim,round(Int, L*k*solver.dim_scaling_factor/(2*pi)))
    new_basis = resize_basis(basis,billiard,dim,k)
    pts = evaluate_points(solver, billiard, k)
    ks, ts = solve(solver,new_basis,pts,k,dk; multithreaded)
    idx = findmin(abs.(ks.-k))[2]
    return ks[idx], ts[idx]
end


"""
    solve_spectrum(solver::AcceleratedBasisSolver, basis::AbsBasis, billiard::AbsBilliard, k, dk; multithreaded::Bool = true) → (ks::Vector, ts::Vector)

Computes all eigenvalue candidates and their tensions found by an accelerated
`solver` in a single diagonalization within the wavenumber window `dk` of `k`.

## Description
The basis dimension is scaled with the boundary length and wavenumber via
`solver.dim_scaling_factor` (bounded below by `solver.min_dim`), the basis is
resized to this dimension with [`resize_basis`](@ref), boundary points are
generated with [`evaluate_points`](@ref), and the accelerated eigenvalue problem is
solved with [`solve`](@ref).

## Arguments
* `solver`: The [`AcceleratedBasisSolver`](@ref) used to solve the eigenvalue problem.
* `basis`: The basis used to approximate the eigenstates.
* `billiard`: The billiard whose boundary is discretized.
* `k`: The wavenumber around which the sweep is performed.
* `dk`: Half-width of the wavenumber window around `k`.

## Keyword arguments
* `multithreaded::Bool = true`: Whether the matrix construction is multithreaded.

## Returns
* `ks`: Vector of candidate wavenumbers found within the window.
* `ts`: Vector of tensions associated with `ks`.
"""
function solve_spectrum(solver::AcceleratedBasisSolver,basis::AbsBasis, billiard::AbsBilliard, k, dk; multithreaded = true)
    L = CompositeCurve(get_boundary_curves(billiard)).length
    dim = max(solver.min_dim,round(Int, L*k*solver.dim_scaling_factor/(2*pi)))
    new_basis = resize_basis(basis,billiard,dim,k)
    pts = evaluate_points(solver, billiard, k)
    ks, ts = solve(solver,new_basis,pts,k,dk; multithreaded)
    return ks, ts
end