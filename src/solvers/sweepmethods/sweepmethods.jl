#include("../../abstracttypes.jl")
#include("../../utils/billiardutils.jl")
#include("decompositions.jl")
#include("../samplers.jl")
include("particularsolutionsmethod.jl")
include("decompositionmethod.jl")
using LinearAlgebra, Optim

"""
Solve for the wavenumber using the sweep solver around an initial value of k in the interval [k - 0.5*dk, k + 0.5*dk]  .

# Description
This function finds the wavenumber around an initial `k` that minimizes the the SVD value given by the sweep solver.

# Arguments
- `solver::SweepSolver`: The sweep solver used to solve the problem.
- `basis::AbsBasis`: The basis functions used.
- `billiard::AbsBilliard`: The billiard geometry on which the problem is defined.
- `k`: Initial wavenumber guess.
- `dk`: Range around `k` within which to optimize.

# Returns
- `k0`: The optimized wavenumber.
- `t0`: The corresponding minimum value of the solution at `k0`.
"""
function solve_wavenumber(solver::SweepSolver,basis::AbsBasis, billiard::AbsBilliard, k, dk)
    dim = max(solver.min_dim,round(Int, billiard.length*k*solver.dim_scaling_factor/(2*pi)))
    new_basis = resize_basis(basis,billiard,dim,k)
    pts = evaluate_points(solver, billiard, k)
    function f(k)
        return solve(solver,new_basis,pts,k)
    end
    res =  optimize(f, k-0.5*dk, k+0.5*dk)
    k0,t0 = res.minimizer, res.minimum
    return k0, t0
end

"""
Perform a sweep over a range of wavenumbers to evaluate the corresponding SVD values using the sweep solver.

# Description
This function evaluates the SVD values at each wavenumber in a provided set `ks` using the sweep solver. For each wavenumber, the solver computes the corresponding solution, which represents the minimum singular value for that wavenumber.

# Arguments
- `solver::SweepSolver`: The sweep solver used to compute the SVD values.
- `basis::AbsBasis`: The basis functions used for the solution.
- `billiard::AbsBilliard`: The billiard geometry on which the problem is defined.
- `ks`: A vector of wavenumbers at which to evaluate the SVD values.

# Returns
- `res`: A vector of SVD values corresponding to each wavenumber in `ks`.
"""
function k_sweep(solver::SweepSolver, basis::AbsBasis, billiard::AbsBilliard, ks)
    k = maximum(ks)
    dim = max(solver.min_dim,round(Int, billiard.length*k*solver.dim_scaling_factor/(2*pi)))
    new_basis = resize_basis(basis,billiard,dim,k)
    pts = evaluate_points(solver, billiard, k)
    res = similar(ks)
    for (i,k) in enumerate(ks)
        res[i] = solve(solver,new_basis,pts,k)
    end
    return res
end
