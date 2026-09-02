"""
    solve_wavenumber(solver::SweepSolver,basis::AbsBasis,billiard::BilliardGeometry.AbsBilliard,k,dk;multithreaded::Bool=true) → (k0::Real,t0::Real)

Finds the wavenumber `k0` within `[k-dk/2,k+dk/2]` that minimizes the tension
computed by the sweep `solver`, together with the minimal tension `t0`.

## Description
The basis dimension is scaled with the boundary length and wavenumber via
`solver.dim_scaling_factor` (bounded below by `solver.min_dim`), the basis is
resized to this dimension with [`resize_basis`](@ref), boundary points are
generated once with [`evaluate_points`](@ref), and the tension
`solve(solver,new_basis,pts,k;multithreaded)` is minimized over `k` in the
search window with `Optim.optimize`.

## Arguments
* `solver`: The [`SweepSolver`](@ref) used to solve for the tension at each wavenumber.
* `basis`: The basis used to approximate the eigenstate.
* `billiard`: The billiard whose boundary is discretized.
* `k`: The center of the wavenumber search window.
* `dk`: Width of the wavenumber search window, `[k-dk/2,k+dk/2]`.

## Keyword Arguments
* `multithreaded::Bool=true`: Whether the matrix construction is multithreaded.

## Returns
* `k0`: The wavenumber minimizing the tension within the search window.
* `t0`: The minimal tension found at `k0`.
"""
function solve_wavenumber(solver::SweepSolver,basis::AbsBasis,billiard::BilliardGeometry.AbsBilliard,k,dk;multithreaded::Bool=true)
    L=CompositeCurve(BilliardGeometry.get_boundary_curves(billiard)).length
    dim=max(solver.min_dim,round(Int,L*k*solver.dim_scaling_factor/(2*pi)))
    new_basis=resize_basis(basis,billiard,dim,k)
    pts=evaluate_points(solver,billiard,k)
    function f(k)
        return solve(solver,new_basis,pts,k;multithreaded=multithreaded)
    end
    res=optimize(f,k-0.5*dk,k+0.5*dk)
    k0,t0=res.minimizer,res.minimum
    return k0,t0
end

"""
    k_sweep(solver::SweepSolver,basis::AbsBasis,billiard::BilliardGeometry.AbsBilliard,ks;multithreaded::Bool=true) → res::Vector

Computes the tension of the sweep `solver` at every wavenumber in `ks`, using a
single basis resized to the largest wavenumber in `ks`.

## Description
The basis dimension is scaled with the boundary length and `maximum(ks)` via
`solver.dim_scaling_factor` (bounded below by `solver.min_dim`), the basis is
resized to this dimension with [`resize_basis`](@ref), boundary points are
generated once with [`evaluate_points`](@ref), and [`solve`](@ref) is called for
every wavenumber in `ks`.

## Arguments
* `solver`: The [`SweepSolver`](@ref) used to solve for the tension at each wavenumber.
* `basis`: The basis used to approximate the eigenstate.
* `billiard`: The billiard whose boundary is discretized.
* `ks`: Vector or range of wavenumbers at which the tension is evaluated.

## Keyword Arguments
* `multithreaded::Bool=true`: Whether the matrix construction is multithreaded.

## Returns
* `res`: Vector of tensions, one for each wavenumber in `ks`.
"""
function k_sweep(solver::SweepSolver,basis::AbsBasis,billiard::BilliardGeometry.AbsBilliard,ks;multithreaded::Bool=true)
    k=maximum(ks)
    L=CompositeCurve(BilliardGeometry.get_boundary_curves(billiard)).length
    dim=max(solver.min_dim,round(Int,L*k*solver.dim_scaling_factor/(2*pi)))
    new_basis=resize_basis(basis,billiard,dim,k)
    pts=evaluate_points(solver,billiard,k)
    res=similar(ks)
    @showprogress for (i,k) in enumerate(ks)
        res[i]=solve(solver,new_basis,pts,k;multithreaded=multithreaded)
    end
    return res
end