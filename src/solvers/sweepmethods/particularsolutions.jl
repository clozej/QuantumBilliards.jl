"""
    ParticularSolutionsMethod{T} <: SweepBasisSolver

`ParticularSolutionsMethod` is a concrete [`SweepBasisSolver`](@ref)
implementing the particular solutions method (PSM) for computing quantum
billiard spectra by sweeping over individual wavenumbers.

## Description
For a fixed wavenumber `k`, the method constructs the basis matrices `B`
(boundary) and `B_int` (interior) (see [`construct_matrices`](@ref)) from a
boundary quadrature and a set of random interior points (see
[`evaluate_points`](@ref)), and defines the tension from the smallest singular
value of `B` normalized against `B_int` (see [`solve`](@ref)). A sequence of
tensions over a range of wavenumbers is minimized/scanned by
[`solve_wavenumber`](@ref) or [`k_sweep`](@ref), inherited from
[`SweepBasisSolver`](@ref), to locate the eigenvalues of the billiard.

## Attributes
* `dim_scaling_factor`: Scaling factor used to determine the basis dimension from the boundary length and wavenumber.
* `pts_scaling_factor`: Vector of scaling factors, one per fundamental boundary curve, used to determine the number of boundary sampling points.
* `int_pts_scaling_factor`: Scaling factor used to determine the number of interior sampling points.
* `sampler`: Vector of samplers, one per fundamental boundary curve, used to generate boundary points.
* `eps`: Relative tolerance used to filter small singular values.
* `min_dim`: Minimum basis dimension.
* `min_pts`: Minimum number of boundary sampling points.
* `min_int_pts`: Minimum number of interior sampling points.

## API
The following functions can be evaluated for this type:
- [`evaluate_points`](@ref)
- [`construct_matrices`](@ref)
- [`solve_full`](@ref)
- [`solve_with_rank_reduction`](@ref)
- [`solve`](@ref)
- [`solve_vect`](@ref)
- [`solve_wavenumber`](@ref)
- [`k_sweep`](@ref)

!!! note "Migration status"
    API scaffold only (struct, constructor and method signatures) ported from
    `QuantumBilliards-develop`'s `ParticularSolutionsMethod`. The
    matrix-assembly bodies are not yet implemented; every method below raises
    an `error` until Step 2 of the migration plan lands. `solve_wavenumber` and
    `k_sweep` are already usable once `evaluate_points`/`solve` are
    implemented, since they are inherited for free from the shared
    [`SweepBasisSolver`](@ref) generics.
"""
struct ParticularSolutionsMethod{T} <: SweepBasisSolver where {T<:Real}
    dim_scaling_factor::T
    pts_scaling_factor::Vector{T}
    int_pts_scaling_factor::T
    sampler::Vector
    eps::T
    min_dim::Int64
    min_pts::Int64
    min_int_pts::Int64
end

"""
    ParticularSolutionsMethod(dim_scaling_factor::T, pts_scaling_factor::Union{T,Vector{T}}, int_pts_scaling_factor::T; min_dim::Int = 100, min_pts::Int = 500, min_int_pts::Int = 500) where {T<:Real} → solver::ParticularSolutionsMethod{T}

Constructs a [`ParticularSolutionsMethod`](@ref) with a single
`GaussLegendreNodes` sampler shared by every fundamental boundary curve.

## Arguments
* `dim_scaling_factor`: Scaling factor used to determine the basis dimension.
* `pts_scaling_factor`: Scaling factor, or vector thereof (one per fundamental boundary curve), used to determine the number of boundary sampling points.
* `int_pts_scaling_factor`: Scaling factor used to determine the number of interior sampling points.

## Keyword arguments
* `min_dim::Int = 100`: Minimum basis dimension.
* `min_pts::Int = 500`: Minimum number of boundary sampling points.
* `min_int_pts::Int = 500`: Minimum number of interior sampling points.

## Returns
* `solver`: A [`ParticularSolutionsMethod{T}`](@ref) instance.
"""
function ParticularSolutionsMethod(dim_scaling_factor::T, pts_scaling_factor::Union{T,Vector{T}},
                                    int_pts_scaling_factor::T;
                                    min_dim::Int=100, min_pts::Int=500, min_int_pts::Int=500) where {T<:Real}
    d = dim_scaling_factor
    bs = pts_scaling_factor isa T ? [pts_scaling_factor] : pts_scaling_factor
    sampler = [GaussLegendreNodes()]
    return ParticularSolutionsMethod(d, bs, int_pts_scaling_factor, sampler, eps(T), min_dim, min_pts, min_int_pts)
end

"""
    ParticularSolutionsMethod(dim_scaling_factor::T, pts_scaling_factor::Union{T,Vector{T}}, int_pts_scaling_factor::T, samplers::Vector{<:AbsSampler}; min_dim::Int = 100, min_pts::Int = 500, min_int_pts::Int = 500) where {T<:Real} → solver::ParticularSolutionsMethod{T}

Constructs a [`ParticularSolutionsMethod`](@ref) with a user-supplied sampler
for each fundamental boundary curve.

## Arguments
* `dim_scaling_factor`: Scaling factor used to determine the basis dimension.
* `pts_scaling_factor`: Scaling factor, or vector thereof (one per fundamental boundary curve), used to determine the number of boundary sampling points.
* `int_pts_scaling_factor`: Scaling factor used to determine the number of interior sampling points.
* `samplers`: Vector of samplers, one per fundamental boundary curve.

## Keyword arguments
* `min_dim::Int = 100`: Minimum basis dimension.
* `min_pts::Int = 500`: Minimum number of boundary sampling points.
* `min_int_pts::Int = 500`: Minimum number of interior sampling points.

## Returns
* `solver`: A [`ParticularSolutionsMethod{T}`](@ref) instance.
"""
function ParticularSolutionsMethod(dim_scaling_factor::T, pts_scaling_factor::Union{T,Vector{T}},
                                    int_pts_scaling_factor::T, samplers::Vector{<:AbsSampler};
                                    min_dim::Int=100, min_pts::Int=500, min_int_pts::Int=500) where {T<:Real}
    d = dim_scaling_factor
    bs = pts_scaling_factor isa T ? [pts_scaling_factor] : pts_scaling_factor
    return ParticularSolutionsMethod(d, bs, int_pts_scaling_factor, samplers, eps(T), min_dim, min_pts, min_int_pts)
end

const _PSM_NOT_IMPLEMENTED = "ParticularSolutionsMethod matrix assembly is not yet implemented (API scaffold only, see the QuantumBilliardsTests migration plan)."

"""
    evaluate_points(solver::ParticularSolutionsMethod, billiard::Bi, k) where {Bi<:AbsBilliard} → pts::BoundaryPoints

Samples the boundary of `billiard`, together with a set of random interior
points, needed to construct the matrices in [`construct_matrices`](@ref).

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function evaluate_points(solver::ParticularSolutionsMethod, billiard::Bi, k) where {Bi<:AbsBilliard}
    error(_PSM_NOT_IMPLEMENTED)
end

"""
    construct_matrices(solver::ParticularSolutionsMethod, basis::Ba, pts::BoundaryPoints, k; multithreaded::Bool = true) where {Ba<:AbsBasis} → (B::Matrix, B_int::Matrix)

Constructs the boundary basis matrix `B` and the interior basis matrix
`B_int` at wavenumber `k`.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function construct_matrices(solver::ParticularSolutionsMethod, basis::Ba, pts::BoundaryPoints, k; multithreaded::Bool=true) where {Ba<:AbsBasis}
    error(_PSM_NOT_IMPLEMENTED)
end

"""
    solve_full(solver::ParticularSolutionsMethod, basis::Ba, pts::BoundaryPoints, k; multithreaded::Bool = true) where {Ba<:AbsBasis} → t::Real

Computes the particular solutions method tension at wavenumber `k` from the
full (non rank-reduced) singular value decomposition of `B`/`B_int`.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function solve_full(solver::ParticularSolutionsMethod, basis::Ba, pts::BoundaryPoints, k; multithreaded::Bool=true) where {Ba<:AbsBasis}
    error(_PSM_NOT_IMPLEMENTED)
end

"""
    solve_with_rank_reduction(solver::ParticularSolutionsMethod, basis::Ba, pts::BoundaryPoints, k; multithreaded::Bool = true) where {Ba<:AbsBasis} → t::Real

Computes the particular solutions method tension at wavenumber `k` using a
rank-reduced singular value decomposition of `B`/`B_int` for improved
performance at large basis dimension.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function solve_with_rank_reduction(solver::ParticularSolutionsMethod, basis::Ba, pts::BoundaryPoints, k; multithreaded::Bool=true) where {Ba<:AbsBasis}
    error(_PSM_NOT_IMPLEMENTED)
end

"""
    solve(solver::ParticularSolutionsMethod, basis::Ba, pts::BoundaryPoints, k; multithreaded::Bool = true) where {Ba<:AbsBasis} → t::Real

Computes the particular solutions method tension at wavenumber `k` for `basis`
on the boundary/interior points `pts`.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function solve(solver::ParticularSolutionsMethod, basis::Ba, pts::BoundaryPoints, k; multithreaded::Bool=true) where {Ba<:AbsBasis}
    error(_PSM_NOT_IMPLEMENTED)
end

"""
    solve_vect(solver::ParticularSolutionsMethod, basis::Ba, pts::BoundaryPoints, k; multithreaded::Bool = true) where {Ba<:AbsBasis} → (t::Real, x::Vector)

Computes the particular solutions method tension and the corresponding
eigenvector (expressed in the original basis) at wavenumber `k`.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function solve_vect(solver::ParticularSolutionsMethod, basis::Ba, pts::BoundaryPoints, k; multithreaded::Bool=true) where {Ba<:AbsBasis}
    error(_PSM_NOT_IMPLEMENTED)
end
