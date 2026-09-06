"""
    BeynSolver{T,K} <: AcceleratedBIMSolver

`BeynSolver` is a concrete [`AcceleratedBIMSolver`](@ref) implementing Beyn's
contour-integral method for recovering every root of the nonlinear
eigenproblem `A(k)v = 0` within a wavenumber window from a single contour
solve.

## Description
`BeynSolver` wraps an inner `kernel::SweepBIMSolver` (a
[`DoubleLayerPotentialSolver`](@ref), [`CombinedFieldIntegralEquationSolver`](@ref)
or [`CompositeBIMSolver`](@ref)) supplying the Fredholm operator `A(k)`. On a
circular contour of center `k0` and radius `R`, the contour moments

    A0 = (1/2πi) ∮ A(z)⁻¹V dz,
    A1 = (1/2πi) ∮ z A(z)⁻¹V dz,

are formed against a random probing matrix `V` of rank `r` (see
[`construct_matrices`](@ref)); a rank-revealing SVD of `A0` projects the
problem onto a small dense generalized eigenproblem whose eigenpairs
approximate the roots of `A(k)v = 0` inside the contour (see [`solve`](@ref),
[`solve_vectors`](@ref)). Candidates are filtered by contour containment and
residual norm using `svd_tol`/`res_tol`.

Reference: W.-J. Beyn, "An integral method for solving nonlinear eigenvalue
problems", Linear Algebra Appl. 436 (2012).

## Attributes
* `kernel`: The wrapped [`SweepBIMSolver`](@ref) supplying `A(k)`.
* `m`: Target eigenvalue count per contour window, used for Weyl-window planning.
* `nq`: Number of contour quadrature nodes.
* `r`: Random probing rank.
* `svd_tol`: Singular-value cutoff used for rank detection in `A0`.
* `res_tol`: Residual threshold used to reject spurious roots.
* `auto_discard_spurious`: Whether candidates with residual above `res_tol` are automatically rejected.
* `use_chebyshev`: Whether Chebyshev-accelerated kernel evaluation is used.
* `n_panels_h`: Hankel-function Chebyshev panel count.
* `M_h`: Hankel-function Chebyshev polynomial degree.
* `n_panels_j`: Bessel-J-function Chebyshev panel count.
* `M_j`: Bessel-J-function Chebyshev polynomial degree.

## API
The following functions can be evaluated for this type:
- [`evaluate_points`](@ref)
- [`construct_matrices`](@ref)
- [`solve`](@ref)
- [`solve_vectors`](@ref)
- [`solve_wavenumber`](@ref)
- [`solve_spectrum`](@ref)

!!! note "Migration status"
    API scaffold only (struct, constructor and method signatures) ported from
    `QuantumBilliards-develop`'s `BeynSolver`/`solve_wavenumber_beyn`. The
    matrix-assembly/contour-solve bodies are not yet implemented; every
    method below raises an `error` until Step 2 of the migration plan lands.
"""
struct BeynSolver{T<:Real,K<:SweepBIMSolver} <: AcceleratedBIMSolver
    kernel::K
    m::Int
    nq::Int
    r::Int
    svd_tol::T
    res_tol::T
    auto_discard_spurious::Bool
    use_chebyshev::Bool
    n_panels_h::Int
    M_h::Int
    n_panels_j::Int
    M_j::Int
end

"""
    BeynSolver(kernel::K; m::Int = 10, nq::Int = 48, r::Int = 48, svd_tol::Real = 1e-12, res_tol::Real = 1e-9, auto_discard_spurious::Bool = true, use_chebyshev::Bool = true, n_panels_h::Int = 15000, M_h::Int = 5, n_panels_j::Int = 10000, M_j::Int = 5) where {K<:SweepBIMSolver} → solver::BeynSolver

Constructs a [`BeynSolver`](@ref) wrapping the boundary-integral `kernel`.

## Arguments
* `kernel`: The [`SweepBIMSolver`](@ref) supplying the Fredholm operator `A(k)`.

## Keyword arguments
* `m::Int = 10`: Target eigenvalue count per contour window.
* `nq::Int = 48`: Number of contour quadrature nodes.
* `r::Int = 48`: Random probing rank.
* `svd_tol::Real = 1e-12`: Singular-value cutoff for rank detection.
* `res_tol::Real = 1e-9`: Residual threshold for spurious-root rejection.
* `auto_discard_spurious::Bool = true`: Whether to automatically reject high-residual candidates.
* `use_chebyshev::Bool = true`: Whether to use Chebyshev-accelerated kernel evaluation.
* `n_panels_h::Int = 15000`: Hankel-function Chebyshev panel count.
* `M_h::Int = 5`: Hankel-function Chebyshev polynomial degree.
* `n_panels_j::Int = 10000`: Bessel-J-function Chebyshev panel count.
* `M_j::Int = 5`: Bessel-J-function Chebyshev polynomial degree.

## Returns
* `solver`: A [`BeynSolver`](@ref) instance.
"""
function BeynSolver(kernel::K; m::Int=10, nq::Int=48, r::Int=48,
                     svd_tol::Real=1e-12, res_tol::Real=1e-9,
                     auto_discard_spurious::Bool=true, use_chebyshev::Bool=true,
                     n_panels_h::Int=15000, M_h::Int=5, n_panels_j::Int=10000, M_j::Int=5) where {K<:SweepBIMSolver}
    T = _bim_numeric_type(kernel)
    return BeynSolver{T,K}(kernel, m, nq, r, T(svd_tol), T(res_tol), auto_discard_spurious, use_chebyshev, n_panels_h, M_h, n_panels_j, M_j)
end

const _BEYN_NOT_IMPLEMENTED = "BeynSolver contour assembly is not yet implemented (API scaffold only, see the QuantumBilliardsTests migration plan)."

"""
    construct_matrices(solver::BeynSolver, pts::BoundaryPoints, k0, R; multithreaded::Bool = true) → (A0::Matrix, A1::Matrix)

Assembles the two Beyn contour moments `A0`, `A1` on a circular contour of
center `k0` and radius `R`.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function construct_matrices(solver::BeynSolver, pts::BoundaryPoints, k0, R; multithreaded::Bool=true)
    error(_BEYN_NOT_IMPLEMENTED)
end

"""
    solve(solver::BeynSolver, pts::BoundaryPoints, k0, dk; multithreaded::Bool = true) → (ks::Vector, ts::Vector)

Solves the Beyn contour eigenproblem centered at `k0` with radius `dk`,
returning every retained candidate wavenumber and its tension.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function solve(solver::BeynSolver, pts::BoundaryPoints, k0, dk; multithreaded::Bool=true)
    error(_BEYN_NOT_IMPLEMENTED)
end

"""
    solve_vectors(solver::BeynSolver, pts::BoundaryPoints, k0, dk; multithreaded::Bool = true) → (ks::Vector, ts::Vector, X::Matrix)

Solves the Beyn contour eigenproblem centered at `k0` with radius `dk`,
returning every retained candidate wavenumber, its tension, and the
corresponding boundary density eigenvector.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function solve_vectors(solver::BeynSolver, pts::BoundaryPoints, k0, dk; multithreaded::Bool=true)
    error(_BEYN_NOT_IMPLEMENTED)
end

"""
    solve_wavenumber(solver::BeynSolver, billiard::Bi, k, dk; multithreaded::Bool = true) where {Bi<:AbsBilliard} → (k0::Real, t0::Real)

Finds the Beyn eigenvalue candidate closest to the target wavenumber `k`
within a contour of radius `dk`, together with its tension.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function solve_wavenumber(solver::BeynSolver, billiard::Bi, k, dk; multithreaded::Bool=true) where {Bi<:AbsBilliard}
    error(_BEYN_NOT_IMPLEMENTED)
end

"""
    solve_spectrum(solver::BeynSolver, billiard::Bi, k, dk; multithreaded::Bool = true) where {Bi<:AbsBilliard} → (ks::Vector, ts::Vector)

Computes every Beyn eigenvalue candidate and its tension within a contour of
center `k` and radius `dk`.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function solve_spectrum(solver::BeynSolver, billiard::Bi, k, dk; multithreaded::Bool=true) where {Bi<:AbsBilliard}
    error(_BEYN_NOT_IMPLEMENTED)
end