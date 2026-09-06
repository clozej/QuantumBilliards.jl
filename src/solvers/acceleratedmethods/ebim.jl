"""
    ExpandedBIMSolver{T,K} <: AcceleratedBIMSolver

`ExpandedBIMSolver` is a concrete [`AcceleratedBIMSolver`](@ref) implementing
the expanded boundary integral method (EBIM): a single locally-corrected root
of the nonlinear eigenproblem `A(k)v = 0` obtained from a second-order local
Taylor expansion of `A(k)`.

## Description
`ExpandedBIMSolver` wraps an inner `kernel::SweepBIMSolver` (a
[`DoubleLayerPotentialSolver`](@ref), [`CombinedFieldIntegralEquationSolver`](@ref)
or [`CompositeBIMSolver`](@ref)) supplying the Fredholm operator and its first
two `k`-derivatives,

    A(k+ε) = A(k) + ε A'(k) + (1/2) ε² A''(k) + O(ε³).

The generalized eigenproblem `A(k)v = λ A'(k)v` gives the first-order root
correction `ε₁ = -λ`; with the corresponding left generalized eigenvector `u`,
the second-order correction is `ε₂ = -(1/2) ε₁² [u'A''(k)v]/[u'A'(k)v]`, giving
the corrected wavenumber `k_corr = k + ε₁ + ε₂` (see [`construct_matrices`](@ref),
[`solve`](@ref)). Unlike [`BeynSolver`](@ref), EBIM produces a single locally
refined root per call rather than every root in a window.

## Attributes
* `kernel`: The wrapped [`SweepBIMSolver`](@ref) supplying `A(k)`, `A'(k)`, `A''(k)`.
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
- [`solve_wavenumber`](@ref)
- [`solve_spectrum`](@ref)

!!! note "Migration status"
    API scaffold only (struct, constructor and method signatures) ported from
    `QuantumBilliards-develop`'s `EBIMSolver`. The matrix-assembly/root-correction
    bodies are not yet implemented; every method below raises an `error` until
    Step 2 of the migration plan lands.
"""
struct ExpandedBIMSolver{T<:Real,K<:SweepBIMSolver} <: AcceleratedBIMSolver
    kernel::K
    use_chebyshev::Bool
    n_panels_h::Int
    M_h::Int
    n_panels_j::Int
    M_j::Int
end

"""
    ExpandedBIMSolver(kernel::K; use_chebyshev::Bool = true, n_panels_h::Int = 15000, M_h::Int = 5, n_panels_j::Int = 10000, M_j::Int = 5) where {K<:SweepBIMSolver} → solver::ExpandedBIMSolver

Constructs an [`ExpandedBIMSolver`](@ref) wrapping the boundary-integral
`kernel`.

## Arguments
* `kernel`: The [`SweepBIMSolver`](@ref) supplying the Fredholm operator `A(k)` and its `k`-derivatives.

## Keyword arguments
* `use_chebyshev::Bool = true`: Whether to use Chebyshev-accelerated kernel evaluation.
* `n_panels_h::Int = 15000`: Hankel-function Chebyshev panel count.
* `M_h::Int = 5`: Hankel-function Chebyshev polynomial degree.
* `n_panels_j::Int = 10000`: Bessel-J-function Chebyshev panel count.
* `M_j::Int = 5`: Bessel-J-function Chebyshev polynomial degree.

## Returns
* `solver`: An [`ExpandedBIMSolver`](@ref) instance.
"""
function ExpandedBIMSolver(kernel::K; use_chebyshev::Bool=true,
                            n_panels_h::Int=15000, M_h::Int=5, n_panels_j::Int=10000, M_j::Int=5) where {K<:SweepBIMSolver}
    T = _bim_numeric_type(kernel)
    return ExpandedBIMSolver{T,K}(kernel, use_chebyshev, n_panels_h, M_h, n_panels_j, M_j)
end

const _EBIM_NOT_IMPLEMENTED = "ExpandedBIMSolver matrix assembly is not yet implemented (API scaffold only, see the QuantumBilliardsTests migration plan)."

"""
    construct_matrices(solver::ExpandedBIMSolver, pts::BoundaryPoints, k; multithreaded::Bool = true) → (A::Matrix, dA::Matrix, ddA::Matrix)

Assembles the Fredholm matrix `A(k)` and its first two `k`-derivatives
`A'(k)`, `A''(k)`.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function construct_matrices(solver::ExpandedBIMSolver, pts::BoundaryPoints, k; multithreaded::Bool=true)
    error(_EBIM_NOT_IMPLEMENTED)
end

"""
    solve(solver::ExpandedBIMSolver, pts::BoundaryPoints, k; multithreaded::Bool = true) → (k_corr::Real, t0::Real)

Computes the second-order locally-corrected root `k_corr` near `k` and its
tension `t0`.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function solve(solver::ExpandedBIMSolver, pts::BoundaryPoints, k; multithreaded::Bool=true)
    error(_EBIM_NOT_IMPLEMENTED)
end

"""
    solve_wavenumber(solver::ExpandedBIMSolver, billiard::Bi, k, dk; multithreaded::Bool = true) where {Bi<:AbsBilliard} → (k0::Real, t0::Real)

Computes the second-order locally-corrected root nearest `k` (`dk` is retained
for API parity with [`solve_wavenumber(::BeynSolver, ...)`](@ref) but is
unused by the local expansion).

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function solve_wavenumber(solver::ExpandedBIMSolver, billiard::Bi, k, dk; multithreaded::Bool=true) where {Bi<:AbsBilliard}
    error(_EBIM_NOT_IMPLEMENTED)
end

"""
    solve_spectrum(solver::ExpandedBIMSolver, billiard::Bi, k, dk; multithreaded::Bool = true) where {Bi<:AbsBilliard} → (ks::Vector, ts::Vector)

Computes the second-order locally-corrected roots for every wavenumber in a
sweep, one call per target `k`.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function solve_spectrum(solver::ExpandedBIMSolver, billiard::Bi, k, dk; multithreaded::Bool=true) where {Bi<:AbsBilliard}
    error(_EBIM_NOT_IMPLEMENTED)
end