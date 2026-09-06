"""
    CombinedFieldIntegralEquationSolver{T,G,Sy} <: SweepBIMSolver

`CombinedFieldIntegralEquationSolver` is a concrete [`SweepBIMSolver`](@ref)
implementing the (optionally Kress-corrected) combined-field boundary integral
method for computing quantum billiard spectra from the coupled Helmholtz
double-layer/single-layer Fredholm operator.

## Description
The assembled Fredholm operator is

    A(k) = I - (D(k) + i k S(k)),

where `D(k)` is the Nyström discretization of the Helmholtz double-layer
operator and `S(k)` of the single-layer operator. The combined-field
coupling removes the spurious interior-nullspace problem that a pure
double-layer formulation ([`DoubleLayerPotentialSolver`](@ref)) can suffer on
some geometries. The tension at a fixed wavenumber `k` is a function of the
smallest singular value / nullspace residual of `A(k)` (see
[`construct_matrices`](@ref), [`solve`](@ref)). The boundary discretization
strategy is controlled by `grading`, see [`BoundaryGrading`](@ref)
(`CFIE` additionally supports [`CornerGrading`](@ref), a single-curve
parametric-corner variant that [`DoubleLayerPotentialSolver`](@ref) does not
need); an optional discrete `symmetry` folds the discretization onto a
fundamental domain via a [`BilliardGeometry.SymmetryOrbitMap`](@ref).

## Attributes
* `pts_scaling_factor`: Vector of scaling factors, one per fundamental boundary curve, used to determine the number of boundary sampling points.
* `min_pts`: Minimum number of boundary sampling points per component.
* `grading`: [`BoundaryGrading`](@ref) strategy used to discretize the boundary.
* `symmetry`: Optional `AbsSymmetry` used to fold the discretization onto a fundamental domain.
* `eps`: Relative tolerance used to determine the tension from the smallest singular value / nullspace residual.

## API
The following functions can be evaluated for this type:
- [`evaluate_points`](@ref)
- [`boundary_matrix_size`](@ref)
- [`construct_matrices`](@ref)
- [`solve`](@ref)
- [`solve_vect`](@ref)
- [`solve_wavenumber`](@ref)
- [`k_sweep`](@ref)

!!! note "Migration status"
    API scaffold only (struct, constructor and method signatures) ported from
    `QuantumBilliards-develop`'s `CFIE_kress`/`CFIE_kress_corners`/
    `CFIE_kress_global_corners`. The matrix-assembly bodies are not yet
    implemented; every method below raises an `error` until Step 2 of the
    migration plan lands.
"""
struct CombinedFieldIntegralEquationSolver{T<:Real,G<:BoundaryGrading,Sy<:Union{AbsSymmetry,Nothing}} <: SweepBIMSolver
    pts_scaling_factor::Vector{T}
    min_pts::Int64
    grading::G
    symmetry::Sy
    eps::T
end

"""
    CombinedFieldIntegralEquationSolver(pts_scaling_factor::Union{T,Vector{T}}; min_pts::Int = 200, grading::BoundaryGrading = SmoothPeriodicGrading(), symmetry::Union{Nothing,AbsSymmetry} = nothing, eps::T = T(1e-15)) where {T<:Real} → solver::CombinedFieldIntegralEquationSolver{T}

Constructs a [`CombinedFieldIntegralEquationSolver`](@ref).

## Arguments
* `pts_scaling_factor`: Scaling factor, or vector thereof (one per fundamental boundary curve), used to determine the number of boundary sampling points.

## Keyword arguments
* `min_pts::Int = 200`: Minimum number of boundary sampling points per component.
* `grading::BoundaryGrading = SmoothPeriodicGrading()`: Boundary discretization/grading strategy, see [`BoundaryGrading`](@ref).
* `symmetry::Union{Nothing,AbsSymmetry} = nothing`: Optional discrete symmetry used to fold the discretization onto a fundamental domain.
* `eps::T = T(1e-15)`: Relative tolerance used to determine the tension.

## Returns
* `solver`: A [`CombinedFieldIntegralEquationSolver{T}`](@ref) instance.
"""
function CombinedFieldIntegralEquationSolver(pts_scaling_factor::Union{T,Vector{T}}; min_pts::Int=200,
                                              grading::BoundaryGrading=SmoothPeriodicGrading(),
                                              symmetry::Union{Nothing,AbsSymmetry}=nothing,
                                              eps::T=T(1e-15)) where {T<:Real}
    bs = pts_scaling_factor isa T ? [pts_scaling_factor] : pts_scaling_factor
    return CombinedFieldIntegralEquationSolver{T,typeof(grading),typeof(symmetry)}(bs, min_pts, grading, symmetry, eps)
end

_bim_numeric_type(::CombinedFieldIntegralEquationSolver{T}) where {T} = T

const _CFIE_NOT_IMPLEMENTED = "CombinedFieldIntegralEquationSolver matrix assembly is not yet implemented (API scaffold only, see the QuantumBilliardsTests migration plan)."

"""
    evaluate_points(solver::CombinedFieldIntegralEquationSolver, billiard::Bi, k) where {Bi<:AbsBilliard} → pts::BoundaryPoints

Samples the boundary of `billiard` according to `solver.grading`, producing the
boundary discretization needed to assemble the combined-field Fredholm matrix
in [`construct_matrices`](@ref).

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function evaluate_points(solver::CombinedFieldIntegralEquationSolver, billiard::Bi, k) where {Bi<:AbsBilliard}
    error(_CFIE_NOT_IMPLEMENTED)
end

"""
    boundary_matrix_size(solver::CombinedFieldIntegralEquationSolver, pts::BoundaryPoints) → N::Int

Returns the dimension of the assembled Fredholm matrix, accounting for any
symmetry-orbit folding onto a fundamental domain.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function boundary_matrix_size(solver::CombinedFieldIntegralEquationSolver, pts::BoundaryPoints)
    error(_CFIE_NOT_IMPLEMENTED)
end

"""
    construct_matrices(solver::CombinedFieldIntegralEquationSolver, pts::BoundaryPoints, k; multithreaded::Bool = true) → A::Matrix{Complex}

Assembles the combined-field Fredholm matrix `A(k) = I - (D(k) + i k S(k))`.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function construct_matrices(solver::CombinedFieldIntegralEquationSolver, pts::BoundaryPoints, k; multithreaded::Bool=true)
    error(_CFIE_NOT_IMPLEMENTED)
end

"""
    solve(solver::CombinedFieldIntegralEquationSolver, pts::BoundaryPoints, k; multithreaded::Bool = true, use_krylov::Bool = true) → t::Real

Computes the combined-field tension at wavenumber `k`, defined from the
smallest singular value / Krylov nullspace residual of `A(k)` (see
[`construct_matrices`](@ref)).

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function solve(solver::CombinedFieldIntegralEquationSolver, pts::BoundaryPoints, k; multithreaded::Bool=true, use_krylov::Bool=true)
    error(_CFIE_NOT_IMPLEMENTED)
end

"""
    solve_vect(solver::CombinedFieldIntegralEquationSolver, pts::BoundaryPoints, k; multithreaded::Bool = true) → (t::Real, x::Vector)

Computes the combined-field tension and the associated boundary density
eigenvector at wavenumber `k`.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function solve_vect(solver::CombinedFieldIntegralEquationSolver, pts::BoundaryPoints, k; multithreaded::Bool=true)
    error(_CFIE_NOT_IMPLEMENTED)
end