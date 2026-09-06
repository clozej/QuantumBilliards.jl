"""
    CompositeBIMSolver{T,CS,Sy} <: SweepBIMSolver

`CompositeBIMSolver` is a concrete [`SweepBIMSolver`](@ref) for multiply
connected geometries whose connected boundary components require different
boundary-integral discretizations.

## Description
`CompositeBIMSolver` assigns one existing [`SweepBIMSolver`](@ref) component
solver (a [`DoubleLayerPotentialSolver`](@ref) or
[`CombinedFieldIntegralEquationSolver`](@ref)) to each connected physical
boundary component, assembling one globally coupled Fredholm operator. The
component solvers control only the discretization and same-component
Kress/grading quadrature of their own boundary component; inter-component
interactions are evaluated with ordinary Nyström quadrature. The first
component solver is interpreted as the outer boundary; the remaining
component solvers (`2:end`) are interpreted as holes and are
orientation-reversed after discretization.

## Attributes
* `component_solvers`: Tuple with one [`SweepBIMSolver`](@ref) per connected boundary component, outer boundary first.
* `symmetry`: Optional discrete symmetry shared by every component solver.

## API
The following functions can be evaluated for this type:
- [`evaluate_points`](@ref)
- [`construct_matrices`](@ref)
- [`solve`](@ref)
- [`solve_vect`](@ref)
- [`solve_wavenumber`](@ref)
- [`k_sweep`](@ref)

!!! note "Migration status"
    API scaffold only (struct, constructor and method signatures), generalizing
    `QuantumBilliards-develop`'s `CFIE_kress_composite_solver` to both the
    [`DoubleLayerPotentialSolver`](@ref) and
    [`CombinedFieldIntegralEquationSolver`](@ref) families. The matrix-assembly
    bodies are not yet implemented; every method below raises an `error` until
    Step 2 of the migration plan lands.
"""
struct CompositeBIMSolver{T<:Real,CS<:Tuple,Sy<:Union{AbsSymmetry,Nothing}} <: SweepBIMSolver
    component_solvers::CS
    symmetry::Sy
end

"""
    CompositeBIMSolver(component_solvers::SweepBIMSolver...) → solver::CompositeBIMSolver

Constructs a [`CompositeBIMSolver`](@ref) from one component solver per
connected boundary component, outer boundary first.

## Arguments
* `component_solvers`: One [`SweepBIMSolver`](@ref) instance per connected boundary component. Every component solver must share the same `symmetry`.

## Returns
* `solver`: A [`CompositeBIMSolver`](@ref) instance.
"""
function CompositeBIMSolver(component_solvers::Vararg{SweepBIMSolver})
    isempty(component_solvers) && throw(ArgumentError("CompositeBIMSolver requires at least one component solver"))
    symmetry = component_solvers[1].symmetry
    all(cs -> cs.symmetry == symmetry, component_solvers) || throw(ArgumentError("All component solvers passed to CompositeBIMSolver must share the same symmetry"))
    T = _bim_numeric_type(component_solvers[1])
    return CompositeBIMSolver{T,typeof(component_solvers),typeof(symmetry)}(component_solvers, symmetry)
end

_bim_numeric_type(solver::CompositeBIMSolver{T}) where {T} = T

const _COMPOSITE_BIM_NOT_IMPLEMENTED = "CompositeBIMSolver matrix assembly is not yet implemented (API scaffold only, see the QuantumBilliardsTests migration plan)."

"""
    evaluate_points(solver::CompositeBIMSolver, billiard::Bi, k) where {Bi<:AbsBilliard} → pts::BoundaryPoints

Samples every connected boundary component of `billiard` with its assigned
component solver, concatenating the results (holes orientation-reversed) into
one composite [`BoundaryPoints`](@ref) discretization.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function evaluate_points(solver::CompositeBIMSolver, billiard::Bi, k) where {Bi<:AbsBilliard}
    error(_COMPOSITE_BIM_NOT_IMPLEMENTED)
end

"""
    construct_matrices(solver::CompositeBIMSolver, pts::BoundaryPoints, k; multithreaded::Bool = true) → A::Matrix{Complex}

Assembles the globally coupled composite Fredholm matrix `A(k)`.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function construct_matrices(solver::CompositeBIMSolver, pts::BoundaryPoints, k; multithreaded::Bool=true)
    error(_COMPOSITE_BIM_NOT_IMPLEMENTED)
end

"""
    solve(solver::CompositeBIMSolver, pts::BoundaryPoints, k; multithreaded::Bool = true) → t::Real

Computes the composite tension at wavenumber `k`.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function solve(solver::CompositeBIMSolver, pts::BoundaryPoints, k; multithreaded::Bool=true)
    error(_COMPOSITE_BIM_NOT_IMPLEMENTED)
end

"""
    solve_vect(solver::CompositeBIMSolver, pts::BoundaryPoints, k; multithreaded::Bool = true) → (t::Real, x::Vector)

Computes the composite tension and the associated boundary density eigenvector
at wavenumber `k`.

!!! note "Migration status"
    Not yet implemented; raises an `error`.
"""
function solve_vect(solver::CompositeBIMSolver, pts::BoundaryPoints, k; multithreaded::Bool=true)
    error(_COMPOSITE_BIM_NOT_IMPLEMENTED)
end
