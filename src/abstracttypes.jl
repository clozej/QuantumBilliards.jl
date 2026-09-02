abstract type AbsPoints end
"""
CoordinateSystem

`CoordinateSystem` is the abstract supertype for local coordinate frames used
to evaluate basis functions and their gradients.

## Description
A `CoordinateSystem` bundles an origin, a rotation angle and the corresponding
affine maps (and their inverses) needed to transform points between the global
Cartesian frame and a local frame in which basis functions are naturally
expressed, e.g. a Cartesian frame aligned with a symmetry axis, or a polar
frame centered at a corner. Concrete subtypes are [`CartesianCS`](@ref) and
[`PolarCS`](@ref).
"""
abstract type CoordinateSystem end

"""
AbsBasis

`AbsBasis` is the abstract supertype for all basis representations used to
approximate eigenstates and boundary solutions of a quantum billiard.

## Description
A concrete subtype of `AbsBasis` stores the parameters defining a family of
basis functions (e.g. plane waves or corner-adapted Fourier-Bessel functions)
at a fixed dimension `dim` and wavenumber. Solvers combine an `AbsBasis` with
a boundary discretization to construct the matrices used to determine
eigenvalues and eigenvectors. Concrete subtypes include
[`RealPlaneWaves`](@ref) and [`CornerAdaptedFourierBessel`](@ref).

## API
Every concrete subtype of `AbsBasis` is expected to implement:
- `resize_basis`
- `basis_fun`
- `gradient`
- `basis_and_gradient`
"""
abstract type AbsBasis end

"""
AbsSolver

`AbsSolver` is the abstract supertype for all algorithms that determine
quantum billiard eigenvalues (wavenumbers) and eigenvectors from a boundary
discretization.

## Description
Concrete solvers hold the parameters controlling boundary sampling and basis
dimension scaling (e.g. `dim_scaling_factor`, `pts_scaling_factor`, `sampler`,
`min_dim`, `min_pts`) as well as the numerical tolerance `eps` used when
filtering the generalized eigenvalue problem. `AbsSolver` has two direct
concrete-algorithm branches, [`SweepSolver`](@ref) and
[`AcceleratedSolver`](@ref), distinguished by whether the spectrum is scanned
one wavenumber at a time or obtained in windows via a single diagonalization.

## API
The following functions can be evaluated for any `AbsSolver`:
- `evaluate_points`
- [`adjust_scaling_and_samplers`](@ref)
- `compute_spectrum`
"""
abstract type AbsSolver end

"""
SweepSolver <: AbsSolver

`SweepSolver` is the abstract supertype for solvers that locate quantum
billiard eigenvalues by sweeping over a range of individual wavenumbers and
minimizing a tension function at each one.

## Description
At each wavenumber `k`, a `SweepSolver` constructs matrices from a boundary
quadrature and solves a generalized eigenvalue problem whose smallest
eigenvalue defines a tension quantifying how well the boundary condition is
satisfied; scanning this tension over a range of wavenumbers locates the
billiard's eigenvalues. The concrete implementation is
[`DecompositionMethodSolver`](@ref).

## API
The following functions can be evaluated for any `SweepSolver`:
- `construct_matrices`
- `solve`
- `solve_vect`
- [`solve_wavenumber`](@ref)
- [`k_sweep`](@ref)
- [`compute_eigenstate`](@ref)
"""
abstract type SweepSolver <: AbsSolver end

"""
AcceleratedSolver <: AbsSolver

`AcceleratedSolver` is the abstract supertype for solvers that recover every
eigenvalue within a wavenumber window `dk` of a target wavenumber `k` from a
single diagonalization.

## Description
An `AcceleratedSolver` constructs a generalized eigenvalue problem whose
spectrum, restricted to the window around `k`, approximates the tensions of
all billiard eigenstates in that window, avoiding the need to scan
wavenumber-by-wavenumber as [`SweepSolver`](@ref) does. The concrete
implementation is [`VerginiSaracenoSolver`](@ref).

## API
The following functions can be evaluated for any `AcceleratedSolver`:
- `construct_matrices`
- `solve`
- `solve_vectors`
- [`solve_wavenumber`](@ref)
- `solve_spectrum`
- [`compute_eigenstate`](@ref)
"""
abstract type AcceleratedSolver <: AbsSolver end

"""
AbsState

`AbsState` is the abstract supertype for all representations of a quantum
billiard state expressed through a coefficient vector `vec` at a wavenumber
`k`.

## Description
Every concrete subtype of `AbsState` stores at least a wavenumber `k`, the
wavenumber `k_basis` at which the associated basis coefficients `vec` were
evaluated, the dimension `dim` of `vec`, and a numerical precision threshold
`eps` (see `set_precision`) below which coefficients are treated as zero.
[`StationaryState`](@ref) is the branch of `AbsState` for states expressed in
a genuine [`AbsBasis`](@ref); [`GaussianRandomState`](@ref) is a direct
`AbsState` subtype not tied to any specific basis, used as a random-wave
reference ensemble.

## API
The following functions can be evaluated for any `AbsState`:
- `boundary_function`
- `momentum_function`
- `wavefunction`
- `husimi_function`
"""
abstract type AbsState end

"""
StationaryState <: AbsState

`StationaryState` is the abstract supertype for states expressed as expansion
coefficients in a concrete [`AbsBasis`](@ref).

## Description
In addition to the fields required by [`AbsState`](@ref), a `StationaryState`
stores the `basis` (resized/evaluated at `k_basis`) in which its coefficient
vector `vec` is expressed, so that the state can be evaluated pointwise via
the basis' evaluation functions. Concrete subtypes are [`Eigenstate`](@ref),
a numerically computed billiard eigenstate, and [`BasisState`](@ref), a
single unmixed basis function viewed as a stationary state.
"""
abstract type StationaryState <: AbsState end
