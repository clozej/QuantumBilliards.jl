# BIM solver (DLP/CFIE/CompositeBIM) migration notes

- Step 3 (`DoubleLayerPotentialSolver`, solvers/sweepmethods/dlp.jl) is
  fully implemented and numerically verified (not a stub) as of 2026-09-07.
  Do not re-implement; if bodies look like stubs again, something regressed.
- `_global_t_to_segment_u`/`_eval_composite_geom_global_t` live in
  BilliardGeometry.jl/src/geometry/boundarycomponents.jl (ported, unexported).
- `BIMEigenstate{K,T,S,Bi}` struct has `k::K`/`vec::Vector{K}` sharing the
  SAME type param K. Since BIM boundary densities are complex, `k` must be
  promoted to `Complex{T}` in the `BIMEigenstate(k,vec,ten,solver,billiard)`
  constructor (`kK = eltype(vec)(k)`), and `eps = set_precision(real(vec[1]))`
  (NOT `set_precision(vec[1])`) so `eps::T` matches `ten::T` (real). Otherwise
  you get a MethodError (eps ends up Complex, ten stays real -> no match).
- `KrylovKit` is a genuine dependency of both QuantumBilliards.jl and
  QuantumBilliards-develop (confirmed in both Project.toml `[deps]`), and
  -develop's own `DLP_kress` `solve`/`solve_vect` do use Krylov nullspace
  methods (`smallest_nullvec_krylov!`, `@svd_or_det_solve` dispatching on
  `use_krylov`) — so `dlp.jl`'s `KrylovKit.svdsolve(A,1,:SR)` in `solve`/
  `solve_vect` is a justified match, not a new/guessed dependency.

## Important gotcha: `get_boundary_curves` is NOT the full physical boundary
  for every billiard fixture — it filters to only `SpecularReflection`-typed
  curves in the `fundamental_domain`, dropping `Transparent`/
  `ReflectionSymmetry`/`QuantumSolverIgnore` curves. This means:
  - `StadiumBilliard(half_width)` (`BilliardGeometry.jl/src/geometry/billiards/stadium.jl`)
    is built with a D2-symmetry-reduced quarter fundamental domain (one
    `CircleSegment` + one `LineSegment`, total perimeter ≈ quarter of the
    real stadium). `get_boundary_curves(StadiumBilliard(0.5))` returns just
    that quarter, artificially closed — closing that loop creates a FAKE
    90°-angle "corner" at the seam (confirmed via `_component_corner_locations`
    returning `[0.0]` with junction angle ≈ π/2). This is fine for
    `RealPlaneWaves`-based basis solvers (plane waves already satisfy the
    reflection BCs across the omitted symmetry lines) but is **wrong** for
    any BIM solver (DLP/CFIE) — do not use `StadiumBilliard` as a BIM test
    fixture without first building/using its true full un-reduced boundary.
  - `make_triangle_and_basis(gamma, chi; edge_i)` (`utils/billiardutils.jl`)
    sets 2 of 3 triangle edges to `QuantumSolverIgnore()` (only `edge_i` is
    real) — it's built for the Veech corner-adapted basis method, and its
    `get_boundary_curves` returns only ONE edge. Also not usable directly
    for BIM.
  - For a genuine full closed billiard boundary (BIM fixture), use
    `BilliardGeometry.TriangleBilliard(gamma, chi)` **without** the `bcs`
    keyword (default is all 3 edges `SpecularReflection()`, no symmetry) —
    same physical shape/spectrum as `make_triangle_and_basis(gamma,chi)`'s
    triangle (translation doesn't change the spectrum), so it cross-validates
    against `test/solvertests.jl`'s Veech-triangle reference k0's. Confirmed:
    low-k ground state k0≈6.06509 matches VS reference 6.065082959967892 to
    ~8e-6 with `DoubleLayerPotentialSolver(5.0; grading=GlobalCornerGrading())`.
  - `PolarBilliard(coef)` (`r(φ)=1+Σaₙcos+Σbₙsin`, `coef=[b1,a1,b2,a2,...]`)
    with `coef=[0.0,0.0]` gives a genuine full-boundary unit circle (no
    symmetry, no corners) — great `SmoothPeriodicGrading` analytic-eigenvalue
    fixture (Bessel zeros): DLP matched J0/J1 first zeros to ~1e-8/~1e-9.
    `coef=[0,0,0,0.3]` (i.e. `r=1+0.3cos(2φ)`) gives a genuinely D2-symmetric
    FULL boundary (no pre-reduction) — good fixture for verifying
    `solver.symmetry=BilliardGeometry.XAxisReflection()` folding: confirmed
    the symmetric-sector sweep reproduces symmetric minima exactly and
    correctly OMITS antisymmetric ones present in the unrestricted sweep.

## `solve_wavenumber`/`Optim.optimize` gotcha at high k
  At high k (e.g. k~110 on a triangle of modest area), Weyl-law eigenvalue
  spacing can be smaller than a "reasonable-looking" search window `dk=0.1-0.2`.
  `solve_wavenumber`'s bounded `Optim.optimize` can converge to the WRONG
  nearby local tension minimum if two eigenvalues both fall inside the
  window (confirmed: window [109.99,110.19] contains minima at both k≈110.09
  AND k≈110.14, and Brent's method picked the k≈110.14 one first). Not a bug
  in the solver — narrow the window (`dk` small enough to isolate a single
  local minimum, verified via a `k_sweep` scan first) when cross-validating
  a specific reference eigenvalue at high k.

## QBPlotting.jl environment
  `QBPlotting.jl` has no `Manifest.toml` in this workspace (never
  instantiated) — `using QBPlotting` fails on missing `StaticArrays` etc.
  This is pre-existing/unrelated to any BIM work; `get_errors` (static) is
  the only available signal for QBPlotting.jl changes unless someone runs
  `Pkg.instantiate()` there (not done automatically — avoid unless asked,
  it's a broad environment change).
