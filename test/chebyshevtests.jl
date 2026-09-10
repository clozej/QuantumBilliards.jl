# Regression tests for the BIM-solver-migration-plan Step 10 audit (Chebyshev-
# accelerated Hankel/Bessel evaluation). Covers two real fidelity bugs found
# and fixed while auditing the (previously undertested) `use_chebyshev=true`
# path, plus cross-checks that `use_chebyshev=true` agrees with the already-
# verified `use_chebyshev=false` reference values in solvertests.jl.
#
# `plan_h`/`panel_t`/`eval_h`/`hankel_z_chebyshev_cutoff` (solvers/chebyshev/
# bessels.jl) are not exported, so they're accessed here as
# `QuantumBilliards.<name>` — do not redefine them here.
using SpecialFunctions

# Unit regression test for the scalar `eval_h` mid-range-cutoff fix.
# `_cheb_geom_rminmax` floors a *batch* of Chebyshev plans' shared `rmin` using
# the LARGEST `|k|` in the batch (e.g. Beyn's contour nodes, which span a
# range of `|k|` around the circle). A plan built for a much SMALLER `|k|`
# sharing that same `rmin` can then have `|k*rmin| < hankel_z_chebyshev_cutoff`
# even though `r=rmin` still lies inside the "panelized" domain (`pidx != 0`,
# so the near-zero/`pidx==0` fallback never triggers). Before this fix,
# `eval_h` had no protection against this and evaluated the Chebyshev panel
# fit directly in the numerically sensitive near-origin band where the
# Hankel function's log singularity is hardest to interpolate; the fix adds
# the same mid-range direct-evaluation fallback that
# `QuantumBilliards-develop`'s combined H₀/H₁ evaluators already use.
@testset "Chebyshev eval_h mid-range cutoff fidelity (Step 10 fix)" begin
    k_big = 60.0
    k_small = 2.0
    rmin = QuantumBilliards.hankel_z_chebyshev_cutoff/k_big # floored for k_big, not k_small
    rmax = 5.0
    @test k_small*rmin<QuantumBilliards.hankel_z_chebyshev_cutoff # invariant violated for k_small

    plan1_small = QuantumBilliards.plan_h(1, 1, ComplexF64(k_small), rmin, rmax; npanels=200, M=8)
    pidx, t = QuantumBilliards.panel_t(plan1_small, rmin)
    @test pidx!=0 # r=rmin is inside the panelized domain, not the pidx==0 fallback

    val = QuantumBilliards.eval_h(plan1_small, pidx, t, rmin)
    ref = SpecialFunctions.besselh(1, 1, k_small*rmin)
    @test isapprox(val, ref; rtol=1e-9)
end

# Regression test for the `ExpandedBIMSolver` NaN bug (migration plan Step
# 10, implementation step 7, previously left uninvestigated). A plain
# `DoubleLayerPotentialSolver` kernel's k-derivative `A'(k)` is severely
# rank-deficient (roughly half its singular values are near machine epsilon,
# unlike `CombinedFieldIntegralEquationSolver`'s well-conditioned derivative,
# whose extra undifferentiated `i*S(k)` term adds full-rank structure).
# `LinearAlgebra.eigen`'s underlying `ggev` therefore reports some `NaN`
# generalized eigenvalues for the pencil `(A,A'(k))`. Before the
# `_argmin_finite` fix, plain `argmin(abs.(λ))` locked onto one of these
# spurious `NaN` eigenvalues (Julia's `min`-based reduction propagates `NaN`
# regardless of position), making `solve(::ExpandedBIMSolver,...)` return
# `NaN` for every trial `k` when wrapping a plain DLP kernel.
@testset "ExpandedBIMSolver (DLP kernel) - NaN regression - Full Triangle" begin
    billiard = QuantumBilliards.make_veech_right_triangle(5)
    kernel = DoubleLayerPotentialSolver(5.0; grading=GlobalCornerGrading())
    solver = ExpandedBIMSolver(kernel; use_chebyshev=false)
    k0 = 6.065090963607035 # Beyn-verified ground state (see solvertests.jl)
    pts = evaluate_points(solver, billiard, k0)
    kcorr, t0 = solve(solver, pts, k0)

    k_test = 6.065090994323157
    atol = 1e-6
    @test !isnan(kcorr)
    @test !isnan(t0)
    @test isapprox(kcorr, k_test; atol=atol)
    @test t0<1e-6
end

# solver: Beyn (contour-integral accelerated method, use_chebyshev=true)
# kernel: Double Layer Potential (Kress-corrected, GlobalCornerGrading)
# billiard: Triangle (full, un-reduced boundary, all edges SpecularReflection)
# symmetry: None
# functions to test: solve_wavenumber, compute_eigenstate
# Cross-checks against the already-verified `use_chebyshev=false` reference
# values in solvertests.jl's "Beyn (DLP kernel) - Full Triangle - Ground
# State" testset.
@testset "Beyn (DLP kernel, Chebyshev) - Full Triangle - Ground State" begin
    billiard = QuantumBilliards.make_veech_right_triangle(5)
    kernel = DoubleLayerPotentialSolver(5.0; grading=GlobalCornerGrading())
    solver_direct = BeynSolver(kernel; use_chebyshev=false)
    solver_cheb = BeynSolver(kernel; use_chebyshev=true)
    k0 = 6.1
    dk = 0.2
    kd, td = solve_wavenumber(solver_direct, billiard, k0, dk)
    kc, tc = solve_wavenumber(solver_cheb, billiard, k0, dk)
    state = compute_eigenstate(solver_cheb, billiard, kc)

    k_test = 6.065090994323157
    ten_test = 2.130269457223247e-9
    dim_test = 200
    atol = 1e-3
    @test isapprox(kc, k_test; atol=atol)
    @test isapprox(kc, kd; atol=1e-9) # Chebyshev agrees with direct evaluation to near machine precision
    @test isapprox(state.ten, ten_test; atol=atol)
    @test state.dim==dim_test
end

# solver: Beyn (contour-integral accelerated method, use_chebyshev=true)
# kernel: Combined Field Integral Equation (ungraded periodic, SmoothPeriodicGrading)
# billiard: Circle (PolarBilliard, no true corners)
# symmetry: None
# functions to test: solve_wavenumber, compute_eigenstate
@testset "Beyn (CFIE kernel, Chebyshev) - Circle - Ground State" begin
    billiard = BilliardGeometry.PolarBilliard([0.0, 0.0])
    kernel = CombinedFieldIntegralEquationSolver(5.0; grading=SmoothPeriodicGrading())
    solver_direct = BeynSolver(kernel; use_chebyshev=false)
    solver_cheb = BeynSolver(kernel; use_chebyshev=true)
    k0 = 2.4
    dk = 0.2
    kd, td = solve_wavenumber(solver_direct, billiard, k0, dk)
    kc, tc = solve_wavenumber(solver_cheb, billiard, k0, dk)
    state = compute_eigenstate(solver_cheb, billiard, kc)

    k_test = 2.404825557695771 # matches the first zero of J0, ≈ 2.404825557695772
    dim_test = 200
    atol = 1e-3
    @test isapprox(kc, k_test; atol=atol)
    @test isapprox(kc, kd; atol=1e-9) # Chebyshev agrees with direct evaluation to near machine precision
    @test state.dim==dim_test
end

# solver: ExpandedBIMSolver (local Taylor-expansion accelerated method, use_chebyshev=true)
# kernel: Double Layer Potential (Kress-corrected, GlobalCornerGrading)
# billiard: Triangle (full, un-reduced boundary, all edges SpecularReflection)
# symmetry: None
# functions to test: solve
# Also exercises the with-derivatives Chebyshev path
# (`_dlp_kernel_entry_with_derivatives_cheb`) together with the NaN-safe
# `solve` fix above (the plain-DLP NaN bug applies identically regardless of
# `use_chebyshev`, since both share the same `solve` eigenvalue selection).
@testset "ExpandedBIMSolver (DLP kernel, Chebyshev) - Full Triangle - Ground State" begin
    billiard = QuantumBilliards.make_veech_right_triangle(5)
    kernel = DoubleLayerPotentialSolver(5.0; grading=GlobalCornerGrading())
    solver_direct = ExpandedBIMSolver(kernel; use_chebyshev=false)
    solver_cheb = ExpandedBIMSolver(kernel; use_chebyshev=true)
    k0 = 6.065090963607035
    pts = evaluate_points(solver_direct, billiard, k0)
    kd, td = solve(solver_direct, pts, k0)
    kc, tc = solve(solver_cheb, pts, k0)

    k_test = 6.065090994323157
    atol = 1e-3
    @test !isnan(kc)
    @test !isnan(tc)
    @test isapprox(kc, k_test; atol=atol)
    @test isapprox(kc, kd; atol=1e-6) # Chebyshev agrees with direct evaluation
    @test tc<1e-6
end

# solver: ExpandedBIMSolver (local Taylor-expansion accelerated method, use_chebyshev=true)
# kernel: Combined Field Integral Equation (ungraded periodic, SmoothPeriodicGrading)
# billiard: Circle (PolarBilliard, no true corners)
# symmetry: None
# functions to test: solve
@testset "ExpandedBIMSolver (CFIE kernel, Chebyshev) - Circle - Ground State" begin
    billiard = BilliardGeometry.PolarBilliard([0.0, 0.0])
    kernel = CombinedFieldIntegralEquationSolver(5.0; grading=SmoothPeriodicGrading())
    solver_direct = ExpandedBIMSolver(kernel; use_chebyshev=false)
    solver_cheb = ExpandedBIMSolver(kernel; use_chebyshev=true)
    k0 = 2.4048255576957738
    pts = evaluate_points(solver_direct, billiard, k0)
    kd, td = solve(solver_direct, pts, k0)
    kc, tc = solve(solver_cheb, pts, k0)

    k_test = 2.404825557695771
    atol = 1e-3
    @test !isnan(kc)
    @test !isnan(tc)
    @test isapprox(kc, k_test; atol=atol)
    @test isapprox(kc, kd; atol=1e-9) # Chebyshev agrees with direct evaluation to near machine precision
    @test tc<1e-6
end
