################################################################################
# Auto-tuning of Chebyshev radial panel counts / polynomial degrees for the
# Hankel/Bessel-J plans used by the accelerated BIM solvers (BeynSolver,
# ExpandedBIMSolver), bundled as a `ChebyshevConfig` (see
# solvers/acceleratedmethods/chebyshevconfig.jl).
#
# Adapted from QuantumBilliards-develop/src/chebyshev/chebyshev_optimal_panelization.jl's
# `chebyshev_params` tuning loop, restructured to:
#   - read/return a `ChebyshevConfig` instead of a bare
#     `(n_panels_h,M_h,n_panels_j,M_j)` tuple (see the Step 9.5 struct and the
#     migration plan's Step 10 implementation notes),
#   - determine the radial interpolation interval directly from a
#     `BoundaryGeomCache`'s pairwise distance matrix `G.R` (main's boundary
#     geometry cache already holds every pairwise distance needed; no
#     separate `estimate_rmin_rmax`-style symmetry-aware helper is required
#     here since we conservatively scan the full, unreduced pairwise-distance
#     matrix),
#   - validate against `SpecialFunctions.besselh`/`besselj` directly (the same
#     functions `_bim_hankelh1`/`_bim_besselj` dispatch to for `Complex` k, see
#     boundarygeomcache.jl).
################################################################################

# Off-diagonal pairwise-distance extrema of a boundary geometry cache, used as
# the Chebyshev radial interpolation interval (conservative: scans the full
# `G.R`, even when the caller only needs a symmetry-reduced subset of pairs).
# The lower bound is floored at `hankel_z_chebyshev_cutoff/max|k|` (mirroring
# `-develop`'s `rmin_cheb`), since Hankel functions are singular at `z=kr=0`:
# without this floor, a corner-graded discretization's very small nearest-
# neighbor spacing would force the interpolation domain deep into the
# near-singular region, where a uniform low-degree panelization needs an
# impractically large panel count to converge (radii below the floor are
# already handled correctly by `eval_h`/`eval_j`'s `pidx==0` direct/small-z
# series fallback, so nothing is lost by excluding them from the plan).
function _cheb_geom_rminmax(G::BoundaryGeomCache{T}, ks::Vector{ComplexF64}) where {T<:Real}
    Rm = G.R
    n = size(Rm, 1)
    rmin = Inf
    rmax = 0.0
    @inbounds for j in 1:n, i in 1:n
        i==j && continue
        r = Float64(Rm[i,j])
        r<rmin && (rmin = r)
        r>rmax && (rmax = r)
    end
    isfinite(rmin) && rmax>0.0 || throw(ArgumentError("Unable to determine a nonzero Chebyshev radial interpolation interval from the boundary geometry cache"))
    rmin_cheb = hankel_z_chebyshev_cutoff/maximum(abs, ks)
    rmin = max(rmin, rmin_cheb)
    rmin<rmax || throw(ArgumentError("Empty Chebyshev radial interpolation interval after flooring rmin at the near-zero cutoff: rmin=$rmin, rmax=$rmax"))
    return rmin, rmax
end

# Validate a batch of `ChebHankelPlanH`/`ChebJPlan` (one per wavenumber)
# against direct `SpecialFunctions` evaluation at sampled radii, writing the
# per-wavenumber maximum absolute error into `err`.
function _cheb_check_h_errors!(err::Vector{Float64}, plans::Vector{ChebHankelPlanH}, rs::Vector{Float64})
    Threads.@threads for j in eachindex(plans)
        pl = plans[j]
        e = 0.0
        @inbounds for r in rs
            p, t = panel_t(pl, r)
            e = max(e, abs(eval_h(pl, p, t, r)-SpecialFunctions.besselh(pl.ν, pl.κ, pl.k*r)))
        end
        err[j] = e
    end
    return err
end
function _cheb_check_j_errors!(err::Vector{Float64}, plans::Vector{ChebJPlan}, rs::Vector{Float64})
    Threads.@threads for j in eachindex(plans)
        pl = plans[j]
        e = 0.0
        @inbounds for r in rs
            p, t = panel_t(pl, r)
            e = max(e, abs(eval_j(pl, p, t, r)-SpecialFunctions.besselj(pl.ν, pl.k*r)))
        end
        err[j] = e
    end
    return err
end

"""
    tune_dlp_cheb_plans(rmin, rmax, ks, cfg::ChebyshevConfig) → (plans1, plansj1, cfg_used)

Builds (and, unless `cfg.param_strategy===:manual`, auto-tunes) the `H₁^(1)`
and `J₁` Chebyshev plans, one per complex wavenumber in `ks`, needed by the
value-only (Beyn) DLP Chebyshev assembly (dlp.jl in this directory).

When `cfg.param_strategy===:manual`, the plans are built once from
`cfg.n_panels_h`/`cfg.M_h`/`cfg.n_panels_j`/`cfg.M_j` with no validation loop.
Otherwise the panel count/polynomial degree are grown
(`cfg.grow_panels`/`cfg.grow_M`, alternating every 5th iteration as in
`-develop`) until the worst-case error over `cfg.sampling_points` validation
radii drops below `cfg.tol`, or `cfg.max_iter` is reached (a `@warn` is
raised and the best-effort plans are returned).

Returns the plans together with a `ChebyshevConfig` recording the
panel/degree parameters actually used (so callers with `param_strategy in
(:global,:segment)` can reuse the tuned config without re-tuning, see
`compute_spectrum(::ExpandedBIMSolver,...)`).
"""
function tune_dlp_cheb_plans(rmin::Float64, rmax::Float64, ks::Vector{ComplexF64}, cfg::ChebyshevConfig{T}) where {T<:Real}
    nz = length(ks)
    if cfg.param_strategy===:manual
        plans1 = Vector{ChebHankelPlanH}(undef, nz)
        plansj1 = Vector{ChebJPlan}(undef, nz)
        @inbounds Threads.@threads for j in 1:nz
            plans1[j] = plan_h(1, 1, ks[j], rmin, rmax; npanels=cfg.n_panels_h, M=cfg.M_h)
            plansj1[j] = plan_j(1, ks[j], rmin, rmax; npanels=cfg.n_panels_j, M=cfg.M_j)
        end
        return plans1, plansj1, cfg
    end
    rs = collect(range(rmin, rmax; length=cfg.sampling_points))
    nh, Mh = cfg.n_panels_h, cfg.M_h
    nj, Mj = cfg.n_panels_j, cfg.M_j
    tol = Float64(cfg.tol)
    plans1 = Vector{ChebHankelPlanH}(undef, nz)
    plansj1 = Vector{ChebJPlan}(undef, nz)
    errh = fill(Inf, nz)
    errj = fill(Inf, nz)
    for it in 1:cfg.max_iter
        @inbounds Threads.@threads for j in 1:nz
            plans1[j] = plan_h(1, 1, ks[j], rmin, rmax; npanels=nh, M=Mh)
            plansj1[j] = plan_j(1, ks[j], rmin, rmax; npanels=nj, M=Mj)
        end
        _cheb_check_h_errors!(errh, plans1, rs)
        _cheb_check_j_errors!(errj, plansj1, rs)
        okh = all(<(tol), errh)
        okj = all(<(tol), errj)
        if okh && okj
            cfg_used = ChebyshevConfig(T; n_panels_h=nh, M_h=Mh, n_panels_j=nj, M_j=Mj, tol=cfg.tol, max_iter=cfg.max_iter, sampling_points=cfg.sampling_points, grow_panels=cfg.grow_panels, grow_M=cfg.grow_M, param_strategy=cfg.param_strategy)
            return plans1, plansj1, cfg_used
        end
        okh || (it%5==0 ? (Mh += cfg.grow_M) : (nh = ceil(Int, cfg.grow_panels*nh)))
        okj || (it%5==0 ? (Mj += cfg.grow_M) : (nj = ceil(Int, cfg.grow_panels*nj)))
    end
    @warn "DLP Chebyshev tuning did not reach tol=$tol after $(cfg.max_iter) iterations. Using best-effort panels/degrees." maximum(errh) maximum(errj) nh Mh nj Mj
    cfg_used = ChebyshevConfig(T; n_panels_h=nh, M_h=Mh, n_panels_j=nj, M_j=Mj, tol=cfg.tol, max_iter=cfg.max_iter, sampling_points=cfg.sampling_points, grow_panels=cfg.grow_panels, grow_M=cfg.grow_M, param_strategy=cfg.param_strategy)
    return plans1, plansj1, cfg_used
end

"""
    tune_cfie_cheb_plans(rmin, rmax, ks, cfg::ChebyshevConfig) → (plans0, plans1, plansj0, plansj1, cfg_used)

Same as [`tune_dlp_cheb_plans`](@ref), but additionally builds/tunes the `H₀^(1)`
and `J₀` plans needed by the CFIE's single-layer `S(k)` term (see cfie.jl in
this directory).
"""
function tune_cfie_cheb_plans(rmin::Float64, rmax::Float64, ks::Vector{ComplexF64}, cfg::ChebyshevConfig{T}) where {T<:Real}
    nz = length(ks)
    if cfg.param_strategy===:manual
        plans0 = Vector{ChebHankelPlanH}(undef, nz)
        plans1 = Vector{ChebHankelPlanH}(undef, nz)
        plansj0 = Vector{ChebJPlan}(undef, nz)
        plansj1 = Vector{ChebJPlan}(undef, nz)
        @inbounds Threads.@threads for j in 1:nz
            plans0[j] = plan_h(0, 1, ks[j], rmin, rmax; npanels=cfg.n_panels_h, M=cfg.M_h)
            plans1[j] = plan_h(1, 1, ks[j], rmin, rmax; npanels=cfg.n_panels_h, M=cfg.M_h)
            plansj0[j] = plan_j(0, ks[j], rmin, rmax; npanels=cfg.n_panels_j, M=cfg.M_j)
            plansj1[j] = plan_j(1, ks[j], rmin, rmax; npanels=cfg.n_panels_j, M=cfg.M_j)
        end
        return plans0, plans1, plansj0, plansj1, cfg
    end
    rs = collect(range(rmin, rmax; length=cfg.sampling_points))
    nh, Mh = cfg.n_panels_h, cfg.M_h
    nj, Mj = cfg.n_panels_j, cfg.M_j
    tol = Float64(cfg.tol)
    plans0 = Vector{ChebHankelPlanH}(undef, nz)
    plans1 = Vector{ChebHankelPlanH}(undef, nz)
    plansj0 = Vector{ChebJPlan}(undef, nz)
    plansj1 = Vector{ChebJPlan}(undef, nz)
    errh0 = fill(Inf, nz); errh1 = fill(Inf, nz)
    errj0 = fill(Inf, nz); errj1 = fill(Inf, nz)
    for it in 1:cfg.max_iter
        @inbounds Threads.@threads for j in 1:nz
            plans0[j] = plan_h(0, 1, ks[j], rmin, rmax; npanels=nh, M=Mh)
            plans1[j] = plan_h(1, 1, ks[j], rmin, rmax; npanels=nh, M=Mh)
            plansj0[j] = plan_j(0, ks[j], rmin, rmax; npanels=nj, M=Mj)
            plansj1[j] = plan_j(1, ks[j], rmin, rmax; npanels=nj, M=Mj)
        end
        _cheb_check_h_errors!(errh0, plans0, rs)
        _cheb_check_h_errors!(errh1, plans1, rs)
        _cheb_check_j_errors!(errj0, plansj0, rs)
        _cheb_check_j_errors!(errj1, plansj1, rs)
        okh = all(<(tol), errh0) && all(<(tol), errh1)
        okj = all(<(tol), errj0) && all(<(tol), errj1)
        if okh && okj
            cfg_used = ChebyshevConfig(T; n_panels_h=nh, M_h=Mh, n_panels_j=nj, M_j=Mj, tol=cfg.tol, max_iter=cfg.max_iter, sampling_points=cfg.sampling_points, grow_panels=cfg.grow_panels, grow_M=cfg.grow_M, param_strategy=cfg.param_strategy)
            return plans0, plans1, plansj0, plansj1, cfg_used
        end
        okh || (it%5==0 ? (Mh += cfg.grow_M) : (nh = ceil(Int, cfg.grow_panels*nh)))
        okj || (it%5==0 ? (Mj += cfg.grow_M) : (nj = ceil(Int, cfg.grow_panels*nj)))
    end
    @warn "CFIE Chebyshev tuning did not reach tol=$tol after $(cfg.max_iter) iterations. Using best-effort panels/degrees." maximum(errh0) maximum(errh1) maximum(errj0) maximum(errj1) nh Mh nj Mj
    cfg_used = ChebyshevConfig(T; n_panels_h=nh, M_h=Mh, n_panels_j=nj, M_j=Mj, tol=cfg.tol, max_iter=cfg.max_iter, sampling_points=cfg.sampling_points, grow_panels=cfg.grow_panels, grow_M=cfg.grow_M, param_strategy=cfg.param_strategy)
    return plans0, plans1, plansj0, plansj1, cfg_used
end
