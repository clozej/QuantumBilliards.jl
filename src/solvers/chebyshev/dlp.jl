################################################################################
# Chebyshev-accelerated DLP Fredholm matrix assembly.
#
# New code (not a verbatim port of QuantumBilliards-develop's
# chebyshev_dlp.jl/chebyshev_dlp_kress.jl): main's `DoubleLayerPotentialSolver`
# already unifies every `-develop` DLP grading variant behind a single
# `_dlp_fredholm_full!`/`_dlp_fredholm_reduced!`/`_dlp_kernel_entry`
# assembly (solvers/sweepmethods/dlp.jl) and a single
# `_dlp_kernel_entry_with_derivatives` (solvers/acceleratedmethods/ebim.jl),
# so the Chebyshev acceleration is written as drop-in replacements for those
# exact functions (same Kress-split kernel algebra, same call sites), rather
# than as separate per-grading files. Every arithmetic expression below is
# copied unchanged from its direct-evaluation counterpart; only the
# `_bim_hankelh1`/`_bim_besselj` calls are replaced with Chebyshev plan
# evaluations (`eval_h`/`eval_j` from bessels.jl in this directory).
#
# `plan1`/`planj1` (value-only) and `plan0`/`plan1`/`planj0`/`planj1`
# (with-derivatives) are built once per `construct_matrices` call by
# `tune_dlp_cheb_plans` (optimalpanelization.jl) and reused across the full
# O(N²) pairwise loop.
################################################################################

# Full (unfolded) Chebyshev-accelerated Fredholm matrix F(k) = I - D(k), value only.
function _dlp_fredholm_full_cheb!(F::AbstractMatrix{ComplexF64}, pts::BoundaryPoints{T}, Rmat::AbstractMatrix{T}, G::BoundaryGeomCache{T}, k::ComplexF64, plan1::ChebHankelPlanH, planj1::ChebJPlan; multithreaded::Bool=true) where {T<:Real}
    invtwopi = inv(2*pi)
    αL1 = -k*invtwopi
    αL2 = im*k/2
    N = length(pts)
    fill!(F, zero(ComplexF64))
    @inbounds for i in 1:N
        F[i,i] = one(ComplexF64)-Complex{Float64}(pts.ws[i]*G.kappa[i], 0.0)
    end
    @use_threads multithreading=(multithreaded && N>=32) for j in 2:N
        @inbounds for i in 1:j-1
            r = Float64(G.R[i,j])
            invr = Float64(G.invR[i,j])
            lt = Float64(G.logterm[i,j])
            inn_ij = Float64(G.inner[i,j])
            inn_ji = Float64(G.inner[j,i])
            pidx_h, t_h = panel_t(plan1, r)
            pidx_j, t_j = panel_t(planj1, r)
            h1 = eval_h(plan1, pidx_h, t_h, r)
            j1 = eval_j(planj1, pidx_j, t_j, r)
            l1_ij = αL1*inn_ij*j1*invr
            l2_ij = αL2*inn_ij*h1*invr-l1_ij*lt
            F[i,j] = -(Rmat[i,j]*l1_ij+pts.ws[j]*l2_ij)
            l1_ji = αL1*inn_ji*j1*invr
            l2_ji = αL2*inn_ji*h1*invr-l1_ji*lt
            F[j,i] = -(Rmat[j,i]*l1_ji+pts.ws[i]*l2_ji)
        end
    end
    return F
end

# Single Kress-corrected DLP kernel entry D[i,j] at full-boundary indices, Chebyshev-evaluated (mirrors `_dlp_kernel_entry`).
@inline function _dlp_kernel_entry_cheb(pts::BoundaryPoints{T}, Rmat::AbstractMatrix{T}, G::BoundaryGeomCache{T}, k::ComplexF64, plan1::ChebHankelPlanH, planj1::ChebJPlan, i::Int, j::Int) where {T<:Real}
    i==j && return Complex{Float64}(pts.ws[i]*G.kappa[i], 0.0)
    invtwopi = inv(2*pi)
    r = Float64(G.R[i,j])
    invr = Float64(G.invR[i,j])
    lt = Float64(G.logterm[i,j])
    inn = Float64(G.inner[i,j])
    pidx_h, t_h = panel_t(plan1, r)
    pidx_j, t_j = panel_t(planj1, r)
    h1 = eval_h(plan1, pidx_h, t_h, r)
    j1 = eval_j(planj1, pidx_j, t_j, r)
    l1 = -k*invtwopi*inn*j1*invr
    l2 = im*k/2*inn*h1*invr-l1*lt
    return Rmat[i,j]*l1+pts.ws[j]*l2
end

# Symmetry-reduced Chebyshev-accelerated Fredholm matrix, value only (mirrors `_dlp_fredholm_reduced!`).
function _dlp_fredholm_reduced_cheb!(F::AbstractMatrix{ComplexF64}, pts::BoundaryPoints{T}, Rmat::AbstractMatrix{T}, G::BoundaryGeomCache{T}, orbits::SymmetryOrbitMap{T}, k::ComplexF64, plan1::ChebHankelPlanH, planj1::ChebJPlan; multithreaded::Bool=true) where {T<:Real}
    m = fundamental_size(orbits)
    N = length(orbits)
    fund = orbits.fundamental_indices
    orbit_of = orbits.orbit_of
    phase = orbits.phase
    images = [Int[] for _ in 1:m]
    @inbounds for j in 1:N
        push!(images[orbit_of[j]], j)
    end
    fill!(F, zero(ComplexF64))
    @use_threads multithreading=(multithreaded && m>=32) for b in 1:m
        @inbounds for a in 1:m
            i = fund[a]
            acc = zero(ComplexF64)
            for j in images[b]
                acc += phase[j]*_dlp_kernel_entry_cheb(pts, Rmat, G, k, plan1, planj1, i, j)
            end
            F[a,b] = -acc
        end
        F[b,b] += one(ComplexF64)
    end
    return F
end

# `D(k)` kernel entry plus its first two `k`-derivatives, Chebyshev-evaluated
# (mirrors `_dlp_kernel_entry_with_derivatives` in ebim.jl; `plan0`/`plan1`
# share panelization, as do `planj0`/`planj1`, since they are tuned together
# by `tune_dlp_cheb_plans`/`tune_cfie_cheb_plans`).
@inline function _dlp_kernel_entry_with_derivatives_cheb(pts::BoundaryPoints{T}, Rmat::AbstractMatrix{T}, G::BoundaryGeomCache{T}, k::ComplexF64, plan0::ChebHankelPlanH, plan1::ChebHankelPlanH, planj0::ChebJPlan, planj1::ChebJPlan, i::Int, j::Int) where {T<:Real}
    if i==j
        return Complex{Float64}(pts.ws[i]*G.kappa[i], 0.0), zero(ComplexF64), zero(ComplexF64)
    end
    invtwopi = inv(2*pi)
    r = Float64(G.R[i,j])
    invr = Float64(G.invR[i,j])
    lt = Float64(G.logterm[i,j])
    inn = Float64(G.inner[i,j])
    pidx_h, t_h = panel_t(plan1, r)
    pidx_j, t_j = panel_t(planj1, r)
    h1 = eval_h(plan1, pidx_h, t_h, r)
    h0 = eval_h(plan0, pidx_h, t_h, r)
    j1 = eval_j(planj1, pidx_j, t_j, r)
    j0 = eval_j(planj0, pidx_j, t_j, r)
    αL1 = -k*invtwopi
    αL2 = im*k/2
    l1, dl1, ddl1 = _ebim_lin1_deriv(αL1*inn, r, invr, j0, j1, k)
    l2a, dl2a, ddl2a = _ebim_lin1_deriv(αL2*inn, r, invr, h0, h1, k)
    l2 = l2a-l1*lt
    dl2 = dl2a-dl1*lt
    ddl2 = ddl2a-ddl1*lt
    val = Rmat[i,j]*l1+pts.ws[j]*l2
    dval = Rmat[i,j]*dl1+pts.ws[j]*dl2
    ddval = Rmat[i,j]*ddl1+pts.ws[j]*ddl2
    return val, dval, ddval
end
