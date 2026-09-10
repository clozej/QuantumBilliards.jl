################################################################################
# Chebyshev-accelerated CFIE Fredholm matrix assembly.
#
# New code mirroring solvers/sweepmethods/cfie.jl's
# `_cfie_fredholm_full!`/`_cfie_fredholm_reduced!`/`_cfie_kernel_entry` and
# solvers/acceleratedmethods/ebim.jl's `_cfie_kernel_entry_with_derivatives`,
# with `_bim_hankelh1`/`_bim_besselj` calls replaced by Chebyshev plan
# evaluations — see dlp.jl in this directory for the DLP counterpart and the
# module-level design note there.
################################################################################

# Full (unfolded) Chebyshev-accelerated CFIE Fredholm matrix A(k) = I - (D(k)+ikS(k)), value only.
function _cfie_fredholm_full_cheb!(F::AbstractMatrix{ComplexF64}, pts::BoundaryPoints{T}, Rmat::AbstractMatrix{T}, G::BoundaryGeomCache{T}, k::ComplexF64, plan0::ChebHankelPlanH, plan1::ChebHankelPlanH, planj0::ChebJPlan, planj1::ChebJPlan; multithreaded::Bool=true) where {T<:Real}
    invtwopi = inv(2*pi)
    αL1 = -k*invtwopi
    αL2 = im*k/2
    αM1 = -invtwopi
    αM2 = Complex{Float64}(0, 0.5)
    ik = im*k
    euler_over_pi = Float64(Base.MathConstants.eulergamma)/pi
    N = length(pts)
    fill!(F, zero(ComplexF64))
    @inbounds for i in 1:N
        si = Float64(G.speed[i])
        wi = pts.ws[i]
        dval = Complex{Float64}(wi*G.kappa[i], 0.0)
        m1 = αM1*si
        m2 = ((Complex{Float64}(0, 0.5)-euler_over_pi)-invtwopi*log((k^2/4)*si^2))*si
        sval = Complex{Float64}(Rmat[i,i]*m1, 0.0)+wi*m2
        F[i,i] = one(ComplexF64)-(dval+ik*sval)
    end
    @use_threads multithreading=(multithreaded && N>=32) for j in 2:N
        sj = Float64(G.speed[j])
        wj = pts.ws[j]
        @inbounds for i in 1:j-1
            si = Float64(G.speed[i])
            wi = pts.ws[i]
            r = Float64(G.R[i,j])
            invr = Float64(G.invR[i,j])
            lt = Float64(G.logterm[i,j])
            inn_ij = Float64(G.inner[i,j])
            inn_ji = Float64(G.inner[j,i])
            pidx_h, t_h = panel_t(plan1, r)
            pidx_j, t_j = panel_t(planj1, r)
            h1 = eval_h(plan1, pidx_h, t_h, r)
            h0 = eval_h(plan0, pidx_h, t_h, r)
            j1 = eval_j(planj1, pidx_j, t_j, r)
            j0 = eval_j(planj0, pidx_j, t_j, r)
            l1_ij = αL1*inn_ij*j1*invr
            l2_ij = αL2*inn_ij*h1*invr-l1_ij*lt
            dval_ij = Rmat[i,j]*l1_ij+wj*l2_ij
            m1_ij = αM1*j0*sj
            m2_ij = αM2*h0*sj-m1_ij*lt
            sval_ij = Rmat[i,j]*m1_ij+wj*m2_ij
            F[i,j] = -(dval_ij+ik*sval_ij)
            l1_ji = αL1*inn_ji*j1*invr
            l2_ji = αL2*inn_ji*h1*invr-l1_ji*lt
            dval_ji = Rmat[j,i]*l1_ji+wi*l2_ji
            m1_ji = αM1*j0*si
            m2_ji = αM2*h0*si-m1_ji*lt
            sval_ji = Rmat[j,i]*m1_ji+wi*m2_ji
            F[j,i] = -(dval_ji+ik*sval_ji)
        end
    end
    return F
end

# Single Kress-corrected (D(k)+ikS(k)) kernel entry, Chebyshev-evaluated (mirrors `_cfie_kernel_entry`).
@inline function _cfie_kernel_entry_cheb(pts::BoundaryPoints{T}, Rmat::AbstractMatrix{T}, G::BoundaryGeomCache{T}, k::ComplexF64, plan0::ChebHankelPlanH, plan1::ChebHankelPlanH, planj0::ChebJPlan, planj1::ChebJPlan, i::Int, j::Int) where {T<:Real}
    invtwopi = inv(2*pi)
    ik = im*k
    if i==j
        si = Float64(G.speed[i])
        wi = pts.ws[i]
        dval = Complex{Float64}(wi*G.kappa[i], 0.0)
        euler_over_pi = Float64(Base.MathConstants.eulergamma)/pi
        m1 = -invtwopi*si
        m2 = ((Complex{Float64}(0, 0.5)-euler_over_pi)-invtwopi*log((k^2/4)*si^2))*si
        sval = Complex{Float64}(Rmat[i,i]*m1, 0.0)+wi*m2
        return dval+ik*sval
    end
    r = Float64(G.R[i,j])
    invr = Float64(G.invR[i,j])
    lt = Float64(G.logterm[i,j])
    inn = Float64(G.inner[i,j])
    sj = Float64(G.speed[j])
    wj = pts.ws[j]
    pidx_h, t_h = panel_t(plan1, r)
    pidx_j, t_j = panel_t(planj1, r)
    h1 = eval_h(plan1, pidx_h, t_h, r)
    h0 = eval_h(plan0, pidx_h, t_h, r)
    j1 = eval_j(planj1, pidx_j, t_j, r)
    j0 = eval_j(planj0, pidx_j, t_j, r)
    αL1 = -k*invtwopi
    αL2 = im*k/2
    αM1 = -invtwopi
    αM2 = Complex{Float64}(0, 0.5)
    l1 = αL1*inn*j1*invr
    l2 = αL2*inn*h1*invr-l1*lt
    dval = Rmat[i,j]*l1+wj*l2
    m1 = αM1*j0*sj
    m2 = αM2*h0*sj-m1*lt
    sval = Rmat[i,j]*m1+wj*m2
    return dval+ik*sval
end

# Symmetry-reduced Chebyshev-accelerated CFIE Fredholm matrix, value only (mirrors `_cfie_fredholm_reduced!`).
function _cfie_fredholm_reduced_cheb!(F::AbstractMatrix{ComplexF64}, pts::BoundaryPoints{T}, Rmat::AbstractMatrix{T}, G::BoundaryGeomCache{T}, orbits::SymmetryOrbitMap{T}, k::ComplexF64, plan0::ChebHankelPlanH, plan1::ChebHankelPlanH, planj0::ChebJPlan, planj1::ChebJPlan; multithreaded::Bool=true) where {T<:Real}
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
                acc += phase[j]*_cfie_kernel_entry_cheb(pts, Rmat, G, k, plan0, plan1, planj0, planj1, i, j)
            end
            F[a,b] = -acc
        end
        F[b,b] += one(ComplexF64)
    end
    return F
end

# `(D(k)+ikS(k))` kernel entry plus its first two `k`-derivatives, Chebyshev-evaluated (mirrors `_cfie_kernel_entry_with_derivatives` in ebim.jl).
@inline function _cfie_kernel_entry_with_derivatives_cheb(pts::BoundaryPoints{T}, Rmat::AbstractMatrix{T}, G::BoundaryGeomCache{T}, k::ComplexF64, plan0::ChebHankelPlanH, plan1::ChebHankelPlanH, planj0::ChebJPlan, planj1::ChebJPlan, i::Int, j::Int) where {T<:Real}
    invtwopi = inv(2*pi)
    ik = im*k
    if i==j
        si = Float64(G.speed[i])
        wi = pts.ws[i]
        Dval = Complex{Float64}(wi*G.kappa[i], 0.0)
        euler_over_pi = Float64(Base.MathConstants.eulergamma)/pi
        m1 = -invtwopi*si
        m2 = ((Complex{Float64}(0, 0.5)-euler_over_pi)-invtwopi*log((k^2/4)*si^2))*si
        Sval = Complex{Float64}(Rmat[i,i]*m1, 0.0)+wi*m2
        val = Dval+ik*Sval
        dm2 = -si/(pi*k)
        ddm2 = si/(pi*k^2)
        dSval = wi*dm2
        ddSval = wi*ddm2
        dval = im*Sval+ik*dSval
        ddval = 2*im*dSval+ik*ddSval
        return val, dval, ddval
    end
    r = Float64(G.R[i,j])
    invr = Float64(G.invR[i,j])
    lt = Float64(G.logterm[i,j])
    inn = Float64(G.inner[i,j])
    sj = Float64(G.speed[j])
    wj = pts.ws[j]
    pidx_h, t_h = panel_t(plan1, r)
    pidx_j, t_j = panel_t(planj1, r)
    h1 = eval_h(plan1, pidx_h, t_h, r)
    h0 = eval_h(plan0, pidx_h, t_h, r)
    j1 = eval_j(planj1, pidx_j, t_j, r)
    j0 = eval_j(planj0, pidx_j, t_j, r)
    αL1 = -k*invtwopi
    αL2 = im*k/2
    αM1 = -invtwopi
    αM2 = Complex{Float64}(0, 0.5)
    l1, dl1, ddl1 = _ebim_lin1_deriv(αL1*inn, r, invr, j0, j1, k)
    l2a, dl2a, ddl2a = _ebim_lin1_deriv(αL2*inn, r, invr, h0, h1, k)
    l2 = l2a-l1*lt
    dl2 = dl2a-dl1*lt
    ddl2 = ddl2a-ddl1*lt
    Dval = Rmat[i,j]*l1+wj*l2
    dDval = Rmat[i,j]*dl1+wj*dl2
    ddDval = Rmat[i,j]*ddl1+wj*ddl2
    m1, dm1, ddm1 = _ebim_const0_deriv(αM1*sj, r, k, j0, j1)
    m2a, dm2a, ddm2a = _ebim_const0_deriv(αM2*sj, r, k, h0, h1)
    m2 = m2a-m1*lt
    dm2 = dm2a-dm1*lt
    ddm2 = ddm2a-ddm1*lt
    Sval = Rmat[i,j]*m1+wj*m2
    dSval = Rmat[i,j]*dm1+wj*dm2
    ddSval = Rmat[i,j]*ddm1+wj*ddm2
    val = Dval+ik*Sval
    dval = dDval+im*Sval+ik*dSval
    ddval = ddDval+2*im*dSval+ik*ddSval
    return val, dval, ddval
end
