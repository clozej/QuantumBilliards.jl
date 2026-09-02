################################################################################
# Determine suitable radial panel counts and polynomial degrees for the
# piecewise-Chebyshev approximation of Hankel and Bessel functions used by the
# boundary-integral solvers.
# Hankel interpolation: Hν⁽¹⁾(kr), ν=0,1,
# is restricted to radii satisfying |k|r >= hankel_z_chebyshev_cutoff,
# because the small-z Hankel region is evaluated directly. When several
# wavenumbers share one radial interpolation interval, the common lower bound is rmin_H=max(rmin_geom,max_j(zcut/|k_j|)).
# Bessel-J interpolation: Jν(kr), ν=0,1,
# is regular at r=0 and therefore uses the complete radial interval 0 <= r <= rmax.
################################################################################

# so the user does not need to know the actual :Symbol
@inline chebyshev_kind(::BoundaryIntegralMethod)=Val(:dlp)
@inline chebyshev_kind(::DLP)=Val(:dlp_kress)
@inline chebyshev_kind(::CFIE)=Val(:cfie_kress)

# common API so that e.g. Beyn and EBIM can call the same function to construct the boundary matrices for any solver type
@inline function construct_boundary_matrices!(Tbufs::Vector{Matrix{Complex{T}}},solver,pts,zj::AbstractVector{Complex{T}};multithreaded::Bool=true,use_chebyshev::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false) where {T<:Real}
    if use_chebyshev
        construct_matrices_chebyshev!(Tbufs,chebyshev_kind(solver),solver,pts,zj;multithreaded=multithreaded,n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,timeit=timeit)
    else
        construct_matrices!(Tbufs,solver,pts,zj;multithreaded=multithreaded)
    end
    return nothing
end

@inline function construct_boundary_matrices!(Tbufs::Vector{Matrix{Complex{T}}},dTbufs::Vector{Matrix{Complex{T}}},ddTbufs::Vector{Matrix{Complex{T}}},solver,pts,zj::AbstractVector{Complex{T}};multithreaded::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false) where {T<:Real}
    construct_matrices_chebyshev_with_derivatives!(Tbufs,dTbufs,ddTbufs,chebyshev_kind(solver),solver,pts,zj;multithreaded=multithreaded,n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,timeit=timeit)
    return nothing
end

"""
    _check_H0H1_errors!(err0,err1,plans0,plans1,ks,rs)

Check `H₀⁽¹⁾` and `H₁⁽¹⁾` Chebyshev plans against direct special-function
evaluation at sampled radii.

For each wavenumber `ks[j]`, the maximum absolute interpolation error over `rs`
is written to `err0[j]` and `err1[j]`.

## Arguments
- `err0`: Preallocated vector receiving maximum `H₀⁽¹⁾` errors.
- `err1`: Preallocated vector receiving maximum `H₁⁽¹⁾` errors.
- `plans0`: `ChebHankelPlanH` plans for `H₀⁽¹⁾`.
- `plans1`: `ChebHankelPlanH` plans for `H₁⁽¹⁾`.
- `ks`: Wavenumbers associated with the plans.
- `rs`: Sampled radii used for validation.

## Returns
`(err0,err1)`.
"""
function _check_H0H1_errors!(err0,err1,plans0,plans1,ks,rs)
    nz=length(ks)
    Threads.@threads for j in 1:nz
        e0=0.0
        e1=0.0
        k=ComplexF64(ks[j])
        plan0=plans0[j]
        plan1=plans1[j]
        @inbounds for r in rs
            p0,t0=panel_t(plan0,r)
            p1,t1=panel_t(plan1,r)
            z=k*r
            e0=max(e0,abs(eval_h(plan0,p0,t0,r)-SpecialFunctions.besselh(0,1,z)))
            e1=max(e1,abs(eval_h(plan1,p1,t1,r)-SpecialFunctions.besselh(1,1,z)))
        end
        err0[j]=e0
        err1[j]=e1
    end
    return err0,err1
end
"""
    _check_H1_errors!(err1,plans1,ks,rs)

Check `H₁⁽¹⁾` Chebyshev plans against direct special-function evaluation at
sampled radii.

For each wavenumber `ks[j]`, the maximum absolute interpolation error over `rs`
is written to `err1[j]`.

## Arguments
- `err1`: Preallocated vector receiving maximum `H₁⁽¹⁾` errors.
- `plans1`: `ChebHankelPlanH` plans for `H₁⁽¹⁾`.
- `ks`: Wavenumbers associated with the plans.
- `rs`: Sampled radii used for validation.

## Returns
`err1`.
"""
function _check_H1_errors!(err1,plans1,ks,rs)
    nz=length(ks)
    Threads.@threads for j in 1:nz
        e1=0.0
        k=ComplexF64(ks[j])
        plan1=plans1[j]
        @inbounds for r in rs
            p,t=panel_t(plan1,r)
            e1=max(e1,abs(eval_h(plan1,p,t,r)-SpecialFunctions.besselh(1,1,k*r)))
        end
        err1[j]=e1
    end
    return err1
end
"""
    _check_J0J1_errors!(err0,err1,plans0,plans1,ks,rs)

Check `J₀` and `J₁` Chebyshev plans against direct special-function evaluation
at sampled radii.

For each wavenumber `ks[j]`, the maximum absolute interpolation error over `rs`
is written to `err0[j]` and `err1[j]`.

## Arguments
- `err0`: Preallocated vector receiving maximum `J₀` errors.
- `err1`: Preallocated vector receiving maximum `J₁` errors.
- `plans0`: `ChebJPlan` plans for `J₀`.
- `plans1`: `ChebJPlan` plans for `J₁`.
- `ks`: Wavenumbers associated with the plans.
- `rs`: Sampled radii used for validation.

## Returns
`(err0,err1)`.
"""
function _check_J0J1_errors!(err0,err1,plans0,plans1,ks,rs)
    nz=length(ks)
    Threads.@threads for j in 1:nz
        e0=0.0
        e1=0.0
        k=ComplexF64(ks[j])
        plan0=plans0[j]
        plan1=plans1[j]
        @inbounds for r in rs
            p0,t0=panel_t(plan0,r)
            p1,t1=panel_t(plan1,r)
            z=k*r
            e0=max(e0,abs(eval_j(plan0,p0,t0,r)-SpecialFunctions.besselj(0,z)))
            e1=max(e1,abs(eval_j(plan1,p1,t1,r)-SpecialFunctions.besselj(1,z)))
        end
        err0[j]=e0
        err1[j]=e1
    end
    return err0,err1
end
########################################
######## ORDINARY BIM H1 TUNER #########
########################################
"""
    chebyshev_params(solver::BoundaryIntegralMethod,pts,zj;...)

Determine suitable Chebyshev panel count and polynomial degree for the
unscaled `H₁⁽¹⁾` kernel used by `BoundaryIntegralMethod`.

The tuner constructs `H₁⁽¹⁾` plans, validates them against
`SpecialFunctions.besselh`, and increases the radial panel count or polynomial
degree until the requested absolute tolerance is reached.
## Arguments
- `solver`: `BoundaryIntegralMethod` containing the active symmetry.
- `pts`: Boundary discretization.
- `zj`: Complex wavenumbers to be represented by the Chebyshev plans.

## Keyword Arguments
- `npanels_h_init`: Initial number of Hankel radial panels.
- `M_h_init`: Initial Hankel polynomial degree per panel.
- `npanels_j_init`: Unused compatibility keyword.
- `M_j_init`: Unused compatibility keyword.
- `tol`: Target maximum absolute interpolation error.
- `sampling_points`: Number of validation radii.
- `max_iter`: Maximum number of refinement iterations.
- `grow_panels`: Multiplicative growth factor for the panel count.
- `grow_M`: Additive increase in polynomial degree.
- `verbose`: Print tuning diagnostics.

## Returns
A tuple

    (n,M,0,0,plans1,err1)

where `n` is the final Hankel panel count, `M` is the final polynomial degree,
`plans1` contains the `H₁⁽¹⁾` plans, and `err1` contains the maximum validation
error for each wavenumber.
"""
function chebyshev_params(solver::BoundaryIntegralMethod,pts::BoundaryPoints{T},zj::AbstractVector{Complex{T}};npanels_h_init::Int=15_000,M_h_init::Int=5,npanels_j_init::Int=10_000,M_j_init::Int=5,tol::Real=1e-10,sampling_points::Int=50_000,max_iter::Int=20,grow_panels::Real=1.5,grow_M::Int=2,verbose::Bool=false) where {T<:Real}
    rmin_raw,rmax=estimate_rmin_rmax(pts,solver.symmetry)
    rmin_cheb=maximum(hankel_z_chebyshev_cutoff./abs.(zj))
    rmin_interp=max(rmin_raw,Float64(rmin_cheb))
    rmin_interp<rmax||throw(ArgumentError("Empty Hankel interpolation interval: rmin=$rmin_interp, rmax=$rmax"))
    verbose&&@info "Estimated BIM Chebyshev radial bounds" rmin_raw rmax rmin_cheb rmin_interp
    rs=collect(range(rmin_interp,rmax;length=sampling_points))
    nz=length(zj)
    n=npanels_h_init
    M=M_h_init
    plans1=Vector{ChebHankelPlanH}(undef,nz)
    err1=fill(Inf,nz)
    for it in 1:max_iter
        Threads.@threads for j in eachindex(zj)
            plans1[j]=plan_h(1,1,ComplexF64(zj[j]),rmin_interp,rmax;npanels=n,M=M)
        end
        _check_H1_errors!(err1,plans1,zj,rs)
        verbose&&@info "Worst BIM H1 | n_panels M" maximum(err1) n M
        all(<(tol),err1)&&return n,M,0,0,plans1,err1
        it%5==0 ? (M+=grow_M) : (n=ceil(Int,grow_panels*n))
    end
    @warn "BIM Chebyshev tuning did not reach tol=$tol after $max_iter iterations. Returning best effort."
    return n,M,0,0,plans1,err1
end
########################################
######## KRESS H0/H1/J0/J1 TUNER #######
########################################
"""
    chebyshev_params(solver::Union{CFIE,DLP},pts,zj;...)

Determine suitable Chebyshev panel counts and polynomial degrees for the
`H₀⁽¹⁾`, `H₁⁽¹⁾`, `J₀`, and `J₁` functions used by the Kress backends.

The tuner constructs all four plan families, validates them against direct
`SpecialFunctions` values, and independently refines the Hankel and Bessel-J
panel parameters until the requested tolerance is reached.

## Arguments
- `solver`: Kress-based CFIE or DLP solver.
- `pts`: One boundary discretization or a vector of boundary components.
- `zj`: Complex wavenumbers to be represented by the Chebyshev plans.

## Keyword Arguments
- `npanels_h_init`: Initial Hankel panel count.
- `M_h_init`: Initial Hankel polynomial degree.
- `npanels_j_init`: Initial Bessel-J panel count.
- `M_j_init`: Initial Bessel-J polynomial degree.
- `tol`: Target maximum absolute interpolation error.
- `sampling_points`: Number of validation radii in each radial interval.
- `max_iter`: Maximum number of refinement iterations.
- `grow_panels`: Multiplicative panel-count growth factor.
- `grow_M`: Additive polynomial-degree increment.
- `verbose`: Print tuning diagnostics.

## Returns
A tuple

    (nh,Mh,nj,Mj,plans0,plans1,plansj0,plansj1,
     errH0,errH1,errJ0,errJ1)

containing the final Hankel and Bessel-J interpolation parameters, the four plan
collections, and the corresponding per-wavenumber maximum validation errors.
"""
function chebyshev_params(solver::Union{CFIE,DLP},pts::Union{Vector{BoundaryPoints{T}},BoundaryPoints{T}},zj::AbstractVector{Complex{T}};npanels_h_init::Int=15_000,M_h_init::Int=5,npanels_j_init::Int=10_000,M_j_init::Int=5,tol::Real=1e-10,sampling_points::Int=50_000,max_iter::Int=20,grow_panels::Real=1.5,grow_M::Int=2,verbose::Bool=false) where {T<:Real}
    rmin_cheb=maximum(hankel_z_chebyshev_cutoff./abs.(zj))
    if solver isa CFIE
        ptsv=pts isa Vector ? pts : [pts]
        block_cache=build_cfie_kress_block_caches(solver,ptsv;npanels_h=16,M_h=4,npanels_j=16,M_j=4,rmin_cheb=rmin_cheb)
    else
        pts1=pts isa Vector ? (length(pts)==1 ? pts[1] : error("DLP-Kress expects one BoundaryPoints component.")) : pts
        block_cache=build_dlp_kress_block_cache(solver,pts1;npanels_h=16,M_h=4,npanels_j=16,M_j=4,rmin_cheb=rmin_cheb)
    end
    rmin_h=Float64(block_cache.rmin);rmax=Float64(block_cache.rmax)
    rmin_h<rmax||throw(ArgumentError("Empty Kress Hankel interpolation interval: rmin=$rmin_h, rmax=$rmax"))
    rsH=collect(range(rmin_h,rmax;length=sampling_points));rsJ=collect(range(0.0,rmax;length=sampling_points))
    nz=length(zj);nh=npanels_h_init;nj=npanels_j_init;Mh=M_h_init;Mj=M_j_init
    plans0=Vector{ChebHankelPlanH}(undef,nz);plans1=Vector{ChebHankelPlanH}(undef,nz)
    plansj0=Vector{ChebJPlan}(undef,nz);plansj1=Vector{ChebJPlan}(undef,nz)
    errH0=fill(Inf,nz);errH1=fill(Inf,nz);errJ0=fill(Inf,nz);errJ1=fill(Inf,nz)
    for it in 1:max_iter
        plans0,plans1,plansj0,plansj1=build_cfie_kress_plans(zj,rmin_h,rmax;npanels_h=nh,M_h=Mh,npanels_j=nj,M_j=Mj)
        _check_H0H1_errors!(errH0,errH1,plans0,plans1,zj,rsH)
        _check_J0J1_errors!(errJ0,errJ1,plansj0,plansj1,zj,rsJ)
        okH=all(<(tol),errH0)&&all(<(tol),errH1);okJ=all(<(tol),errJ0)&&all(<(tol),errJ1)
        verbose&&@info "Worst Kress H0 H1 J0 J1 | nh Mh nj Mj" maximum(errH0) maximum(errH1) maximum(errJ0) maximum(errJ1) nh Mh nj Mj
        okH&&okJ&&return nh,Mh,nj,Mj,plans0,plans1,plansj0,plansj1,errH0,errH1,errJ0,errJ1
        !okH&&(it%5==0 ? (Mh+=grow_M) : (nh=ceil(Int,grow_panels*nh)))
        !okJ&&(it%5==0 ? (Mj+=grow_M) : (nj=ceil(Int,grow_panels*nj)))
    end
    @warn "Kress Chebyshev tuning did not reach tol=$tol after $max_iter iterations. Returning best effort."
    return nh,Mh,nj,Mj,plans0,plans1,plansj0,plansj1,errH0,errH1,errJ0,errJ1
end