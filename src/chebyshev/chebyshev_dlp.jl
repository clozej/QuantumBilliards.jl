################################################################################
# DLP Chebyshev matrix assembly
#
# Fast construction of the doubled 2D Helmholtz double-layer boundary-integral
# operator for complex wavenumbers, with optional exact symmetry reduction and
# first two k-derivatives.
#
# DLP convention:
#
#     K(k;x,y)=(i k/2)H₁⁽¹⁾(kr)c(x,y),
#
# where
#
#     r=|x-y|,
#     c(x,y)=((x-y)⋅n_y)/r.
#
# The smooth diagonal limit is
#
#     K(x,x;k)=-κ(x)/(2π).
#
# The first two wavenumber derivatives are
#
#     K'(k;x,y)=(i k r/2)H₀⁽¹⁾(kr)c(x,y),
#
#     K''(k;x,y)
#       =(i/2)[rH₀⁽¹⁾(kr)-k r²H₁⁽¹⁾(kr)]c(x,y).
#
# Equivalently,
#
#     K''(k;x,y)
#       =(i/(2k))[(2-(kr)²)H₁⁽¹⁾(kr)-krH₂⁽¹⁾(kr)]c(x,y).
#
# The Fredholm operator is
#
#     A(k)=I-K(k)W,
#
# with
#
#     W=diag(ds).
#
# Hence
#
#     A'(k)=-K'(k)W,
#     A''(k)=-K''(k)W.
#
# SYMMETRY
#
# The complete physical boundary is always discretized. When symmetry is
# active, an exact SymmetryOrbitMap folds the complete source sum onto the
# fundamental boundary indices.
#
# For reduced source column b let
#
#     j=Ifund[b],
#     q_l=fund_to_full[l,b],
#     χ_l=fund_to_scale[l,b].
#
# The raw reduced kernel is defined as
#
#     K_red[a,b]
#       =Σ_l χ_l K(i_a,q_l) ds[q_l]/ds[j].
#
# Therefore multiplying reduced column b by ds[j] gives
#
#     K_red[a,b]ds[j]
#       =Σ_l χ_l K(i_a,q_l)ds[q_l],
#
# which is exactly the full-boundary Nyström source sum restricted to the
# selected irrep.
################################################################################
#
# Reference:
#   R. Kress, "Boundary Integral Equations in Time-Harmonic Acoustic
#   Scattering," Mathl. Comput. Modelling 15(3-5), 229-243 (1991).
################################################################################

const _DLP_CHEB_INV_TWO_PI=inv(2*pi)

########################################
######## DIRECT COMPLEX-k DLP ##########
########################################

@inline function _dlp_regular_complex_entry(xi::T,yi::T,xj::T,yj::T,nxj::T,nyj::T,k::Complex{T}) where {T<:Real}
    dx=xi-xj
    dy=yi-yj
    r=hypot(dx,dy)
    iszero(r)&&throw(ArgumentError("Regular complex-k DLP kernel received coincident target and source points"))
    c=(nxj*dx+nyj*dy)/r
    return Complex{T}(0,one(T)/2)*k*c*SpecialFunctions.hankelh1(1,k*r)
end

@inline function _h0_h1_one_k_at_r!(h0vals::Vector{ComplexF64},h1vals::Vector{ComplexF64},plans0::Vector{ChebHankelPlanH},plans1::Vector{ChebHankelPlanH},p,t,r)
    h0_h1_multi_ks_at_r!(h0vals,h1vals,plans0,plans1,p,t,r)
    return h0vals[1],h1vals[1]
end

"""
    compute_kernel_matrix_complex_k!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},k::Complex{T};multithreaded::Bool=true) where {T<:Real} → Nothing

Assemble the complete raw doubled Helmholtz DLP kernel at one complex
wavenumber without symmetry reduction.

## Arguments
* `K::AbstractMatrix{Complex{T}}`: Preallocated `N×N` raw kernel matrix.
* `bp::BoundaryPoints{T}`: Full boundary discretization.
* `k::Complex{T}`: Complex wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Enable threaded kernel assembly.

## Returns
* `nothing`.
"""
function compute_kernel_matrix_complex_k!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},k::Complex{T};multithreaded::Bool=true) where {T<:Real}
    N=length(bp)
    @assert size(K)==(N,N)
    fill!(K,zero(Complex{T}))
    xy=bp.xy
    nrm=bp.normal
    κ=bp.curvature
    tol2=eps(T)^2
    pref=Complex{T}(0,one(T)/2)*k
    @use_threads multithreading=multithreaded for i in 1:N
        point_i=xy[i]
        normal_i=nrm[i]
        @inbounds begin
            K[i,i]=-Complex{T}(κ[i]*_DLP_CHEB_INV_TWO_PI,zero(T))
            for j in 1:i-1
                point_j=xy[j]
                normal_j=nrm[j]
                dx=point_i[1]-point_j[1]
                dy=point_i[2]-point_j[2]
                d2=muladd(dx,dx,dy*dy)
                d2<=tol2&&throw(ArgumentError("Distinct full-boundary DLP nodes coincide at indices ($i,$j)"))
                r=sqrt(d2)
                h=pref*SpecialFunctions.hankelh1(1,k*r)
                invr=inv(r)
                K[i,j]=(normal_j[1]*dx+normal_j[2]*dy)*invr*h
                K[j,i]=(normal_i[1]*(-dx)+normal_i[2]*(-dy))*invr*h
            end
        end
    end
    return nothing
end

"""
    compute_kernel_matrix_complex_k!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},k::Complex{T};multithreaded::Bool=true) where {T<:Real} → Nothing

Assemble the exact symmetry-reduced raw DLP kernel at one complex wavenumber.

For representative source `j=Ifund[b]`,

    K_red[a,b]=Σ_l χ_l K(i_a,q_l)ds[q_l]/ds[j].

## Arguments
* `K::AbstractMatrix{Complex{T}}`: Preallocated reduced kernel matrix.
* `bp::BoundaryPoints{T}`: Full boundary discretization.
* `orbits::SymmetryOrbitMap{T}`: Exact symmetry-orbit map.
* `k::Complex{T}`: Complex wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Enable threaded reduced assembly.

## Returns
* `nothing`.
"""
function compute_kernel_matrix_complex_k!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},k::Complex{T};multithreaded::Bool=true) where {T<:Real}
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    ng=orbit_size(orbits)
    @assert size(K)==(m,m)
    fill!(K,zero(Complex{T}))
    xy=bp.xy
    nrm=bp.normal
    κ=bp.curvature
    ds=bp.ds
    tol2=eps(T)^2
    pref=Complex{T}(0,one(T)/2)*k
    @use_threads multithreading=multithreaded for b in 1:m
        j=Ifund[b]
        wj=ds[j]
        @inbounds for a in 1:m
            i=Ifund[a]
            point_i=xy[i]
            val=zero(Complex{T})
            for l in 1:ng
                q=orbits.fund_to_full[l,b]
                scale=orbits.fund_to_scale[l,b]
                weight_ratio=ds[q]/wj
                if l==1&&i==j
                    val+=scale*(-Complex{T}(κ[i]*_DLP_CHEB_INV_TWO_PI,zero(T)))*weight_ratio
                    continue
                end
                point_q=xy[q]
                dx=point_i[1]-point_q[1]
                dy=point_i[2]-point_q[2]
                d2=muladd(dx,dx,dy*dy)
                d2<=tol2&&throw(ArgumentError("A nonidentity DLP symmetry image coincides with reduced target index $i"))
                r=sqrt(d2)
                normal_q=nrm[q]
                c=(normal_q[1]*dx+normal_q[2]*dy)/r
                val+=scale*weight_ratio*c*pref*SpecialFunctions.hankelh1(1,k*r)
            end
            K[a,b]=val
        end
    end
    return nothing
end

"""
    compute_kernel_matrix_complex_k!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},symmetry,k::Complex{T};multithreaded::Bool=true) where {T<:Real} → Nothing

Assemble a full or symmetry-reduced raw complex-k DLP kernel.

## Arguments
* `K::AbstractMatrix{Complex{T}}`: Destination matrix.
* `bp::BoundaryPoints{T}`: Full boundary discretization.
* `symmetry`: Active symmetry or `nothing`.
* `k::Complex{T}`: Complex wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Enable threaded assembly.

## Returns
* `nothing`.
"""
function compute_kernel_matrix_complex_k!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},symmetry,k::Complex{T};multithreaded::Bool=true) where {T<:Real}
    isnothing(symmetry)&&return compute_kernel_matrix_complex_k!(K,bp,k;multithreaded=multithreaded)
    orbits=symmetry_index_orbits(T,bp,symmetry)
    return compute_kernel_matrix_complex_k!(K,bp,orbits,k;multithreaded=multithreaded)
end

"""
    fredholm_matrix_complex_k!(A::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},symmetry,k::Complex{T};multithreaded::Bool=true) where {T<:Real} → Nothing

Assemble the full or symmetry-reduced BIM Fredholm matrix

    A(k)=I-K(k)W

using direct Hankel evaluation.

## Arguments
* `A::AbstractMatrix{Complex{T}}`: Destination Fredholm matrix.
* `bp::BoundaryPoints{T}`: Full boundary discretization.
* `symmetry`: Active symmetry or `nothing`.
* `k::Complex{T}`: Complex wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Enable threaded kernel assembly.

## Returns
* `nothing`.
"""
function fredholm_matrix_complex_k!(A::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},symmetry,k::Complex{T};multithreaded::Bool=true) where {T<:Real}
    if isnothing(symmetry)
        compute_kernel_matrix_complex_k!(A,bp,k;multithreaded=multithreaded)
        assemble_fredholm_matrices!(A,bp)
    else
        orbits=symmetry_index_orbits(T,bp,symmetry)
        compute_kernel_matrix_complex_k!(A,bp,orbits,k;multithreaded=multithreaded)
        assemble_fredholm_matrices!(A,bp,orbits)
    end
    return nothing
end

########################################
########### KERNEL FORMULAS ############
########################################

@inline function _kernel_triplet_from_hankels(k,r,H0,H1)
    kr=k*r
    hK=0.5im*k*H1
    hdK=0.5im*kr*H0
    hddK=0.5im*(r*H0-k*r*r*H1)
    return hK,hdK,hddK
end

#############################
#### ACCUMULATION HELPERS ####
#############################

@inline function _accum_dlp_default_nosym!(K::AbstractMatrix{Complex{T}},i::Int,j::Int,nxi::T,nyi::T,nxj::T,nyj::T,dx::T,dy::T,invr::T,h) where {T<:Real}
    @inbounds K[i,j]+=(nxj*dx+nyj*dy)*invr*h
    if i!=j
        @inbounds K[j,i]+=(nxi*(-dx)+nyi*(-dy))*invr*h
    end
    return nothing
end

@inline function _accum_dlp_triplet_nosym!(K::AbstractMatrix{Complex{T}},dK::AbstractMatrix{Complex{T}},ddK::AbstractMatrix{Complex{T}},i::Int,j::Int,nxi::T,nyi::T,nxj::T,nyj::T,dx::T,dy::T,invr::T,hK,hdK,hddK) where {T<:Real}
    c=(nxj*dx+nyj*dy)*invr
    @inbounds begin
        K[i,j]+=c*hK
        dK[i,j]+=c*hdK
        ddK[i,j]+=c*hddK
    end
    if i!=j
        c=(nxi*(-dx)+nyi*(-dy))*invr
        @inbounds begin
            K[j,i]+=c*hK
            dK[j,i]+=c*hdK
            ddK[j,i]+=c*hddK
        end
    end
    return nothing
end

########################################
######## DERIVATIVE WORKSPACE ##########
########################################

"""
    DLPDerivChebWorkspace

Thread-local Hankel buffers used by the derivative Chebyshev pathway.

## Attributes
* `h0_tls::Vector{Vector{ComplexF64}}`: Thread-local `H₀⁽¹⁾` values.
* `h1_tls::Vector{Vector{ComplexF64}}`: Thread-local `H₁⁽¹⁾` values.
"""
struct DLPDerivChebWorkspace
    h0_tls::Vector{Vector{ComplexF64}}
    h1_tls::Vector{Vector{ComplexF64}}
end

"""
    DLPDerivChebWorkspace(Mk::Int,nth::Int=Threads.nthreads()) → DLPDerivChebWorkspace

Allocate thread-local Hankel buffers used by derivative DLP assembly.

`nth` must be at least `Threads.nthreads()` because threaded matrix assembly
indexes these buffers using `Threads.threadid()`.

## Arguments
* `Mk::Int`: Number of simultaneously evaluated wavenumbers.
* `nth::Int`: Number of thread-local buffer sets.

## Returns
* `ws::DLPDerivChebWorkspace`: Allocated derivative workspace.
"""
function DLPDerivChebWorkspace(Mk::Int,nth::Int=Threads.nthreads())
    nth>=Threads.nthreads()||throw(ArgumentError("nth must be at least Threads.nthreads()=$(Threads.nthreads()); received nth=$nth"))
    return DLPDerivChebWorkspace([Vector{ComplexF64}(undef,Mk) for _ in 1:nth],[Vector{ComplexF64}(undef,Mk) for _ in 1:nth])
end

########################################
######## RADIAL INTERPOLATION RANGE ####
########################################

@inline function _update_dlp_radius_bounds(rmin::Float64,rmax::Float64,r::Float64)
    r>0.0||return rmin,rmax
    return min(rmin,r),max(rmax,r)
end

function _dlp_cheb_rmin_rmax(bp::BoundaryPoints{T},::Nothing) where {T<:Real}
    N=length(bp)
    xy=bp.xy
    tol2=eps(T)^2
    rmin=Inf
    rmax=0.0
    @inbounds for j in 2:N
        point_j=xy[j]
        for i in 1:j-1
            point_i=xy[i]
            dx=point_i[1]-point_j[1]
            dy=point_i[2]-point_j[2]
            d2=muladd(dx,dx,dy*dy)
            d2<=tol2&&continue
            r=sqrt(Float64(d2))
            rmin,rmax=_update_dlp_radius_bounds(rmin,rmax,r)
        end
    end
    isfinite(rmin)&&rmax>0.0||throw(ArgumentError("Unable to determine a nonzero DLP Chebyshev radial interval"))
    return rmin,rmax
end

function _dlp_cheb_rmin_rmax(bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T}) where {T<:Real}
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    ng=orbit_size(orbits)
    xy=bp.xy
    tol2=eps(T)^2
    rmin=Inf
    rmax=0.0
    @inbounds for b in 1:m
        j=Ifund[b]
        for a in 1:m
            i=Ifund[a]
            point_i=xy[i]
            for l in 1:ng
                q=orbits.fund_to_full[l,b]
                l==1&&i==j&&continue
                point_q=xy[q]
                dx=point_i[1]-point_q[1]
                dy=point_i[2]-point_q[2]
                d2=muladd(dx,dx,dy*dy)
                d2<=tol2&&continue
                r=sqrt(Float64(d2))
                rmin,rmax=_update_dlp_radius_bounds(rmin,rmax,r)
            end
        end
    end
    isfinite(rmin)&&rmax>0.0||throw(ArgumentError("Unable to determine a nonzero reduced DLP Chebyshev radial interval"))
    return rmin,rmax
end

####################################
######## NO-SYMMETRY CHEBYSHEV #####
####################################

function _all_k_nosymm_DLP_chebyshev!(Ks::Vector{Matrix{Complex{T}}},bp::BoundaryPoints{T},plans1::Vector{ChebHankelPlanH};multithreaded::Bool=true) where {T<:Real}
    Mk=length(plans1)
    N=length(bp)
    @assert length(Ks)==Mk
    tol2=eps(T)^2
    for K in Ks
        @assert size(K)==(N,N)
        fill!(K,zero(eltype(K)))
    end
    pref=Vector{Complex{T}}(undef,Mk)
    @inbounds for m in 1:Mk
        pref[m]=Complex{T}(0,one(T)/2)*Complex{T}(plans1[m].k)
    end
    h1_tls=[Vector{ComplexF64}(undef,Mk) for _ in 1:Threads.nthreads()]
    xy=bp.xy
    nrm=bp.normal
    κ=bp.curvature
    @use_threads multithreading=multithreaded for i in 1:N
        point_i=xy[i]
        normal_i=nrm[i]
        tid=Threads.threadid()
        h1vals=h1_tls[tid]
        @inbounds begin
            val=-Complex{T}(κ[i]*_DLP_CHEB_INV_TWO_PI,zero(T))
            for m in 1:Mk
                Ks[m][i,i]=val
            end
            for j in 1:i-1
                point_j=xy[j]
                normal_j=nrm[j]
                dx=point_i[1]-point_j[1]
                dy=point_i[2]-point_j[2]
                d2=muladd(dx,dx,dy*dy)
                d2<=tol2&&throw(ArgumentError("Distinct DLP nodes coincide at indices ($i,$j)"))
                r=sqrt(d2)
                invr=inv(r)
                p,t=panel_t(plans1[1],Float64(r))
                h1_multi_ks_at_r!(h1vals,plans1,p,t,Float64(r))
                for m in 1:Mk
                    _accum_dlp_default_nosym!(Ks[m],i,j,normal_i[1],normal_i[2],normal_j[1],normal_j[2],dx,dy,invr,pref[m]*h1vals[m])
                end
            end
        end
    end
    return nothing
end

function _one_k_nosymm_DLP_chebyshev!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},plan1::ChebHankelPlanH;multithreaded::Bool=true) where {T<:Real}
    N=length(bp)
    @assert size(K)==(N,N)
    tol2=eps(T)^2
    k=Complex{T}(plan1.k)
    pref=Complex{T}(0,one(T)/2)*k
    fill!(K,zero(eltype(K)))
    xy=bp.xy
    nrm=bp.normal
    κ=bp.curvature
    @use_threads multithreading=multithreaded for i in 1:N
        point_i=xy[i]
        normal_i=nrm[i]
        @inbounds begin
            K[i,i]=-Complex{T}(κ[i]*_DLP_CHEB_INV_TWO_PI,zero(T))
            for j in 1:i-1
                point_j=xy[j]
                normal_j=nrm[j]
                dx=point_i[1]-point_j[1]
                dy=point_i[2]-point_j[2]
                d2=muladd(dx,dx,dy*dy)
                d2<=tol2&&throw(ArgumentError("Distinct DLP nodes coincide at indices ($i,$j)"))
                r=sqrt(d2)
                invr=inv(r)
                p,t=panel_t(plan1,Float64(r))
                h=pref*h1_at_r(plan1,p,t,Float64(r))
                _accum_dlp_default_nosym!(K,i,j,normal_i[1],normal_i[2],normal_j[1],normal_j[2],dx,dy,invr,h)
            end
        end
    end
    return nothing
end

function _all_k_nosymm_DLP_chebyshev_derivatives!(Ks::Vector{Matrix{Complex{T}}},dKs::Vector{Matrix{Complex{T}}},ddKs::Vector{Matrix{Complex{T}}},bp::BoundaryPoints{T},plans0::Vector{ChebHankelPlanH},plans1::Vector{ChebHankelPlanH};multithreaded::Bool=true,ws::Union{Nothing,DLPDerivChebWorkspace}=nothing) where {T<:Real}
    Mk=length(plans0)
    N=length(bp)
    @assert length(plans1)==Mk
    @assert length(Ks)==Mk
    @assert length(dKs)==Mk
    @assert length(ddKs)==Mk
    for m in 1:Mk
        @assert size(Ks[m])==(N,N)
        @assert size(dKs[m])==(N,N)
        @assert size(ddKs[m])==(N,N)
        fill!(Ks[m],zero(eltype(Ks[m])))
        fill!(dKs[m],zero(eltype(dKs[m])))
        fill!(ddKs[m],zero(eltype(ddKs[m])))
    end
    kvec=ComplexF64[plans0[m].k for m in 1:Mk]
    local_ws=isnothing(ws) ? DLPDerivChebWorkspace(Mk) : ws
    length(local_ws.h0_tls)>=Threads.nthreads()||throw(ArgumentError("Derivative DLP workspace has $(length(local_ws.h0_tls)) thread-local buffers but Threads.nthreads()=$(Threads.nthreads())"))
    length(local_ws.h1_tls)>=Threads.nthreads()||throw(ArgumentError("Derivative DLP workspace has $(length(local_ws.h1_tls)) thread-local buffers but Threads.nthreads()=$(Threads.nthreads())"))
    xy=bp.xy
    nrm=bp.normal
    κ=bp.curvature
    tol2=eps(T)^2
    @use_threads multithreading=multithreaded for i in 1:N
        point_i=xy[i]
        normal_i=nrm[i]
        tid=Threads.threadid()
        h0vals=local_ws.h0_tls[tid]
        h1vals=local_ws.h1_tls[tid]
        @inbounds begin
            val=-Complex{T}(κ[i]*_DLP_CHEB_INV_TWO_PI,zero(T))
            for m in 1:Mk
                Ks[m][i,i]=val
            end
            for j in 1:i-1
                point_j=xy[j]
                normal_j=nrm[j]
                dx=point_i[1]-point_j[1]
                dy=point_i[2]-point_j[2]
                d2=muladd(dx,dx,dy*dy)
                d2<=tol2&&throw(ArgumentError("Distinct DLP nodes coincide at indices ($i,$j)"))
                r=sqrt(d2)
                invr=inv(r)
                p,t=panel_t(plans0[1],Float64(r))
                h0_h1_multi_ks_at_r!(h0vals,h1vals,plans0,plans1,p,t,Float64(r))
                for m in 1:Mk
                    hK,hdK,hddK=_kernel_triplet_from_hankels(kvec[m],r,h0vals[m],h1vals[m])
                    _accum_dlp_triplet_nosym!(Ks[m],dKs[m],ddKs[m],i,j,normal_i[1],normal_i[2],normal_j[1],normal_j[2],dx,dy,invr,hK,hdK,hddK)
                end
            end
        end
    end
    return nothing
end

function _one_k_nosymm_DLP_chebyshev_derivatives!(K::AbstractMatrix{Complex{T}},dK::AbstractMatrix{Complex{T}},ddK::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},plan0::ChebHankelPlanH,plan1::ChebHankelPlanH;multithreaded::Bool=true) where {T<:Real}
    N=length(bp)
    @assert size(K)==(N,N)
    @assert size(dK)==(N,N)
    @assert size(ddK)==(N,N)
    k=ComplexF64(plan0.k)
    fill!(K,zero(eltype(K)))
    fill!(dK,zero(eltype(dK)))
    fill!(ddK,zero(eltype(ddK)))
    xy=bp.xy
    nrm=bp.normal
    κ=bp.curvature
    tol2=eps(T)^2
    plans0=[plan0]
    plans1=[plan1]
    h0_tls=[Vector{ComplexF64}(undef,1) for _ in 1:Threads.nthreads()]
    h1_tls=[Vector{ComplexF64}(undef,1) for _ in 1:Threads.nthreads()]
    @use_threads multithreading=multithreaded for i in 1:N
        tid=Threads.threadid()
        h0vals=h0_tls[tid]
        h1vals=h1_tls[tid]
        point_i=xy[i]
        normal_i=nrm[i]
        @inbounds begin
            K[i,i]=-Complex{T}(κ[i]*_DLP_CHEB_INV_TWO_PI,zero(T))
            for j in 1:i-1
                point_j=xy[j]
                normal_j=nrm[j]
                dx=point_i[1]-point_j[1]
                dy=point_i[2]-point_j[2]
                d2=muladd(dx,dx,dy*dy)
                d2<=tol2&&throw(ArgumentError("Distinct DLP nodes coincide at indices ($i,$j)"))
                r=sqrt(d2)
                invr=inv(r)
                p,t=panel_t(plan0,Float64(r))
                H0,H1=_h0_h1_one_k_at_r!(h0vals,h1vals,plans0,plans1,p,t,Float64(r))
                hK,hdK,hddK=_kernel_triplet_from_hankels(k,r,H0,H1)
                _accum_dlp_triplet_nosym!(K,dK,ddK,i,j,normal_i[1],normal_i[2],normal_j[1],normal_j[2],dx,dy,invr,hK,hdK,hddK)
            end
        end
    end
    return nothing
end

####################################
######## REDUCED CHEBYSHEV #########
####################################

function _all_k_reduced_DLP_chebyshev!(Ks::Vector{Matrix{Complex{T}}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},plans1::Vector{ChebHankelPlanH};multithreaded::Bool=true) where {T<:Real}
    Mk=length(plans1)
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    ng=orbit_size(orbits)
    @assert length(Ks)==Mk
    for K in Ks
        @assert size(K)==(m,m)
        fill!(K,zero(eltype(K)))
    end
    pref=Vector{Complex{T}}(undef,Mk)
    @inbounds for q in 1:Mk
        pref[q]=Complex{T}(0,one(T)/2)*Complex{T}(plans1[q].k)
    end
    xy=bp.xy
    nrm=bp.normal
    κ=bp.curvature
    ds=bp.ds
    tol2=eps(T)^2
    h1_tls=[Vector{ComplexF64}(undef,Mk) for _ in 1:Threads.nthreads()]
    acc_tls=[Vector{Complex{T}}(undef,Mk) for _ in 1:Threads.nthreads()]
    @use_threads multithreading=multithreaded for b in 1:m
        tid=Threads.threadid()
        h1vals=h1_tls[tid]
        acc=acc_tls[tid]
        j=Ifund[b]
        wj=ds[j]
        @inbounds for a in 1:m
            fill!(acc,zero(Complex{T}))
            i=Ifund[a]
            point_i=xy[i]
            for l in 1:ng
                qimg=orbits.fund_to_full[l,b]
                scale=orbits.fund_to_scale[l,b]
                weight_ratio=ds[qimg]/wj
                if l==1&&i==j
                    d0=scale*(-Complex{T}(κ[i]*_DLP_CHEB_INV_TWO_PI,zero(T)))*weight_ratio
                    for q in 1:Mk
                        acc[q]+=d0
                    end
                    continue
                end
                point_q=xy[qimg]
                dx=point_i[1]-point_q[1]
                dy=point_i[2]-point_q[2]
                d2=muladd(dx,dx,dy*dy)
                d2<=tol2&&throw(ArgumentError("A nonidentity DLP symmetry image coincides with reduced target index $i"))
                r=sqrt(d2)
                invr=inv(r)
                normal_q=nrm[qimg]
                c=(normal_q[1]*dx+normal_q[2]*dy)*invr
                p,t=panel_t(plans1[1],Float64(r))
                h1_multi_ks_at_r!(h1vals,plans1,p,t,Float64(r))
                s=scale*weight_ratio*c
                for q in 1:Mk
                    acc[q]+=s*pref[q]*h1vals[q]
                end
            end
            for q in 1:Mk
                Ks[q][a,b]=acc[q]
            end
        end
    end
    return nothing
end

function _one_k_reduced_DLP_chebyshev!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},plan1::ChebHankelPlanH;multithreaded::Bool=true) where {T<:Real}
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    ng=orbit_size(orbits)
    @assert size(K)==(m,m)
    fill!(K,zero(eltype(K)))
    k=Complex{T}(plan1.k)
    pref=Complex{T}(0,one(T)/2)*k
    xy=bp.xy
    nrm=bp.normal
    κ=bp.curvature
    ds=bp.ds
    tol2=eps(T)^2
    @use_threads multithreading=multithreaded for b in 1:m
        j=Ifund[b]
        wj=ds[j]
        @inbounds for a in 1:m
            i=Ifund[a]
            point_i=xy[i]
            val=zero(Complex{T})
            for l in 1:ng
                qimg=orbits.fund_to_full[l,b]
                scale=orbits.fund_to_scale[l,b]
                weight_ratio=ds[qimg]/wj
                if l==1&&i==j
                    val+=scale*(-Complex{T}(κ[i]*_DLP_CHEB_INV_TWO_PI,zero(T)))*weight_ratio
                    continue
                end
                point_q=xy[qimg]
                dx=point_i[1]-point_q[1]
                dy=point_i[2]-point_q[2]
                d2=muladd(dx,dx,dy*dy)
                d2<=tol2&&throw(ArgumentError("A nonidentity DLP symmetry image coincides with reduced target index $i"))
                r=sqrt(d2)
                normal_q=nrm[qimg]
                c=(normal_q[1]*dx+normal_q[2]*dy)/r
                p,t=panel_t(plan1,Float64(r))
                h1=h1_at_r(plan1,p,t,Float64(r))
                val+=scale*weight_ratio*c*pref*h1
            end
            K[a,b]=val
        end
    end
    return nothing
end

function _all_k_reduced_DLP_chebyshev_derivatives!(Ks::Vector{Matrix{Complex{T}}},dKs::Vector{Matrix{Complex{T}}},ddKs::Vector{Matrix{Complex{T}}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},plans0::Vector{ChebHankelPlanH},plans1::Vector{ChebHankelPlanH};multithreaded::Bool=true,ws::Union{Nothing,DLPDerivChebWorkspace}=nothing) where {T<:Real}
    Mk=length(plans0)
    @assert length(plans1)==Mk
    @assert length(Ks)==Mk
    @assert length(dKs)==Mk
    @assert length(ddKs)==Mk
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    ng=orbit_size(orbits)
    for q in 1:Mk
        @assert size(Ks[q])==(m,m)
        @assert size(dKs[q])==(m,m)
        @assert size(ddKs[q])==(m,m)
        fill!(Ks[q],zero(eltype(Ks[q])))
        fill!(dKs[q],zero(eltype(dKs[q])))
        fill!(ddKs[q],zero(eltype(ddKs[q])))
    end
    kvec=ComplexF64[plans0[q].k for q in 1:Mk]
    local_ws=isnothing(ws) ? DLPDerivChebWorkspace(Mk) : ws
    length(local_ws.h0_tls)>=Threads.nthreads()||throw(ArgumentError("Derivative DLP workspace has $(length(local_ws.h0_tls)) thread-local buffers but Threads.nthreads()=$(Threads.nthreads())"))
    length(local_ws.h1_tls)>=Threads.nthreads()||throw(ArgumentError("Derivative DLP workspace has $(length(local_ws.h1_tls)) thread-local buffers but Threads.nthreads()=$(Threads.nthreads())"))
    ntls=length(local_ws.h0_tls)
    acc_tls=[Vector{Complex{T}}(undef,Mk) for _ in 1:ntls]
    acc1_tls=[Vector{Complex{T}}(undef,Mk) for _ in 1:ntls]
    acc2_tls=[Vector{Complex{T}}(undef,Mk) for _ in 1:ntls]
    xy=bp.xy
    nrm=bp.normal
    κ=bp.curvature
    ds=bp.ds
    tol2=eps(T)^2
    @use_threads multithreading=multithreaded for b in 1:m
        tid=Threads.threadid()
        h0vals=local_ws.h0_tls[tid]
        h1vals=local_ws.h1_tls[tid]
        acc=acc_tls[tid]
        acc1=acc1_tls[tid]
        acc2=acc2_tls[tid]
        j=Ifund[b]
        wj=ds[j]
        @inbounds for a in 1:m
            fill!(acc,zero(Complex{T}))
            fill!(acc1,zero(Complex{T}))
            fill!(acc2,zero(Complex{T}))
            i=Ifund[a]
            point_i=xy[i]
            for l in 1:ng
                qimg=orbits.fund_to_full[l,b]
                scale=orbits.fund_to_scale[l,b]
                weight_ratio=ds[qimg]/wj
                if l==1&&i==j
                    d0=scale*(-Complex{T}(κ[i]*_DLP_CHEB_INV_TWO_PI,zero(T)))*weight_ratio
                    for q in 1:Mk
                        acc[q]+=d0
                    end
                    continue
                end
                point_q=xy[qimg]
                dx=point_i[1]-point_q[1]
                dy=point_i[2]-point_q[2]
                d2=muladd(dx,dx,dy*dy)
                d2<=tol2&&throw(ArgumentError("A nonidentity DLP symmetry image coincides with reduced target index $i"))
                r=sqrt(d2)
                invr=inv(r)
                normal_q=nrm[qimg]
                c=(normal_q[1]*dx+normal_q[2]*dy)*invr
                p,t=panel_t(plans0[1],Float64(r))
                h0_h1_multi_ks_at_r!(h0vals,h1vals,plans0,plans1,p,t,Float64(r))
                s=scale*weight_ratio*c
                for q in 1:Mk
                    hK,hdK,hddK=_kernel_triplet_from_hankels(kvec[q],r,h0vals[q],h1vals[q])
                    acc[q]+=s*hK
                    acc1[q]+=s*hdK
                    acc2[q]+=s*hddK
                end
            end
            for q in 1:Mk
                Ks[q][a,b]=acc[q]
                dKs[q][a,b]=acc1[q]
                ddKs[q][a,b]=acc2[q]
            end
        end
    end
    return nothing
end

function _one_k_reduced_DLP_chebyshev_derivatives!(K::AbstractMatrix{Complex{T}},dK::AbstractMatrix{Complex{T}},ddK::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},plan0::ChebHankelPlanH,plan1::ChebHankelPlanH;multithreaded::Bool=true) where {T<:Real}
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    ng=orbit_size(orbits)
    @assert size(K)==(m,m)
    @assert size(dK)==(m,m)
    @assert size(ddK)==(m,m)
    fill!(K,zero(eltype(K)))
    fill!(dK,zero(eltype(dK)))
    fill!(ddK,zero(eltype(ddK)))
    k=ComplexF64(plan0.k)
    xy=bp.xy
    nrm=bp.normal
    κ=bp.curvature
    ds=bp.ds
    tol2=eps(T)^2
    plans0=[plan0]
    plans1=[plan1]
    h0_tls=[Vector{ComplexF64}(undef,1) for _ in 1:Threads.nthreads()]
    h1_tls=[Vector{ComplexF64}(undef,1) for _ in 1:Threads.nthreads()]
    @use_threads multithreading=multithreaded for b in 1:m
        tid=Threads.threadid()
        h0vals=h0_tls[tid]
        h1vals=h1_tls[tid]
        j=Ifund[b]
        wj=ds[j]
        @inbounds for a in 1:m
            i=Ifund[a]
            point_i=xy[i]
            val=zero(Complex{T})
            val1=zero(Complex{T})
            val2=zero(Complex{T})
            for l in 1:ng
                qimg=orbits.fund_to_full[l,b]
                scale=orbits.fund_to_scale[l,b]
                weight_ratio=ds[qimg]/wj
                if l==1&&i==j
                    val+=scale*(-Complex{T}(κ[i]*_DLP_CHEB_INV_TWO_PI,zero(T)))*weight_ratio
                    continue
                end
                point_q=xy[qimg]
                dx=point_i[1]-point_q[1]
                dy=point_i[2]-point_q[2]
                d2=muladd(dx,dx,dy*dy)
                d2<=tol2&&throw(ArgumentError("A nonidentity DLP symmetry image coincides with reduced target index $i"))
                r=sqrt(d2)
                normal_q=nrm[qimg]
                c=(normal_q[1]*dx+normal_q[2]*dy)/r
                p,t=panel_t(plan0,Float64(r))
                H0,H1=_h0_h1_one_k_at_r!(h0vals,h1vals,plans0,plans1,p,t,Float64(r))
                hK,hdK,hddK=_kernel_triplet_from_hankels(k,r,H0,H1)
                s=scale*weight_ratio*c
                val+=s*hK
                val1+=s*hdK
                val2+=s*hddK
            end
            K[a,b]=val
            dK[a,b]=val1
            ddK[a,b]=val2
        end
    end
    return nothing
end

####################################
######## DLP PUBLIC DISPATCH #######
####################################

"""
    compute_kernel_matrices_DLP_chebyshev!(Ks::Vector{Matrix{Complex{T}}},bp::BoundaryPoints{T},plans::Vector{ChebHankelPlanH};multithreaded::Bool=true) where {T<:Real} → Nothing

Assemble full-boundary raw DLP kernels for several wavenumbers.

## Arguments
* `Ks::Vector{Matrix{Complex{T}}}`: Destination raw kernel matrices.
* `bp::BoundaryPoints{T}`: Full boundary discretization.
* `plans::Vector{ChebHankelPlanH}`: `H₁⁽¹⁾` Chebyshev plans.

## Keyword Arguments
* `multithreaded::Bool`: Enable threaded assembly.

## Returns
* `nothing`.
"""
compute_kernel_matrices_DLP_chebyshev!(Ks::Vector{Matrix{Complex{T}}},bp::BoundaryPoints{T},plans::Vector{ChebHankelPlanH};multithreaded::Bool=true) where {T<:Real}=_all_k_nosymm_DLP_chebyshev!(Ks,bp,plans;multithreaded=multithreaded)

"""
    compute_kernel_matrices_DLP_chebyshev!(Ks::Vector{Matrix{Complex{T}}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},plans::Vector{ChebHankelPlanH};multithreaded::Bool=true) where {T<:Real} → Nothing

Assemble exact symmetry-reduced raw DLP kernels for several wavenumbers.

## Arguments
* `Ks::Vector{Matrix{Complex{T}}}`: Destination reduced kernel matrices.
* `bp::BoundaryPoints{T}`: Full boundary discretization.
* `orbits::SymmetryOrbitMap{T}`: Exact symmetry-orbit map.
* `plans::Vector{ChebHankelPlanH}`: `H₁⁽¹⁾` Chebyshev plans.

## Keyword Arguments
* `multithreaded::Bool`: Enable threaded assembly.

## Returns
* `nothing`.
"""
compute_kernel_matrices_DLP_chebyshev!(Ks::Vector{Matrix{Complex{T}}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},plans::Vector{ChebHankelPlanH};multithreaded::Bool=true) where {T<:Real}=_all_k_reduced_DLP_chebyshev!(Ks,bp,orbits,plans;multithreaded=multithreaded)

compute_kernel_matrices_DLP_chebyshev!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},plan::ChebHankelPlanH;multithreaded::Bool=true) where {T<:Real}=_one_k_nosymm_DLP_chebyshev!(K,bp,plan;multithreaded=multithreaded)

compute_kernel_matrices_DLP_chebyshev!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},plan::ChebHankelPlanH;multithreaded::Bool=true) where {T<:Real}=_one_k_reduced_DLP_chebyshev!(K,bp,orbits,plan;multithreaded=multithreaded)

function compute_kernel_matrices_DLP_chebyshev!(Ks::Vector{Matrix{Complex{T}}},bp::BoundaryPoints{T},symmetry,plans::Vector{ChebHankelPlanH};multithreaded::Bool=true) where {T<:Real}
    isnothing(symmetry)&&return _all_k_nosymm_DLP_chebyshev!(Ks,bp,plans;multithreaded=multithreaded)
    orbits=symmetry_index_orbits(T,bp,symmetry)
    return _all_k_reduced_DLP_chebyshev!(Ks,bp,orbits,plans;multithreaded=multithreaded)
end

function compute_kernel_matrices_DLP_chebyshev!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},symmetry,plan::ChebHankelPlanH;multithreaded::Bool=true) where {T<:Real}
    isnothing(symmetry)&&return _one_k_nosymm_DLP_chebyshev!(K,bp,plan;multithreaded=multithreaded)
    orbits=symmetry_index_orbits(T,bp,symmetry)
    return _one_k_reduced_DLP_chebyshev!(K,bp,orbits,plan;multithreaded=multithreaded)
end

compute_kernel_matrices_DLP_chebyshev_derivatives!(K::AbstractMatrix{Complex{T}},dK::AbstractMatrix{Complex{T}},ddK::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},plan0::ChebHankelPlanH,plan1::ChebHankelPlanH;multithreaded::Bool=true) where {T<:Real}=_one_k_nosymm_DLP_chebyshev_derivatives!(K,dK,ddK,bp,plan0,plan1;multithreaded=multithreaded)

compute_kernel_matrices_DLP_chebyshev_derivatives!(K::AbstractMatrix{Complex{T}},dK::AbstractMatrix{Complex{T}},ddK::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},plan0::ChebHankelPlanH,plan1::ChebHankelPlanH;multithreaded::Bool=true) where {T<:Real}=_one_k_reduced_DLP_chebyshev_derivatives!(K,dK,ddK,bp,orbits,plan0,plan1;multithreaded=multithreaded)

function compute_kernel_matrices_DLP_chebyshev_derivatives!(Ks::Vector{Matrix{Complex{T}}},dKs::Vector{Matrix{Complex{T}}},ddKs::Vector{Matrix{Complex{T}}},bp::BoundaryPoints{T},symmetry,plans0::Vector{ChebHankelPlanH},plans1::Vector{ChebHankelPlanH};multithreaded::Bool=true,ws::Union{Nothing,DLPDerivChebWorkspace}=nothing) where {T<:Real}
    if isnothing(symmetry)
        return _all_k_nosymm_DLP_chebyshev_derivatives!(Ks,dKs,ddKs,bp,plans0,plans1;multithreaded=multithreaded,ws=ws)
    end
    orbits=symmetry_index_orbits(T,bp,symmetry)
    return _all_k_reduced_DLP_chebyshev_derivatives!(Ks,dKs,ddKs,bp,orbits,plans0,plans1;multithreaded=multithreaded,ws=ws)
end

function compute_kernel_matrices_DLP_chebyshev_derivatives!(Ks::Vector{Matrix{Complex{T}}},dKs::Vector{Matrix{Complex{T}}},ddKs::Vector{Matrix{Complex{T}}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},plans0::Vector{ChebHankelPlanH},plans1::Vector{ChebHankelPlanH};multithreaded::Bool=true,ws::Union{Nothing,DLPDerivChebWorkspace}=nothing) where {T<:Real}
    return _all_k_reduced_DLP_chebyshev_derivatives!(Ks,dKs,ddKs,bp,orbits,plans0,plans1;multithreaded=multithreaded,ws=ws)
end

####################################
######## FREDHOLM ASSEMBLY #########
####################################

"""
    assemble_fredholm_matrices!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T}) where {T<:Real} → Nothing

Apply full-boundary Nyström source weights and form

    A=I-KW.

## Arguments
* `K::AbstractMatrix{Complex{T}}`: Raw full-boundary DLP matrix.
* `bp::BoundaryPoints{T}`: Full boundary discretization.

## Returns
* `nothing`.
"""
function assemble_fredholm_matrices!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T}) where {T<:Real}
    @inbounds for j in eachindex(bp.ds)
        @views K[:,j].*=-bp.ds[j]
        K[j,j]+=one(eltype(K))
    end
    filter_matrix!(K)
    return nothing
end

"""
    assemble_fredholm_matrices!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T}) where {T<:Real} → Nothing

Apply representative Nyström source weights to a symmetry-reduced raw DLP
kernel and form

    A_red=I-K_red W_red.

The raw reduced kernel already contains the image-weight ratio
`ds[q]/ds[j]`.

## Arguments
* `K::AbstractMatrix{Complex{T}}`: Raw reduced DLP matrix.
* `bp::BoundaryPoints{T}`: Full boundary discretization.
* `orbits::SymmetryOrbitMap{T}`: Exact orbit map.

## Returns
* `nothing`.
"""
function assemble_fredholm_matrices!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T}) where {T<:Real}
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    @assert size(K)==(m,m)
    @inbounds for b in 1:m
        @views K[:,b].*=-bp.ds[Ifund[b]]
    end
    @inbounds for a in 1:m
        K[a,a]+=one(eltype(K))
    end
    filter_matrix!(K)
    return nothing
end

function assemble_fredholm_matrices!(Ks::Vector{Matrix{Complex{T}}},bp::BoundaryPoints{T}) where {T<:Real}
    Threads.@threads for q in eachindex(Ks)
        assemble_fredholm_matrices!(Ks[q],bp)
    end
    return nothing
end

function assemble_fredholm_matrices!(Ks::Vector{Matrix{Complex{T}}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T}) where {T<:Real}
    Threads.@threads for q in eachindex(Ks)
        assemble_fredholm_matrices!(Ks[q],bp,orbits)
    end
    return nothing
end

"""
    assemble_fredholm_matrices_with_derivatives!(K::AbstractMatrix{Complex{T}},dK::AbstractMatrix{Complex{T}},ddK::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T}) where {T<:Real} → Nothing

Apply full-boundary Nyström weights and form `A`, `A'`, and `A''`.

## Arguments
* `K::AbstractMatrix{Complex{T}}`: Raw DLP matrix.
* `dK::AbstractMatrix{Complex{T}}`: Raw first derivative.
* `ddK::AbstractMatrix{Complex{T}}`: Raw second derivative.
* `bp::BoundaryPoints{T}`: Full boundary discretization.

## Returns
* `nothing`.
"""
function assemble_fredholm_matrices_with_derivatives!(K::AbstractMatrix{Complex{T}},dK::AbstractMatrix{Complex{T}},ddK::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T}) where {T<:Real}
    @inbounds for j in eachindex(bp.ds)
        s=-bp.ds[j]
        @views begin
            K[:,j].*=s
            dK[:,j].*=s
            ddK[:,j].*=s
        end
        K[j,j]+=one(eltype(K))
    end
    filter_matrix!(K)
    filter_matrix!(dK)
    filter_matrix!(ddK)
    return nothing
end

"""
    assemble_fredholm_matrices_with_derivatives!(K::AbstractMatrix{Complex{T}},dK::AbstractMatrix{Complex{T}},ddK::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T}) where {T<:Real} → Nothing

Apply the representative reduced Nyström weights and form the reduced
Fredholm matrix and its first two derivatives.

## Arguments
* `K::AbstractMatrix{Complex{T}}`: Raw reduced DLP matrix.
* `dK::AbstractMatrix{Complex{T}}`: Raw reduced first derivative.
* `ddK::AbstractMatrix{Complex{T}}`: Raw reduced second derivative.
* `bp::BoundaryPoints{T}`: Full boundary discretization.
* `orbits::SymmetryOrbitMap{T}`: Exact symmetry-orbit map.

## Returns
* `nothing`.
"""
function assemble_fredholm_matrices_with_derivatives!(K::AbstractMatrix{Complex{T}},dK::AbstractMatrix{Complex{T}},ddK::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T}) where {T<:Real}
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    @assert size(K)==(m,m)
    @assert size(dK)==(m,m)
    @assert size(ddK)==(m,m)
    @inbounds for b in 1:m
        s=-bp.ds[Ifund[b]]
        @views begin
            K[:,b].*=s
            dK[:,b].*=s
            ddK[:,b].*=s
        end
    end
    @inbounds for a in 1:m
        K[a,a]+=one(eltype(K))
    end
    filter_matrix!(K)
    filter_matrix!(dK)
    filter_matrix!(ddK)
    return nothing
end

function assemble_fredholm_matrices_with_derivatives!(Ks::Vector{Matrix{Complex{T}}},dKs::Vector{Matrix{Complex{T}}},ddKs::Vector{Matrix{Complex{T}}},bp::BoundaryPoints{T}) where {T<:Real}
    @assert length(Ks)==length(dKs)==length(ddKs)
    Threads.@threads for q in eachindex(Ks)
        assemble_fredholm_matrices_with_derivatives!(Ks[q],dKs[q],ddKs[q],bp)
    end
    return nothing
end

function assemble_fredholm_matrices_with_derivatives!(Ks::Vector{Matrix{Complex{T}}},dKs::Vector{Matrix{Complex{T}}},ddKs::Vector{Matrix{Complex{T}}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T}) where {T<:Real}
    @assert length(Ks)==length(dKs)==length(ddKs)
    Threads.@threads for q in eachindex(Ks)
        assemble_fredholm_matrices_with_derivatives!(Ks[q],dKs[q],ddKs[q],bp,orbits)
    end
    return nothing
end

################################################################################
################ COMMON DERIVATIVE CHEBYSHEV WORKSPACE API #####################
################################################################################

"""
    DLPDerivativeChebyshevWorkspace{O}

Reusable derivative-aware Chebyshev workspace for the plain DLP backend.

The workspace owns all solver-specific state needed to construct `A(k)`,
`A'(k)`, and `A''(k)` repeatedly without exposing DLP implementation details to
higher-level algorithms such as EBIM.

## Attributes
* `plans0::Vector{ChebHankelPlanH}`: `H₀⁽¹⁾` Chebyshev plans.
* `plans1::Vector{ChebHankelPlanH}`: `H₁⁽¹⁾` Chebyshev plans.
* `bessel_ws::DLPDerivChebWorkspace`: Thread-local Hankel buffers.
* `orbits::O`: Exact symmetry-orbit map, or `nothing` without symmetry.
"""
struct DLPDerivativeChebyshevWorkspace{O}
    plans0::Vector{ChebHankelPlanH}
    plans1::Vector{ChebHankelPlanH}
    bessel_ws::DLPDerivChebWorkspace
    orbits::O
end

"""
    build_derivative_chebyshev_workspace(solver::BoundaryIntegralMethod,pts::BoundaryPoints{T},ks::AbstractVector{<:Number};n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false) where {T<:Real} → DLPDerivativeChebyshevWorkspace

Build a reusable derivative-aware Chebyshev workspace for the plain DLP
backend.

The radial interpolation range is determined from the same full or reduced
geometry used during matrix construction, and the exact `SymmetryOrbitMap` is
cached once for reuse by all matrix evaluations.

## Arguments
* `solver::BoundaryIntegralMethod`: Plain DLP solver.
* `pts::BoundaryPoints{T}`: Full boundary discretization.
* `ks::AbstractVector{<:Number}`: Wavenumbers represented by the workspace.

## Keyword Arguments
* `n_panels_h::Int`: Hankel Chebyshev panel count.
* `M_h::Int`: Hankel Chebyshev polynomial degree.
* `n_panels_j::Int`: Compatibility keyword; unused by the plain DLP backend.
* `M_j::Int`: Compatibility keyword; unused by the plain DLP backend.
* `timeit::Bool`: Enable timing diagnostics.

## Returns
* `ws::DLPDerivativeChebyshevWorkspace`: Reusable derivative workspace.
"""
function build_derivative_chebyshev_workspace(solver::BoundaryIntegralMethod,pts::BoundaryPoints{T},ks::AbstractVector{<:Number};n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false) where {T<:Real}
    zks=ComplexF64.(ks)
    orbits=_dlp_symmetry_orbits(solver,pts)
    rmin,rmax=_dlp_cheb_rmin_rmax(pts,orbits)
    plans0=Vector{ChebHankelPlanH}(undef,length(zks))
    plans1=Vector{ChebHankelPlanH}(undef,length(zks))
    @benchit timeit=timeit "DLP derivative plans" Threads.@threads for q in eachindex(zks)
        plans0[q]=plan_h(0,1,zks[q],rmin,rmax;npanels=n_panels_h,M=M_h)
        plans1[q]=plan_h(1,1,zks[q],rmin,rmax;npanels=n_panels_h,M=M_h)
    end
    bessel_ws=DLPDerivChebWorkspace(length(zks))
    return DLPDerivativeChebyshevWorkspace(plans0,plans1,bessel_ws,orbits)
end

"""
    construct_matrices_chebyshev_with_derivatives!(As::Vector{Matrix{Complex{T}}},dAs::Vector{Matrix{Complex{T}}},ddAs::Vector{Matrix{Complex{T}}},solver::BoundaryIntegralMethod,pts::BoundaryPoints{T},ws::DLPDerivativeChebyshevWorkspace;multithreaded::Bool=true) where {T<:Real} → Nothing

Construct all plain-DLP Fredholm matrices and their first two wavenumber
derivatives from a reusable derivative Chebyshev workspace.

All DLP-specific symmetry and quadrature logic remains inside this backend.

## Arguments
* `As::Vector{Matrix{Complex{T}}}`: Destination Fredholm matrices.
* `dAs::Vector{Matrix{Complex{T}}}`: Destination first derivatives.
* `ddAs::Vector{Matrix{Complex{T}}}`: Destination second derivatives.
* `solver::BoundaryIntegralMethod`: Plain DLP solver.
* `pts::BoundaryPoints{T}`: Full boundary discretization.
* `ws::DLPDerivativeChebyshevWorkspace`: Reusable derivative workspace.

## Keyword Arguments
* `multithreaded::Bool`: Enable threaded matrix assembly.

## Returns
* `nothing`.
"""
function construct_matrices_chebyshev_with_derivatives!(As::Vector{Matrix{Complex{T}}},dAs::Vector{Matrix{Complex{T}}},ddAs::Vector{Matrix{Complex{T}}},solver::BoundaryIntegralMethod,pts::BoundaryPoints{T},ws::DLPDerivativeChebyshevWorkspace;multithreaded::Bool=true) where {T<:Real}
    Mk=length(ws.plans0)
    @assert length(ws.plans1)==Mk
    @assert length(As)==Mk
    @assert length(dAs)==Mk
    @assert length(ddAs)==Mk
    n=boundary_matrix_size(solver,pts)
    @inbounds for q in 1:Mk
        @assert size(As[q])==(n,n) "As[$q] has size $(size(As[q])), expected ($n,$n)"
        @assert size(dAs[q])==(n,n) "dAs[$q] has size $(size(dAs[q])), expected ($n,$n)"
        @assert size(ddAs[q])==(n,n) "ddAs[$q] has size $(size(ddAs[q])), expected ($n,$n)"
    end
    if isnothing(ws.orbits)
        compute_kernel_matrices_DLP_chebyshev_derivatives!(As,dAs,ddAs,pts,nothing,ws.plans0,ws.plans1;multithreaded=multithreaded,ws=ws.bessel_ws)
        assemble_fredholm_matrices_with_derivatives!(As,dAs,ddAs,pts)
    else
        compute_kernel_matrices_DLP_chebyshev_derivatives!(As,dAs,ddAs,pts,ws.orbits,ws.plans0,ws.plans1;multithreaded=multithreaded,ws=ws.bessel_ws)
        assemble_fredholm_matrices_with_derivatives!(As,dAs,ddAs,pts,ws.orbits)
    end
    return nothing
end

"""
    construct_matrix_chebyshev_with_derivatives_at!(A::AbstractMatrix{Complex{T}},dA::AbstractMatrix{Complex{T}},ddA::AbstractMatrix{Complex{T}},solver::BoundaryIntegralMethod,pts::BoundaryPoints{T},ws::DLPDerivativeChebyshevWorkspace,idx::Int;multithreaded::Bool=true) where {T<:Real} → Nothing

Construct one plain-DLP Fredholm matrix and its first two wavenumber derivatives
from entry `idx` of a reusable derivative Chebyshev workspace.

The single-wavenumber path uses the cached plans directly and does not allocate
a temporary one-element derivative workspace.

## Arguments
* `A::AbstractMatrix{Complex{T}}`: Destination Fredholm matrix.
* `dA::AbstractMatrix{Complex{T}}`: Destination first derivative.
* `ddA::AbstractMatrix{Complex{T}}`: Destination second derivative.
* `solver::BoundaryIntegralMethod`: Plain DLP solver.
* `pts::BoundaryPoints{T}`: Full boundary discretization.
* `ws::DLPDerivativeChebyshevWorkspace`: Reusable derivative workspace.
* `idx::Int`: Cached wavenumber index.

## Keyword Arguments
* `multithreaded::Bool`: Enable threaded matrix assembly.

## Returns
* `nothing`.
"""
function construct_matrix_chebyshev_with_derivatives_at!(A::AbstractMatrix{Complex{T}},dA::AbstractMatrix{Complex{T}},ddA::AbstractMatrix{Complex{T}},solver::BoundaryIntegralMethod,pts::BoundaryPoints{T},ws::DLPDerivativeChebyshevWorkspace,idx::Int;multithreaded::Bool=true) where {T<:Real}
    checkbounds(ws.plans0,idx)
    checkbounds(ws.plans1,idx)
    n=boundary_matrix_size(solver,pts)
    @assert size(A)==(n,n)
    @assert size(dA)==(n,n)
    @assert size(ddA)==(n,n)
    if isnothing(ws.orbits)
        compute_kernel_matrices_DLP_chebyshev_derivatives!(A,dA,ddA,pts,ws.plans0[idx],ws.plans1[idx];multithreaded=multithreaded)
        assemble_fredholm_matrices_with_derivatives!(A,dA,ddA,pts)
    else
        compute_kernel_matrices_DLP_chebyshev_derivatives!(A,dA,ddA,pts,ws.orbits,ws.plans0[idx],ws.plans1[idx];multithreaded=multithreaded)
        assemble_fredholm_matrices_with_derivatives!(A,dA,ddA,pts,ws.orbits)
    end
    return nothing
end

####################################
###### CHEBYSHEV CONSTRUCTION ######
####################################

"""
    construct_matrices_chebyshev!(Tbufs::Vector{Matrix{Complex{T}}},::Val{:dlp},solver::BoundaryIntegralMethod,pts::BoundaryPoints{T},zj::AbstractVector{Complex{T}};multithreaded::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false) where {T<:Real} → Nothing

Construct the full or exact symmetry-reduced DLP Fredholm matrices

    A(k)=I-K(k)W

for all complex wavenumbers in `zj` using piecewise-Chebyshev interpolation of
`H₁⁽¹⁾`.

## Arguments
* `Tbufs::Vector{Matrix{Complex{T}}}`: Preallocated Fredholm matrices.
* `::Val{:dlp}`: Selects the direct DLP Chebyshev backend.
* `solver::BoundaryIntegralMethod`: DLP solver.
* `pts::BoundaryPoints{T}`: Full boundary discretization.
* `zj::AbstractVector{Complex{T}}`: Complex contour wavenumbers.

## Keyword Arguments
* `multithreaded::Bool`: Enable threaded matrix assembly.
* `n_panels_h::Int`: Hankel Chebyshev panel count.
* `M_h::Int`: Hankel Chebyshev degree.
* `n_panels_j::Int`: Compatibility keyword; unused by this backend.
* `M_j::Int`: Compatibility keyword; unused by this backend.
* `timeit::Bool`: Enable timing diagnostics.

## Returns
* `nothing`.
"""
function construct_matrices_chebyshev!(Tbufs::Vector{Matrix{Complex{T}}},::Val{:dlp},solver::BoundaryIntegralMethod,pts::BoundaryPoints{T},zj::AbstractVector{Complex{T}};multithreaded::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false) where {T<:Real}
    @assert length(Tbufs)==length(zj)
    orbits=_dlp_symmetry_orbits(solver,pts)
    n=boundary_matrix_size(solver,pts)
    rmin,rmax=_dlp_cheb_rmin_rmax(pts,orbits)
    plans1=Vector{ChebHankelPlanH}(undef,length(zj))
    @benchit timeit=timeit "DLP H1 plans" Threads.@threads for q in eachindex(zj)
        plans1[q]=plan_h(1,1,ComplexF64(zj[q]),rmin,rmax;npanels=n_panels_h,M=M_h)
    end
    @inbounds for q in eachindex(Tbufs)
        @assert size(Tbufs[q])==(n,n) "Tbufs[$q] has size $(size(Tbufs[q])), expected ($n,$n)"
    end
    @blas_1 begin
        if isnothing(orbits)
            @benchit timeit=timeit "DLP Chebyshev" compute_kernel_matrices_DLP_chebyshev!(Tbufs,pts,plans1;multithreaded=multithreaded)
            @benchit timeit=timeit "DLP Fredholm" assemble_fredholm_matrices!(Tbufs,pts)
        else
            @benchit timeit=timeit "DLP reduced Chebyshev" compute_kernel_matrices_DLP_chebyshev!(Tbufs,pts,orbits,plans1;multithreaded=multithreaded)
            @benchit timeit=timeit "DLP reduced Fredholm" assemble_fredholm_matrices!(Tbufs,pts,orbits)
        end
    end
    return nothing
end

"""
    construct_matrices_chebyshev_with_derivatives!(Tbufs::Vector{Matrix{Complex{T}}},dTbufs::Vector{Matrix{Complex{T}}},ddTbufs::Vector{Matrix{Complex{T}}},::Val{:dlp},solver::BoundaryIntegralMethod,pts::BoundaryPoints{T},zj::AbstractVector{Complex{T}};multithreaded::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false) where {T<:Real} → Nothing

Construct the full or exact symmetry-reduced DLP Fredholm matrices and their
first two wavenumber derivatives,

    A(k)=I-K(k)W,
    A'(k)=-K'(k)W,
    A''(k)=-K''(k)W,

for all complex wavenumbers in `zj`.

The actual plan construction and derivative matrix assembly are delegated to
the reusable derivative Chebyshev workspace API.

## Arguments
* `Tbufs::Vector{Matrix{Complex{T}}}`: Destination matrices for `A(k)`.
* `dTbufs::Vector{Matrix{Complex{T}}}`: Destination matrices for `A'(k)`.
* `ddTbufs::Vector{Matrix{Complex{T}}}`: Destination matrices for `A''(k)`.
* `::Val{:dlp}`: Selects the derivative-aware DLP backend.
* `solver::BoundaryIntegralMethod`: DLP solver.
* `pts::BoundaryPoints{T}`: Full boundary discretization.
* `zj::AbstractVector{Complex{T}}`: Complex contour wavenumbers.

## Keyword Arguments
* `multithreaded::Bool`: Enable threaded matrix assembly.
* `n_panels_h::Int`: Hankel Chebyshev panel count.
* `M_h::Int`: Hankel Chebyshev degree.
* `n_panels_j::Int`: Compatibility keyword; unused by this backend.
* `M_j::Int`: Compatibility keyword; unused by this backend.
* `timeit::Bool`: Enable timing diagnostics.

## Returns
* `nothing`.
"""
function construct_matrices_chebyshev_with_derivatives!(Tbufs::Vector{Matrix{Complex{T}}},dTbufs::Vector{Matrix{Complex{T}}},ddTbufs::Vector{Matrix{Complex{T}}},::Val{:dlp},solver::BoundaryIntegralMethod,pts::BoundaryPoints{T},zj::AbstractVector{Complex{T}};multithreaded::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false) where {T<:Real}
    @assert length(Tbufs)==length(zj)
    @assert length(dTbufs)==length(zj)
    @assert length(ddTbufs)==length(zj)
    n=boundary_matrix_size(solver,pts)
    @inbounds for q in eachindex(Tbufs)
        @assert size(Tbufs[q])==(n,n) "Tbufs[$q] has size $(size(Tbufs[q])), expected ($n,$n)"
        @assert size(dTbufs[q])==(n,n) "dTbufs[$q] has size $(size(dTbufs[q])), expected ($n,$n)"
        @assert size(ddTbufs[q])==(n,n) "ddTbufs[$q] has size $(size(ddTbufs[q])), expected ($n,$n)"
    end
    ws=build_derivative_chebyshev_workspace(solver,pts,zj;n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,timeit=timeit)
    @blas_1 @benchit timeit=timeit "DLP derivative Chebyshev" construct_matrices_chebyshev_with_derivatives!(Tbufs,dTbufs,ddTbufs,solver,pts,ws;multithreaded=multithreaded)
    return nothing
end

########################################
########### SOLVE-VECT BATCH ###########
########################################

"""
    adjoint_fredholm_matrix_from_bim_chebyshev!(A::AbstractMatrix{ComplexF64},K::AbstractMatrix{ComplexF64},pts::BoundaryPoints{T}) where {T<:Real} → A

Construct the full-boundary weighted formal-transpose Fredholm matrix directly
from a raw DLP kernel.

Since

    D=KW,

and

    A*=I-W⁻¹DᵀW,

the raw-kernel form simplifies to

    A*[i,j]=-K[j,i]ds[j]+δ_ij.

No complex conjugation is applied.

## Arguments
* `A::AbstractMatrix{ComplexF64}`: Destination formal-transpose matrix.
* `K::AbstractMatrix{ComplexF64}`: Raw full-boundary DLP kernel.
* `pts::BoundaryPoints{T}`: Full boundary discretization.

## Returns
* `A::AbstractMatrix{ComplexF64}`: Weighted formal-transpose Fredholm matrix.
"""
function adjoint_fredholm_matrix_from_bim_chebyshev!(A::AbstractMatrix{ComplexF64},K::AbstractMatrix{ComplexF64},pts::BoundaryPoints{T}) where {T<:Real}
    N=length(pts)
    @assert size(A)==(N,N)
    @assert size(K)==(N,N)
    fill!(A,0.0+0.0im)
    @inbounds for j in 1:N,i in 1:N
        A[i,j]=-K[j,i]*pts.ds[j]
    end
    @inbounds for i in 1:N
        A[i,i]+=1.0+0.0im
    end
    return A
end

"""
    adjoint_fredholm_matrix_from_bim_chebyshev!(A::AbstractMatrix{ComplexF64},K::AbstractMatrix{ComplexF64},pts::BoundaryPoints{T},orbits::SymmetryOrbitMap{T}) where {T<:Real} → A

Construct the weighted formal-transpose Fredholm matrix from a raw
symmetry-reduced DLP kernel.

For representative reduced source `b`,

    A*[a,b]=-K_red[b,a]ds[Ifund[b]]+δ_ab.

## Arguments
* `A::AbstractMatrix{ComplexF64}`: Destination reduced formal-transpose matrix.
* `K::AbstractMatrix{ComplexF64}`: Raw reduced DLP kernel.
* `pts::BoundaryPoints{T}`: Full boundary discretization.
* `orbits::SymmetryOrbitMap{T}`: Exact symmetry-orbit map.

## Returns
* `A::AbstractMatrix{ComplexF64}`: Reduced weighted formal-transpose Fredholm matrix.
"""
function adjoint_fredholm_matrix_from_bim_chebyshev!(A::AbstractMatrix{ComplexF64},K::AbstractMatrix{ComplexF64},pts::BoundaryPoints{T},orbits::SymmetryOrbitMap{T}) where {T<:Real}
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    @assert size(A)==(m,m)
    @assert size(K)==(m,m)
    fill!(A,0.0+0.0im)
    @inbounds for b in 1:m,a in 1:m
        A[a,b]=-K[b,a]*pts.ds[Ifund[b]]
    end
    @inbounds for a in 1:m
        A[a,a]+=1.0+0.0im
    end
    return A
end

"""
    solve_vect(solver::BoundaryIntegralMethod,billiard::Bi,basis::Ba,ks::Vector{T};batch_size::Int=40,multithreaded::Bool=true,use_chebyshev::Bool=true,cheb_tol::Real=1e-12,npanels_h_init::Int=15000,M_h_init::Int=5,sampling_points::Int=50000,max_iter::Int=20,grow_panels::Real=1.5,grow_M::Int=2,cheb_verbose::Bool=false,tol=1e-12,maxiter::Int=2000,krylovdim::Int=40) where {T<:Real,Ba<:AbstractHankelBasis,Bi<:BilliardGeometry.AbsBilliard} → Tuple{Vector{Vector{ComplexF64}},Vector{BoundaryPoints{T}}}

Compute BIM near-null boundary vectors for several real wavenumbers in
geometry batches.

When symmetry is active, the returned vectors live in the exact reduced
`SymmetryOrbitMap` basis.

The value-only Chebyshev pathway deliberately constructs only `H₁⁽¹⁾` plans;
it does not use the derivative workspace because `H₀⁽¹⁾` is unnecessary here.

## Arguments
* `solver::BoundaryIntegralMethod`: Boundary-integral solver.
* `billiard::Bi`: Billiard geometry.
* `basis::Ba`: Hankel basis compatibility object.
* `ks::Vector{T}`: Real target wavenumbers.

## Keyword Arguments
* `batch_size::Int`: Number of states per shared geometry batch.
* `multithreaded::Bool`: Enable threaded matrix assembly.
* `use_chebyshev::Bool`: Use Chebyshev-Hankel construction.
* `cheb_tol::Real`: Chebyshev validation tolerance.
* `npanels_h_init::Int`: Initial Hankel panel count.
* `M_h_init::Int`: Initial Hankel polynomial degree.
* `sampling_points::Int`: Number of validation samples.
* `max_iter::Int`: Maximum Chebyshev refinement iterations.
* `grow_panels::Real`: Panel-count growth factor.
* `grow_M::Int`: Polynomial-degree increment.
* `cheb_verbose::Bool`: Enable Chebyshev diagnostics.
* `tol`: Krylov eigensolver tolerance.
* `maxiter::Int`: Maximum Krylov iterations.
* `krylovdim::Int`: Krylov subspace dimension.

## Returns
* `us_all::Vector{Vector{ComplexF64}}`: Computed full or reduced near-null vectors.
* `pts_all::Vector{BoundaryPoints{T}}`: Boundary discretization associated with each vector.
"""
function solve_vect(solver::BoundaryIntegralMethod,billiard::Bi,basis::Ba,ks::Vector{T};batch_size::Int=40,multithreaded::Bool=true,use_chebyshev::Bool=true,cheb_tol::Real=1e-12,npanels_h_init::Int=15_000,M_h_init::Int=5,sampling_points::Int=50_000,max_iter::Int=20,grow_panels::Real=1.5,grow_M::Int=2,cheb_verbose::Bool=false,tol=1e-12,maxiter::Int=2000,krylovdim::Int=40) where {T<:Real,Ba<:AbstractHankelBasis,Bi<:BilliardGeometry.AbsBilliard}
    Nk=length(ks)
    us_all=Vector{Vector{ComplexF64}}(undef,Nk)
    pts_all=Vector{BoundaryPoints{T}}(undef,Nk)
    nb=_nbatches(Nk,batch_size)
    @showprogress "solve_vect BoundaryIntegralMethod" for ibatch in 1:nb
        i1=_batch_first(ibatch,batch_size)
        i2=_batch_last(ibatch,batch_size,Nk)
        inds=i1:i2
        kbatch=@view ks[inds]
        pts=evaluate_points(solver,billiard,maximum(kbatch))
        orbits=_dlp_symmetry_orbits(solver,pts)
        Nmat=boundary_matrix_size(solver,pts)
        if use_chebyshev
            zj=ComplexF64.(kbatch)
            cheb_out=chebyshev_params(solver,pts,zj;npanels_h_init=npanels_h_init,M_h_init=M_h_init,tol=cheb_tol,sampling_points=sampling_points,max_iter=max_iter,grow_panels=grow_panels,grow_M=grow_M,verbose=cheb_verbose)
            npanels_h=cheb_out[1]
            M_h=cheb_out[2]
            rmin,rmax=_dlp_cheb_rmin_rmax(pts,orbits)
            plans1=Vector{ChebHankelPlanH}(undef,length(zj))
            Threads.@threads for q in eachindex(zj)
                plans1[q]=plan_h(1,1,zj[q],rmin,rmax;npanels=npanels_h,M=M_h)
            end
            Ks=[zeros(ComplexF64,Nmat,Nmat) for _ in eachindex(zj)]
            A=Matrix{ComplexF64}(undef,Nmat,Nmat)
            if isnothing(orbits)
                compute_kernel_matrices_DLP_chebyshev!(Ks,pts,plans1;multithreaded=multithreaded)
                for (jlocal,jglobal) in enumerate(inds)
                    adjoint_fredholm_matrix_from_bim_chebyshev!(A,Ks[jlocal],pts)
                    _,u,_=smallest_nullvec_krylov!(A;nev=1,tol=tol,maxiter=maxiter,krylovdim=krylovdim)
                    us_all[jglobal]=ComplexF64.(u)
                    pts_all[jglobal]=pts
                end
            else
                compute_kernel_matrices_DLP_chebyshev!(Ks,pts,orbits,plans1;multithreaded=multithreaded)
                for (jlocal,jglobal) in enumerate(inds)
                    adjoint_fredholm_matrix_from_bim_chebyshev!(A,Ks[jlocal],pts,orbits)
                    _,u,_=smallest_nullvec_krylov!(A;nev=1,tol=tol,maxiter=maxiter,krylovdim=krylovdim)
                    us_all[jglobal]=ComplexF64.(u)
                    pts_all[jglobal]=pts
                end
            end
        else
            A=Matrix{Complex{T}}(undef,Nmat,Nmat)
            D=similar(A)
            for jglobal in inds
                @blas_1 adjoint_fredholm_matrix!(A,D,pts,solver.symmetry,ks[jglobal];multithreaded=multithreaded)
                _,u,_=smallest_nullvec_krylov!(A;nev=1,tol=tol,maxiter=maxiter,krylovdim=krylovdim)
                us_all[jglobal]=ComplexF64.(u)
                pts_all[jglobal]=pts
            end
        end
    end
    return us_all,pts_all
end