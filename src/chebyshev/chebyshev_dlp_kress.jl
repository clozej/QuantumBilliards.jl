################################################################################
# DLP-Kress Chebyshev matrix assembly
#
# For G_k(x,y)=(i/4)H₀⁽¹⁾(kr), r=|x-y|, the doubled DLP kernel is
#
#     K(k;x,y)=(i k/2)H₁⁽¹⁾(kr)c(x,y),   c(x,y)=((x-y)⋅n_y)/r.
#
# For boundary nodes x_i,x_j define
#
#     r_ij=|x_i-x_j|,   inner_ij=(x_i-x_j)⋅n_j,   c_ij=inner_ij/r_ij,
#     L_ij=logterm_ij,   R_ij=Rkress_ij,   w_j=quadrature weight,
#     A_ij=R_ij-w_j L_ij.
#
# For i≠j, the Kress-corrected discrete DLP element is
#
#     D_ij(k)=c_ij[-(k/2π)A_ij J₁(k r_ij)+(i k/2)w_j H₁⁽¹⁾(k r_ij)],
#
# while the analytic diagonal limit is
#
#     D_ii=w_i κ_i.
#
# Using d[kJ₁(kr)]/dk=krJ₀(kr) and d[kH₁⁽¹⁾(kr)]/dk=krH₀⁽¹⁾(kr),
#
#     D'_ij(k)=inner_ij k/(2π)[-A_ij J₀(k r_ij)+iπw_j H₀⁽¹⁾(k r_ij)],
#     D'_ii=0.
#
# With
#
#     u_ij=J₀(k r_ij)-k r_ij J₁(k r_ij),
#     v_ij=H₀⁽¹⁾(k r_ij)-k r_ij H₁⁽¹⁾(k r_ij),
#
# the second derivative is
#
#     D''_ij(k)=inner_ij/(2π)[-A_ij u_ij+iπw_j v_ij],   D''_ii=0.
#
# Hence value-only assembly requires H₁⁽¹⁾ and J₁, while derivative assembly
# requires H₀⁽¹⁾, H₁⁽¹⁾, J₀, and J₁. For complex k, Jν must be evaluated
# independently and cannot in general be replaced by real(Hν⁽¹⁾).
#
# The nonlinear Fredholm matrix and its derivatives are
#
#     F(k)=I-D(k),   F'(k)=-D'(k),   F''(k)=-D''(k).
#
# With symmetry, only fundamental-domain degrees of freedom are retained. For
# fundamental target a and reduced source b, let
#
#     i=Ifund[a],
#     j_l=fund_to_full[l,b],
#     χ_l=fund_to_scale[l,b].
#
# Symmetry reduction is performed by folding the complete discrete Kress
# operator:
#
#     Dred_ab(k)=Σ_l χ_l Dfull(i,j_l;k).
#
# Thus every source image is evaluated with exactly the same Kress
# product-integration formula as its corresponding entry of the full discrete
# operator. The same folding is applied to the first two k-derivatives:
#
#     Dred'_ab(k) =Σ_l χ_l Dfull'(i,j_l;k),
#     Dred''_ab(k)=Σ_l χ_l Dfull''(i,j_l;k).
################################################################################
#
# Reference:
#   R. Kress, "Boundary Integral Equations in Time-Harmonic Acoustic
#   Scattering," Mathl. Comput. Modelling 15(3-5), 229-243 (1991).
################################################################################
const INV_TWO_PI=1/(2*pi)
"""
    DLPKressBlockCache{T} where {T<:Real}

Geometry and interpolation cache for one DLP-Kress boundary component.

All fields depend only on the fixed boundary discretization and may therefore
be reused for many wavenumbers.

## Attributes
* `N::Int`: Number of full-boundary nodes.
* `R::Matrix{T}`: Pairwise distances.
* `invR::Matrix{T}`: Pairwise inverse distances.
* `inner::Matrix{T}`: Oriented DLP source-normal numerator matrix.
* `wi::Vector{T}`: Boundary quadrature parameter weights.
* `pidx::Matrix{Int32}`: Hankel Chebyshev panel indices.
* `tloc::Matrix{Float64}`: Hankel local Chebyshev coordinates.
* `pidxj::Matrix{Int32}`: Bessel-J Chebyshev panel indices.
* `tlocj::Matrix{Float64}`: Bessel-J local Chebyshev coordinates.
* `logterm::Matrix{T}`: Kress logarithmic split matrix.
* `kappa::Vector{T}`: Diagonal DLP Kress-limit values.
* `Rkress::Matrix{T}`: Kress logarithmic product-integration matrix.
"""
struct DLPKressBlockCache{T<:Real}
    N::Int
    R::Matrix{T}
    invR::Matrix{T}
    inner::Matrix{T}
    wi::Vector{T}
    pidx::Matrix{Int32}
    tloc::Matrix{Float64}
    pidxj::Matrix{Int32}
    tlocj::Matrix{Float64}
    logterm::Matrix{T}
    kappa::Vector{T}
    Rkress::Matrix{T}
end

"""
    DLPKressSystemCache{T} where {T<:Real}

Complete k-independent DLP-Kress Chebyshev geometry cache.

## Attributes
* `block::DLPKressBlockCache{T}`: Boundary geometry and interpolation-index cache.
* `rmin::Float64`: Lower Hankel interpolation radius.
* `rmax::Float64`: Upper Hankel/Bessel interpolation radius.
"""
struct DLPKressSystemCache{T<:Real}
    block::DLPKressBlockCache{T}
    rmin::Float64
    rmax::Float64
end

"""
    build_dlp_kress_block_cache(solver::Union{DLP_kress,DLP_kress_global_corners},pts::BoundaryPoints{T};npanels_h::Int=10000,npanels_j::Int=10000,M_h::Int=5,M_j::Int=5,pad::Tuple{T,T}=(T(0.95),T(1.05)),rmin_cheb::Union{Nothing,Float64}=nothing) where {T<:Real} → DLPKressSystemCache{T}

Build the k-independent geometry and interpolation cache used by the
DLP-Kress Chebyshev assembler.

The pairwise distances are assigned to Hankel and Bessel-J Chebyshev panels
once so repeated evaluations at many wavenumbers do not require repeated
panel searches.

## Arguments
* `solver::Union{DLP_kress,DLP_kress_global_corners}`: Smooth or globally graded DLP-Kress solver.
* `pts::BoundaryPoints{T}`: Full boundary discretization.

## Keyword Arguments
* `npanels_h::Int`: Number of Hankel reference panels.
* `npanels_j::Int`: Number of Bessel-J reference panels.
* `M_h::Int`: Hankel Chebyshev degree.
* `M_j::Int`: Bessel-J Chebyshev degree.
* `pad::Tuple{T,T}`: Multiplicative padding of the geometric radial interval.
* `rmin_cheb::Union{Nothing,Float64}`: Optional lower Hankel interpolation cutoff.

## Returns
* `cache::DLPKressSystemCache{T}`: Geometry and interpolation cache.
"""
function build_dlp_kress_block_cache(solver::Union{DLP_kress,DLP_kress_global_corners},pts::BoundaryPoints{T};npanels_h::Int=10000,npanels_j::Int=10000,M_h::Int=5,M_j::Int=5,pad::Tuple{T,T}=(T(0.95),T(1.05)),rmin_cheb::Union{Nothing,Float64}=nothing) where {T<:Real}
    graded=_is_dlp_kress_graded(solver,pts)
    G=boundary_geom_cache(pts,graded)
    N=length(pts.xy)
    R=copy(G.R)
    invR=copy(G.invR)
    inner=copy(G.inner)
    wi=copy(pts.ws)
    logterm=copy(G.logterm)
    kappa=copy(G.kappa)
    Rkress=zeros(T,N,N)
    kress_R!(Rkress)
    rmin0=typemax(T)
    rmax0=zero(T)
    @inbounds for j in 1:N,i in 1:N
        i==j&&continue
        r=R[i,j]
        if isfinite(r)&&r>eps(T)
            r<rmin0&&(rmin0=r)
            r>rmax0&&(rmax0=r)
        end
    end
    @assert isfinite(rmin0)&&rmax0>zero(T)
    rrmin=Float64(pad[1]*rmin0)
    rrmax=Float64(pad[2]*rmax0)
    rminloc=isnothing(rmin_cheb) ? rrmin : max(Float64(rmin_cheb),rrmin)
    pidx=Matrix{Int32}(undef,N,N)
    tloc=Matrix{Float64}(undef,N,N)
    pidxj=Matrix{Int32}(undef,N,N)
    tlocj=Matrix{Float64}(undef,N,N)
    refh=plan_h(0,1,1.0+0im,rminloc,rrmax;npanels=npanels_h,M=M_h)
    refj=plan_j(1,1.0+0im,0.0,rrmax;npanels=npanels_j,M=M_j)
    pansh=refh.panels
    pansj=refj.panels
    @inbounds for j in 1:N,i in 1:N
        if i==j
            pidx[i,j]=Int32(1)
            tloc[i,j]=0.0
            pidxj[i,j]=Int32(1)
            tlocj[i,j]=0.0
            continue
        end
        r=Float64(R[i,j])
        if r<rminloc
            pidx[i,j]=Int32(0)
            tloc[i,j]=0.0
        else
            p=_find_panel(refh,r)
            P=pansh[p]
            pidx[i,j]=Int32(p)
            tloc[i,j]=(2r-(P.b+P.a))/(P.b-P.a)
        end
        p=_find_panel(refj,r)
        P=pansj[p]
        pidxj[i,j]=Int32(p)
        tlocj[i,j]=(2r-(P.b+P.a))/(P.b-P.a)
    end
    block=DLPKressBlockCache{T}(N,R,invR,inner,wi,pidx,tloc,pidxj,tlocj,logterm,kappa,Rkress)
    return DLPKressSystemCache{T}(block,rminloc,rrmax)
end

####################################
######## CHEBYSHEV PLANS ###########
####################################

"""
    build_dlp_kress_plans_h1_j1(ks::AbstractVector{<:Number},rmin::Float64,rmax::Float64;npanels_h::Int=10000,npanels_j::Int=10000,M_h::Int=5,M_j::Int=5) → Tuple{Vector{ChebHankelPlanH},Vector{ChebJPlan}}

Build value-only Chebyshev plans for `H₁⁽¹⁾(kr)` and `J₁(kr)` for all
wavenumbers in `ks`.
"""
function build_dlp_kress_plans_h1_j1(ks::AbstractVector{<:Number},rmin::Float64,rmax::Float64;npanels_h::Int=10000,npanels_j::Int=10000,M_h::Int=5,M_j::Int=5)
    Mk=length(ks)
    plans1=Vector{ChebHankelPlanH}(undef,Mk)
    plansj1=Vector{ChebJPlan}(undef,Mk)
    if Threads.nthreads()==1||Mk==1
        @inbounds for m in 1:Mk
            k=ComplexF64(ks[m])
            plans1[m]=plan_h(1,1,k,rmin,rmax;npanels=npanels_h,M=M_h)
            plansj1[m]=plan_j(1,k,0.0,rmax;npanels=npanels_j,M=M_j)
        end
    else
        Threads.@threads for m in 1:Mk
            k=ComplexF64(ks[m])
            plans1[m]=plan_h(1,1,k,rmin,rmax;npanels=npanels_h,M=M_h)
            plansj1[m]=plan_j(1,k,0.0,rmax;npanels=npanels_j,M=M_j)
        end
    end
    return plans1,plansj1
end

"""
    build_dlp_kress_plans_h0_h1_j0_j1(ks::AbstractVector{<:Number},rmin::Float64,rmax::Float64;npanels_h::Int=10000,npanels_j::Int=10000,M_h::Int=5,M_j::Int=5) → Tuple{Vector{ChebHankelPlanH},Vector{ChebHankelPlanH},Vector{ChebJPlan},Vector{ChebJPlan}}

Build derivative-aware Chebyshev plans for `H₀⁽¹⁾`, `H₁⁽¹⁾`, `J₀`, and
`J₁` for all wavenumbers in `ks`.
"""
function build_dlp_kress_plans_h0_h1_j0_j1(ks::AbstractVector{<:Number},rmin::Float64,rmax::Float64;npanels_h::Int=10000,npanels_j::Int=10000,M_h::Int=5,M_j::Int=5)
    Mk=length(ks)
    plans0=Vector{ChebHankelPlanH}(undef,Mk)
    plans1=Vector{ChebHankelPlanH}(undef,Mk)
    plansj0=Vector{ChebJPlan}(undef,Mk)
    plansj1=Vector{ChebJPlan}(undef,Mk)
    if Threads.nthreads()==1||Mk==1
        @inbounds for m in 1:Mk
            k=ComplexF64(ks[m])
            plans0[m]=plan_h(0,1,k,rmin,rmax;npanels=npanels_h,M=M_h)
            plans1[m]=plan_h(1,1,k,rmin,rmax;npanels=npanels_h,M=M_h)
            plansj0[m]=plan_j(0,k,0.0,rmax;npanels=npanels_j,M=M_j)
            plansj1[m]=plan_j(1,k,0.0,rmax;npanels=npanels_j,M=M_j)
        end
    else
        Threads.@threads for m in 1:Mk
            k=ComplexF64(ks[m])
            plans0[m]=plan_h(0,1,k,rmin,rmax;npanels=npanels_h,M=M_h)
            plans1[m]=plan_h(1,1,k,rmin,rmax;npanels=npanels_h,M=M_h)
            plansj0[m]=plan_j(0,k,0.0,rmax;npanels=npanels_j,M=M_j)
            plansj1[m]=plan_j(1,k,0.0,rmax;npanels=npanels_j,M=M_j)
        end
    end
    return plans0,plans1,plansj0,plansj1
end

####################################
######## BESSEL WORKSPACES #########
####################################

"""
    DLPKressH1J1BesselWorkspace

Thread-local temporary buffers for value-only DLP-Kress special-function
evaluation.

## Attributes
* `h1_tls::Vector{Vector{ComplexF64}}`: Thread-local `H₁⁽¹⁾` values.
* `j1_tls::Vector{Vector{ComplexF64}}`: Thread-local `J₁` values.
"""
struct DLPKressH1J1BesselWorkspace
    h1_tls::Vector{Vector{ComplexF64}}
    j1_tls::Vector{Vector{ComplexF64}}
end

"""
    DLPKressH1J1BesselWorkspace(Mk::Int;ntls::Int=Threads.nthreads()) → DLPKressH1J1BesselWorkspace

Allocate thread-local buffers for value-only DLP-Kress Chebyshev evaluation.
"""
DLPKressH1J1BesselWorkspace(Mk::Int;ntls::Int=Threads.nthreads())=DLPKressH1J1BesselWorkspace([Vector{ComplexF64}(undef,Mk) for _ in 1:ntls],[Vector{ComplexF64}(undef,Mk) for _ in 1:ntls])

"""
    DLPKressH0H1J0J1BesselWorkspace

Thread-local temporary buffers for derivative-aware DLP-Kress
special-function evaluation.

## Attributes
* `h0_tls::Vector{Vector{ComplexF64}}`: Thread-local `H₀⁽¹⁾` values.
* `h1_tls::Vector{Vector{ComplexF64}}`: Thread-local `H₁⁽¹⁾` values.
* `j0_tls::Vector{Vector{ComplexF64}}`: Thread-local `J₀` values.
* `j1_tls::Vector{Vector{ComplexF64}}`: Thread-local `J₁` values.
"""
struct DLPKressH0H1J0J1BesselWorkspace
    h0_tls::Vector{Vector{ComplexF64}}
    h1_tls::Vector{Vector{ComplexF64}}
    j0_tls::Vector{Vector{ComplexF64}}
    j1_tls::Vector{Vector{ComplexF64}}
end

"""
    DLPKressH0H1J0J1BesselWorkspace(Mk::Int;ntls::Int=Threads.nthreads()) → DLPKressH0H1J0J1BesselWorkspace

Allocate thread-local buffers for derivative-aware DLP-Kress Chebyshev
evaluation.
"""
DLPKressH0H1J0J1BesselWorkspace(Mk::Int;ntls::Int=Threads.nthreads())=DLPKressH0H1J0J1BesselWorkspace([Vector{ComplexF64}(undef,Mk) for _ in 1:ntls],[Vector{ComplexF64}(undef,Mk) for _ in 1:ntls],[Vector{ComplexF64}(undef,Mk) for _ in 1:ntls],[Vector{ComplexF64}(undef,Mk) for _ in 1:ntls])

####################################
######## CHEBYSHEV WORKSPACES ######
####################################

"""
    DLPKressH1J1ChebWorkspace{T,MatT} where {T<:Real,MatT<:AbstractMatrix{T}}

Reusable value-only DLP-Kress Chebyshev workspace for full-boundary assembly.
"""
struct DLPKressH1J1ChebWorkspace{T<:Real,MatT<:AbstractMatrix{T}}
    direct::DLPKressWorkspace{T,MatT}
    block_cache::DLPKressSystemCache{T}
    plans1::Vector{ChebHankelPlanH}
    plansj1::Vector{ChebJPlan}
    bessel_ws::DLPKressH1J1BesselWorkspace
    ks::Vector{ComplexF64}
    Mk::Int
end

"""
    DLPKressH0H1J0J1ChebWorkspace{T,MatT} where {T<:Real,MatT<:AbstractMatrix{T}}

Reusable derivative-aware DLP-Kress Chebyshev workspace for full-boundary
assembly.
"""
struct DLPKressH0H1J0J1ChebWorkspace{T<:Real,MatT<:AbstractMatrix{T}}
    direct::DLPKressWorkspace{T,MatT}
    block_cache::DLPKressSystemCache{T}
    plans0::Vector{ChebHankelPlanH}
    plans1::Vector{ChebHankelPlanH}
    plansj0::Vector{ChebJPlan}
    plansj1::Vector{ChebJPlan}
    bessel_ws::DLPKressH0H1J0J1BesselWorkspace
    ks::Vector{ComplexF64}
    Mk::Int
end

"""
    DLPKressReducedH1J1ChebWorkspace{T,MatT} where {T<:Real,MatT<:AbstractMatrix{T}}

Value-only symmetry-reduced DLP-Kress Chebyshev workspace.
"""
struct DLPKressReducedH1J1ChebWorkspace{T<:Real,MatT<:AbstractMatrix{T}}
    direct::DLPKressReducedWorkspace{T,MatT}
    fullcheb::DLPKressH1J1ChebWorkspace{T,MatT}
end

"""
    DLPKressReducedH0H1J0J1ChebWorkspace{T,MatT} where {T<:Real,MatT<:AbstractMatrix{T}}

Derivative-aware symmetry-reduced DLP-Kress Chebyshev workspace.
"""
struct DLPKressReducedH0H1J0J1ChebWorkspace{T<:Real,MatT<:AbstractMatrix{T}}
    direct::DLPKressReducedWorkspace{T,MatT}
    fullcheb::DLPKressH0H1J0J1ChebWorkspace{T,MatT}
end

const DLPKressValueChebWorkspace=Union{DLPKressH1J1ChebWorkspace,DLPKressReducedH1J1ChebWorkspace}
const DLPKressDerivativeChebWorkspace=Union{DLPKressH0H1J0J1ChebWorkspace,DLPKressReducedH0H1J0J1ChebWorkspace}

@inline _cheb_workspace_dim(ws::DLPKressH1J1ChebWorkspace)=ws.block_cache.block.N
@inline _cheb_workspace_dim(ws::DLPKressH0H1J0J1ChebWorkspace)=ws.block_cache.block.N
@inline _cheb_workspace_dim(ws::DLPKressReducedH1J1ChebWorkspace)=fundamental_size(ws.direct.orbits)
@inline _cheb_workspace_dim(ws::DLPKressReducedH0H1J0J1ChebWorkspace)=fundamental_size(ws.direct.orbits)

@inline _cheb_workspace_length(ws::DLPKressH1J1ChebWorkspace)=ws.Mk
@inline _cheb_workspace_length(ws::DLPKressH0H1J0J1ChebWorkspace)=ws.Mk
@inline _cheb_workspace_length(ws::DLPKressReducedH1J1ChebWorkspace)=ws.fullcheb.Mk
@inline _cheb_workspace_length(ws::DLPKressReducedH0H1J0J1ChebWorkspace)=ws.fullcheb.Mk

####################################
######## WORKSPACE BUILDERS ########
####################################

"""
    build_dlp_kress_h1_j1_cheb_workspace(solver::Union{DLP_kress{T},DLP_kress_global_corners{T}},pts::BoundaryPoints{T},direct::DLPKressWorkspace{T,MatT},ks::Vector{ComplexF64};npanels_h::Int=10000,npanels_j::Int=2000,M_h::Int=5,M_j::Int=5,pad::Tuple{T,T}=(T(0.95),T(1.05)),rmin_cheb::Union{Nothing,Float64}=nothing,ntls::Int=Threads.nthreads()) where {T<:Real,MatT<:AbstractMatrix{T}} → DLPKressH1J1ChebWorkspace{T,MatT}

Build the value-only full-boundary DLP-Kress Chebyshev workspace.
"""
function build_dlp_kress_h1_j1_cheb_workspace(solver::Union{DLP_kress{T},DLP_kress_global_corners{T}},pts::BoundaryPoints{T},direct::DLPKressWorkspace{T,MatT},ks::Vector{ComplexF64};npanels_h::Int=10000,npanels_j::Int=2000,M_h::Int=5,M_j::Int=5,pad::Tuple{T,T}=(T(0.95),T(1.05)),rmin_cheb::Union{Nothing,Float64}=nothing,ntls::Int=Threads.nthreads()) where {T<:Real,MatT<:AbstractMatrix{T}}
    cache=build_dlp_kress_block_cache(solver,pts;npanels_h=npanels_h,npanels_j=npanels_j,M_h=M_h,M_j=M_j,pad=pad,rmin_cheb=rmin_cheb)
    plans1,plansj1=build_dlp_kress_plans_h1_j1(ks,cache.rmin,cache.rmax;npanels_h=npanels_h,npanels_j=npanels_j,M_h=M_h,M_j=M_j)
    bessel_ws=DLPKressH1J1BesselWorkspace(length(ks);ntls=ntls)
    return DLPKressH1J1ChebWorkspace{T,MatT}(direct,cache,plans1,plansj1,bessel_ws,ks,length(ks))
end

"""
    build_dlp_kress_h0_h1_j0_j1_cheb_workspace(solver::Union{DLP_kress{T},DLP_kress_global_corners{T}},pts::BoundaryPoints{T},direct::DLPKressWorkspace{T,MatT},ks::Vector{ComplexF64};npanels_h::Int=10000,npanels_j::Int=2000,M_h::Int=5,M_j::Int=5,pad::Tuple{T,T}=(T(0.95),T(1.05)),rmin_cheb::Union{Nothing,Float64}=nothing,ntls::Int=Threads.nthreads()) where {T<:Real,MatT<:AbstractMatrix{T}} → DLPKressH0H1J0J1ChebWorkspace{T,MatT}

Build the derivative-aware full-boundary DLP-Kress Chebyshev workspace.
"""
function build_dlp_kress_h0_h1_j0_j1_cheb_workspace(solver::Union{DLP_kress{T},DLP_kress_global_corners{T}},pts::BoundaryPoints{T},direct::DLPKressWorkspace{T,MatT},ks::Vector{ComplexF64};npanels_h::Int=10000,npanels_j::Int=2000,M_h::Int=5,M_j::Int=5,pad::Tuple{T,T}=(T(0.95),T(1.05)),rmin_cheb::Union{Nothing,Float64}=nothing,ntls::Int=Threads.nthreads()) where {T<:Real,MatT<:AbstractMatrix{T}}
    cache=build_dlp_kress_block_cache(solver,pts;npanels_h=npanels_h,npanels_j=npanels_j,M_h=M_h,M_j=M_j,pad=pad,rmin_cheb=rmin_cheb)
    plans0,plans1,plansj0,plansj1=build_dlp_kress_plans_h0_h1_j0_j1(ks,cache.rmin,cache.rmax;npanels_h=npanels_h,npanels_j=npanels_j,M_h=M_h,M_j=M_j)
    bessel_ws=DLPKressH0H1J0J1BesselWorkspace(length(ks);ntls=ntls)
    return DLPKressH0H1J0J1ChebWorkspace{T,MatT}(direct,cache,plans0,plans1,plansj0,plansj1,bessel_ws,ks,length(ks))
end

"""
    build_dlp_kress_h1_j1_cheb_workspace(solver::Union{DLP_kress{T},DLP_kress_global_corners{T}},pts::BoundaryPoints{T},direct::DLPKressReducedWorkspace{T,MatT},ks::Vector{ComplexF64};kwargs...) where {T<:Real,MatT<:AbstractMatrix{T}} → DLPKressReducedH1J1ChebWorkspace{T,MatT}

Build the value-only symmetry-reduced DLP-Kress Chebyshev workspace.
"""
function build_dlp_kress_h1_j1_cheb_workspace(solver::Union{DLP_kress{T},DLP_kress_global_corners{T}},pts::BoundaryPoints{T},direct::DLPKressReducedWorkspace{T,MatT},ks::Vector{ComplexF64};kwargs...) where {T<:Real,MatT<:AbstractMatrix{T}}
    fullcheb=build_dlp_kress_h1_j1_cheb_workspace(solver,pts,direct.full,ks;kwargs...)
    return DLPKressReducedH1J1ChebWorkspace{T,MatT}(direct,fullcheb)
end

"""
    build_dlp_kress_h0_h1_j0_j1_cheb_workspace(solver::Union{DLP_kress{T},DLP_kress_global_corners{T}},pts::BoundaryPoints{T},direct::DLPKressReducedWorkspace{T,MatT},ks::Vector{ComplexF64};kwargs...) where {T<:Real,MatT<:AbstractMatrix{T}} → DLPKressReducedH0H1J0J1ChebWorkspace{T,MatT}

Build the derivative-aware symmetry-reduced DLP-Kress Chebyshev workspace.
"""
function build_dlp_kress_h0_h1_j0_j1_cheb_workspace(solver::Union{DLP_kress{T},DLP_kress_global_corners{T}},pts::BoundaryPoints{T},direct::DLPKressReducedWorkspace{T,MatT},ks::Vector{ComplexF64};kwargs...) where {T<:Real,MatT<:AbstractMatrix{T}}
    fullcheb=build_dlp_kress_h0_h1_j0_j1_cheb_workspace(solver,pts,direct.full,ks;kwargs...)
    return DLPKressReducedH0H1J0J1ChebWorkspace{T,MatT}(direct,fullcheb)
end

################################################################################
################ COMMON DERIVATIVE CHEBYSHEV WORKSPACE API #####################
################################################################################

"""
    build_derivative_chebyshev_workspace(solver::Union{DLP_kress{T},DLP_kress_global_corners{T}},pts::BoundaryPoints{T},ks::AbstractVector{<:Number};n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,pad::Tuple{T,T}=(T(0.95),T(1.05)),rmin_cheb::Union{Nothing,Float64}=nothing,timeit::Bool=false) where {T<:Real} → DLPKressDerivativeChebWorkspace

Build the reusable derivative-aware DLP-Kress Chebyshev workspace.

The concrete returned workspace is selected automatically by
`build_dlp_kress_workspace`: a full workspace without symmetry and a reduced
workspace when an exact `SymmetryOrbitMap` is active.

Higher-level algorithms therefore do not need to know whether DLP-Kress uses
full or reduced matrices.

## Arguments
* `solver::Union{DLP_kress{T},DLP_kress_global_corners{T}}`: DLP-Kress solver.
* `pts::BoundaryPoints{T}`: Full boundary discretization.
* `ks::AbstractVector{<:Number}`: Wavenumbers represented by the workspace.

## Keyword Arguments
* `n_panels_h::Int`: Hankel Chebyshev panel count.
* `M_h::Int`: Hankel Chebyshev degree.
* `n_panels_j::Int`: Bessel-J Chebyshev panel count.
* `M_j::Int`: Bessel-J Chebyshev degree.
* `pad::Tuple{T,T}`: Radial interpolation padding.
* `rmin_cheb::Union{Nothing,Float64}`: Optional lower Hankel interpolation cutoff.
* `timeit::Bool`: Enable timing diagnostics.

## Returns
* `ws::DLPKressDerivativeChebWorkspace`: Full or reduced derivative workspace.
"""
function build_derivative_chebyshev_workspace(solver::Union{DLP_kress{T},DLP_kress_global_corners{T}},pts::BoundaryPoints{T},ks::AbstractVector{<:Number};n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,pad::Tuple{T,T}=(T(0.95),T(1.05)),rmin_cheb::Union{Nothing,Float64}=nothing,timeit::Bool=false) where {T<:Real}
    zks=ComplexF64.(ks)
    @benchit timeit=timeit "DLP-Kress direct workspace" directws=build_dlp_kress_workspace(solver,pts)
    @benchit timeit=timeit "DLP-Kress derivative Chebyshev workspace" chebws=build_dlp_kress_h0_h1_j0_j1_cheb_workspace(solver,pts,directws,zks;npanels_h=n_panels_h,npanels_j=n_panels_j,M_h=M_h,M_j=M_j,pad=pad,rmin_cheb=rmin_cheb)
    return chebws
end

"""
    construct_matrices_chebyshev_with_derivatives!(Fs::Vector{<:AbstractMatrix{ComplexF64}},F1s::Vector{<:AbstractMatrix{ComplexF64}},F2s::Vector{<:AbstractMatrix{ComplexF64}},solver::Union{DLP_kress,DLP_kress_global_corners},pts::BoundaryPoints{T},ws::DLPKressDerivativeChebWorkspace;multithreaded::Bool=true) where {T<:Real} → Nothing

Construct all DLP-Kress Fredholm matrices and their first two wavenumber
derivatives from a reusable derivative Chebyshev workspace.

Full/reduced dispatch is determined entirely by the concrete workspace type.
"""
function construct_matrices_chebyshev_with_derivatives!(Fs::Vector{<:AbstractMatrix{ComplexF64}},F1s::Vector{<:AbstractMatrix{ComplexF64}},F2s::Vector{<:AbstractMatrix{ComplexF64}},solver::Union{DLP_kress,DLP_kress_global_corners},pts::BoundaryPoints{T},ws::DLPKressDerivativeChebWorkspace;multithreaded::Bool=true) where {T<:Real}
    Mk=_cheb_workspace_length(ws)
    @assert length(Fs)==Mk
    @assert length(F1s)==Mk
    @assert length(F2s)==Mk
    n=_cheb_workspace_dim(ws)
    @inbounds for q in 1:Mk
        @assert size(Fs[q])==(n,n) "Fs[$q] has size $(size(Fs[q])), expected ($n,$n)"
        @assert size(F1s[q])==(n,n) "F1s[$q] has size $(size(F1s[q])), expected ($n,$n)"
        @assert size(F2s[q])==(n,n) "F2s[$q] has size $(size(F2s[q])), expected ($n,$n)"
    end
    construct_dlp_kress_matrices_derivatives_chebyshev!(Fs,F1s,F2s,pts,ws;multithreaded=multithreaded)
    return nothing
end

@inline function _single_derivative_chebyshev_workspace(ws::DLPKressH0H1J0J1ChebWorkspace{T,MatT},idx::Int) where {T<:Real,MatT<:AbstractMatrix{T}}
    checkbounds(ws.ks,idx)
    ntls=length(ws.bessel_ws.h0_tls)
    return DLPKressH0H1J0J1ChebWorkspace{T,MatT}(ws.direct,ws.block_cache,[ws.plans0[idx]],[ws.plans1[idx]],[ws.plansj0[idx]],[ws.plansj1[idx]],DLPKressH0H1J0J1BesselWorkspace(1;ntls=ntls),[ws.ks[idx]],1)
end

@inline function _single_derivative_chebyshev_workspace(ws::DLPKressReducedH0H1J0J1ChebWorkspace{T,MatT},idx::Int) where {T<:Real,MatT<:AbstractMatrix{T}}
    fullcheb=_single_derivative_chebyshev_workspace(ws.fullcheb,idx)
    return DLPKressReducedH0H1J0J1ChebWorkspace{T,MatT}(ws.direct,fullcheb)
end

"""
    construct_matrix_chebyshev_with_derivatives_at!(F::Matrix{ComplexF64},F1::Matrix{ComplexF64},F2::Matrix{ComplexF64},solver::Union{DLP_kress,DLP_kress_global_corners},pts::BoundaryPoints{T},ws::DLPKressDerivativeChebWorkspace,idx::Int;multithreaded::Bool=true) where {T<:Real} → Nothing

Construct one DLP-Kress Fredholm matrix and its first two wavenumber
derivatives from entry `idx` of a reusable derivative Chebyshev workspace.

No geometry cache or Chebyshev plan is rebuilt. A one-wavenumber wrapper around
the existing cached plans is constructed and passed to the ordinary DLP-Kress
assembler.
"""
function construct_matrix_chebyshev_with_derivatives_at!(F::Matrix{ComplexF64},F1::Matrix{ComplexF64},F2::Matrix{ComplexF64},solver::Union{DLP_kress,DLP_kress_global_corners},pts::BoundaryPoints{T},ws::DLPKressDerivativeChebWorkspace,idx::Int;multithreaded::Bool=true) where {T<:Real}
    n=_cheb_workspace_dim(ws)
    @assert size(F)==(n,n)
    @assert size(F1)==(n,n)
    @assert size(F2)==(n,n)
    ws1=_single_derivative_chebyshev_workspace(ws,idx)
    construct_dlp_kress_matrices_derivatives_chebyshev!([F],[F1],[F2],pts,ws1;multithreaded=multithreaded)
    return nothing
end

####################################
######## INNER EVALUATORS ##########
####################################

@inline function _h0_h1_j0_j1_at_pidx_t!(h0vals::AbstractVector{ComplexF64},h1vals::AbstractVector{ComplexF64},j0vals::AbstractVector{ComplexF64},j1vals::AbstractVector{ComplexF64},pidx_h::Int32,t_h::Float64,pidx_j::Int32,t_j::Float64,r::Float64,plans0::AbstractVector{ChebHankelPlanH},plans1::AbstractVector{ChebHankelPlanH},plansj0::AbstractVector{ChebJPlan},plansj1::AbstractVector{ChebJPlan})
    h0_h1_j0_j1_multi_ks_at_r!(h0vals,h1vals,j0vals,j1vals,plans0,plans1,plansj0,plansj1,pidx_h,t_h,pidx_j,t_j,r)
    return nothing
end

@inline function _h1_j1_at_pidx_t!(h1vals::AbstractVector{ComplexF64},j1vals::AbstractVector{ComplexF64},pidx_h::Int32,t_h::Float64,pidx_j::Int32,t_j::Float64,r::Float64,plans1::AbstractVector{ChebHankelPlanH},plansj1::AbstractVector{ChebJPlan})
    h1_j1_multi_ks_at_r!(h1vals,j1vals,plans1,plansj1,pidx_h,t_h,pidx_j,t_j,r)
    return nothing
end

####################################
######## RAW VALUE ASSEMBLY ########
####################################

function _construct_dlp_kress_matrices_chebyshev!(Ds::Vector{Matrix{ComplexF64}},pts::BoundaryPoints{T},ws::DLPKressH1J1ChebWorkspace{T};multithreaded::Bool=true) where {T<:Real}
    Mk=ws.Mk
    blk=ws.block_cache.block
    N=blk.N
    ks=ws.ks
    @inbounds for q in 1:Mk
        fill!(Ds[q],0.0+0.0im)
    end
    αL1s=Vector{ComplexF64}(undef,Mk)
    αL2s=Vector{ComplexF64}(undef,Mk)
    @inbounds for q in 1:Mk
        k=ks[q]
        αL1s[q]=-k*INV_TWO_PI
        αL2s[q]=0.5im*k
    end
    h1_tls=ws.bessel_ws.h1_tls
    j1_tls=ws.bessel_ws.j1_tls
    plans1=ws.plans1
    plansj1=ws.plansj1
    R=blk.R
    invR=blk.invR
    logterm=blk.logterm
    inner=blk.inner
    Rkress=blk.Rkress
    wi=blk.wi
    pidx=blk.pidx
    tloc=blk.tloc
    pidxj=blk.pidxj
    tlocj=blk.tlocj
    kappa=blk.kappa
    @inbounds for q in 1:Mk,i in 1:N
        Ds[q][i,i]=ComplexF64(wi[i]*kappa[i],0.0)
    end
    @use_threads multithreading=multithreaded for j in 2:N
        tid=Threads.threadid()
        h1vals=h1_tls[tid]
        j1vals=j1_tls[tid]
        @inbounds for i in 1:j-1
            r=R[i,j]
            _h1_j1_at_pidx_t!(h1vals,j1vals,pidx[i,j],tloc[i,j],pidxj[i,j],tlocj[i,j],r,plans1,plansj1)
            rinv=invR[i,j]
            lt=logterm[i,j]
            Rij=Rkress[i,j]
            cDij=inner[i,j]*rinv
            cDji=inner[j,i]*rinv
            c1ij=Rij*cDij
            c2ij=wi[j]*cDij
            c3ij=wi[j]*lt*cDij
            c1ji=Rij*cDji
            c2ji=wi[i]*cDji
            c3ji=wi[i]*lt*cDji
            for q in 1:Mk
                l1=αL1s[q]*j1vals[q]
                h1=αL2s[q]*h1vals[q]
                Ds[q][i,j]=c1ij*l1+c2ij*h1-c3ij*l1
                Ds[q][j,i]=c1ji*l1+c2ji*h1-c3ji*l1
            end
        end
    end
    return nothing
end

function _construct_dlp_kress_matrices_chebyshev!(Ds::Vector{Matrix{ComplexF64}},pts::BoundaryPoints{T},rws::DLPKressReducedH1J1ChebWorkspace{T};multithreaded::Bool=true) where {T<:Real}
    fullws=rws.fullcheb
    blk=fullws.block_cache.block
    orbits=rws.direct.orbits
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    ng=orbit_size(orbits)
    Mk=fullws.Mk
    ks=fullws.ks
    @inbounds for q in 1:Mk
        @assert size(Ds[q])==(m,m)
        fill!(Ds[q],0.0+0.0im)
    end
    αL1s=Vector{ComplexF64}(undef,Mk)
    αL2s=Vector{ComplexF64}(undef,Mk)
    @inbounds for q in 1:Mk
        αL1s[q]=-ks[q]*INV_TWO_PI
        αL2s[q]=0.5im*ks[q]
    end
    ntls=length(fullws.bessel_ws.h1_tls)
    acc_tls=[Vector{ComplexF64}(undef,Mk) for _ in 1:ntls]
    h1_tls=fullws.bessel_ws.h1_tls
    j1_tls=fullws.bessel_ws.j1_tls
    plans1=fullws.plans1
    plansj1=fullws.plansj1
    R=blk.R
    invR=blk.invR
    logterm=blk.logterm
    inner=blk.inner
    Rkress=blk.Rkress
    wi=blk.wi
    pidx=blk.pidx
    tloc=blk.tloc
    pidxj=blk.pidxj
    tlocj=blk.tlocj
    kappa=blk.kappa
    @use_threads multithreading=multithreaded for b in 1:m
        tid=Threads.threadid()
        h1vals=h1_tls[tid]
        j1vals=j1_tls[tid]
        acc=acc_tls[tid]
        @inbounds for a in 1:m
            fill!(acc,0.0+0.0im)
            i=Ifund[a]
            for l in 1:ng
                j=orbits.fund_to_full[l,b]
                χ=ComplexF64(orbits.fund_to_scale[l,b])
                if i==j
                    d0=ComplexF64(wi[i]*kappa[i],0.0)
                    for q in 1:Mk
                        acc[q]+=χ*d0
                    end
                else
                    r=R[i,j]
                    rinv=invR[i,j]
                    lt=logterm[i,j]
                    inn=inner[i,j]
                    Rij=Rkress[i,j]
                    wj=wi[j]
                    _h1_j1_at_pidx_t!(h1vals,j1vals,pidx[i,j],tloc[i,j],pidxj[i,j],tlocj[i,j],r,plans1,plansj1)
                    c1=Rij*inn*rinv
                    c2=wj*inn*rinv
                    c3=c2*lt
                    for q in 1:Mk
                        l1=αL1s[q]*j1vals[q]
                        dval=c1*l1+c2*αL2s[q]*h1vals[q]-c3*l1
                        acc[q]+=χ*dval
                    end
                end
            end
            for q in 1:Mk
                Ds[q][a,b]=acc[q]
            end
        end
    end
    return nothing
end

####################################
###### RAW DERIVATIVE ASSEMBLY #####
####################################

function _construct_dlp_kress_matrices_derivatives_chebyshev!(Ds::Vector{<:AbstractMatrix{ComplexF64}},D1s::Vector{<:AbstractMatrix{ComplexF64}},D2s::Vector{<:AbstractMatrix{ComplexF64}},pts::BoundaryPoints{T},ws::DLPKressH0H1J0J1ChebWorkspace{T};multithreaded::Bool=true) where {T<:Real}
    Mk=ws.Mk
    blk=ws.block_cache.block
    N=blk.N
    ks=ws.ks
    R=blk.R
    invR=blk.invR
    logterm=blk.logterm
    inner=blk.inner
    Rkress=blk.Rkress
    wi=blk.wi
    kappa=blk.kappa
    pidx=blk.pidx
    tloc=blk.tloc
    pidxj=blk.pidxj
    tlocj=blk.tlocj
    plans0=ws.plans0
    plans1=ws.plans1
    plansj0=ws.plansj0
    plansj1=ws.plansj1
    αL1s=Vector{ComplexF64}(undef,Mk)
    αL2s=Vector{ComplexF64}(undef,Mk)
    kcs=Vector{ComplexF64}(undef,Mk)
    @inbounds for q in 1:Mk
        fill!(Ds[q],0.0+0.0im)
        fill!(D1s[q],0.0+0.0im)
        fill!(D2s[q],0.0+0.0im)
        k=ks[q]
        αL1s[q]=-k*INV_TWO_PI
        αL2s[q]=0.5im*k
        kcs[q]=k*INV_TWO_PI
        for i in 1:N
            Ds[q][i,i]=ComplexF64(wi[i]*kappa[i],0.0)
        end
    end
    h0_tls=ws.bessel_ws.h0_tls
    h1_tls=ws.bessel_ws.h1_tls
    j0_tls=ws.bessel_ws.j0_tls
    j1_tls=ws.bessel_ws.j1_tls
    @use_threads multithreading=multithreaded for j in 2:N
        tid=Threads.threadid()
        h0vals=h0_tls[tid]
        h1vals=h1_tls[tid]
        j0vals=j0_tls[tid]
        j1vals=j1_tls[tid]
        @inbounds for i in 1:j-1
            r=R[i,j]
            _h0_h1_j0_j1_at_pidx_t!(h0vals,h1vals,j0vals,j1vals,pidx[i,j],tloc[i,j],pidxj[i,j],tlocj[i,j],r,plans0,plans1,plansj0,plansj1)
            rinv=invR[i,j]
            lt=logterm[i,j]
            Rij=Rkress[i,j]
            innij=inner[i,j]
            innji=inner[j,i]
            cDij=innij*rinv
            cDji=innji*rinv
            c1ij=Rij*cDij
            c1ji=Rij*cDji
            c2ij=wi[j]*cDij
            c2ji=wi[i]*cDji
            c3ij=wi[j]*lt*cDij
            c3ji=wi[i]*lt*cDji
            d1j0ij=-innij
            d1j0ji=-innji
            d1h0ij=im*pi*innij
            d1h0ji=im*pi*innji
            d2prefij=INV_TWO_PI*innij
            d2prefji=INV_TWO_PI*innji
            for q in 1:Mk
                k=ks[q]
                h0=h0vals[q]
                h1=h1vals[q]
                j0=j0vals[q]
                j1=j1vals[q]
                l1=αL1s[q]*j1
                Ds[q][i,j]=c1ij*l1+c2ij*αL2s[q]*h1-c3ij*l1
                Ds[q][j,i]=c1ji*l1+c2ji*αL2s[q]*h1-c3ji*l1
                D1s[q][i,j]=kcs[q]*(Rij*d1j0ij*j0+wi[j]*(lt*innij*j0+d1h0ij*h0))
                D1s[q][j,i]=kcs[q]*(Rij*d1j0ji*j0+wi[i]*(lt*innji*j0+d1h0ji*h0))
                u=j0-k*r*j1
                v=h0-k*r*h1
                D2s[q][i,j]=d2prefij*(Rij*(-u)+wi[j]*(lt*u+im*pi*v))
                D2s[q][j,i]=d2prefji*(Rij*(-u)+wi[i]*(lt*u+im*pi*v))
            end
        end
    end
    return nothing
end

function _construct_dlp_kress_matrices_derivatives_chebyshev!(Ds::Vector{Matrix{ComplexF64}},D1s::Vector{Matrix{ComplexF64}},D2s::Vector{Matrix{ComplexF64}},pts::BoundaryPoints{T},rws::DLPKressReducedH0H1J0J1ChebWorkspace{T};multithreaded::Bool=true) where {T<:Real}
    fullws=rws.fullcheb
    blk=fullws.block_cache.block
    orbits=rws.direct.orbits
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    ng=orbit_size(orbits)
    Mk=fullws.Mk
    ks=fullws.ks
    @inbounds for q in 1:Mk
        @assert size(Ds[q])==(m,m)
        @assert size(D1s[q])==(m,m)
        @assert size(D2s[q])==(m,m)
        fill!(Ds[q],0.0+0.0im)
        fill!(D1s[q],0.0+0.0im)
        fill!(D2s[q],0.0+0.0im)
    end
    αL1s=Vector{ComplexF64}(undef,Mk)
    αL2s=Vector{ComplexF64}(undef,Mk)
    @inbounds for q in 1:Mk
        αL1s[q]=-ks[q]*INV_TWO_PI
        αL2s[q]=0.5im*ks[q]
    end
    ntls=length(fullws.bessel_ws.h0_tls)
    acc_tls=[Vector{ComplexF64}(undef,Mk) for _ in 1:ntls]
    acc1_tls=[Vector{ComplexF64}(undef,Mk) for _ in 1:ntls]
    acc2_tls=[Vector{ComplexF64}(undef,Mk) for _ in 1:ntls]
    h0_tls=fullws.bessel_ws.h0_tls
    h1_tls=fullws.bessel_ws.h1_tls
    j0_tls=fullws.bessel_ws.j0_tls
    j1_tls=fullws.bessel_ws.j1_tls
    plans0=fullws.plans0
    plans1=fullws.plans1
    plansj0=fullws.plansj0
    plansj1=fullws.plansj1
    R=blk.R
    invR=blk.invR
    logterm=blk.logterm
    inner=blk.inner
    Rkress=blk.Rkress
    wi=blk.wi
    pidx=blk.pidx
    tloc=blk.tloc
    pidxj=blk.pidxj
    tlocj=blk.tlocj
    kappa=blk.kappa
    @use_threads multithreading=multithreaded for b in 1:m
        tid=Threads.threadid()
        h0vals=h0_tls[tid]
        h1vals=h1_tls[tid]
        j0vals=j0_tls[tid]
        j1vals=j1_tls[tid]
        acc=acc_tls[tid]
        acc1=acc1_tls[tid]
        acc2=acc2_tls[tid]
        @inbounds for a in 1:m
            fill!(acc,0.0+0.0im)
            fill!(acc1,0.0+0.0im)
            fill!(acc2,0.0+0.0im)
            i=Ifund[a]
            for l in 1:ng
                j=orbits.fund_to_full[l,b]
                χ=ComplexF64(orbits.fund_to_scale[l,b])
                if i==j
                    d0=ComplexF64(wi[i]*kappa[i],0.0)
                    for q in 1:Mk
                        acc[q]+=χ*d0
                    end
                else
                    r=R[i,j]
                    rinv=invR[i,j]
                    lt=logterm[i,j]
                    inn=inner[i,j]
                    Rij=Rkress[i,j]
                    wj=wi[j]
                    _h0_h1_j0_j1_at_pidx_t!(h0vals,h1vals,j0vals,j1vals,pidx[i,j],tloc[i,j],pidxj[i,j],tlocj[i,j],r,plans0,plans1,plansj0,plansj1)
                    cD1=Rij*inn*rinv
                    cD2=wj*inn*rinv
                    cD3=cD2*lt
                    cR=Rij*inn*INV_TWO_PI
                    cW=wj*inn*INV_TWO_PI
                    for q in 1:Mk
                        k=ks[q]
                        h0=h0vals[q]
                        h1=h1vals[q]
                        j0=j0vals[q]
                        j1=j1vals[q]
                        kr=k*r
                        l1=αL1s[q]*j1
                        d0=cD1*l1+cD2*αL2s[q]*h1-cD3*l1
                        d1=cR*(-k*j0)+cW*(k*(lt*j0+im*pi*h0))
                        d2=cR*(kr*j1-j0)+cW*(lt*(j0-kr*j1)+im*pi*(h0-kr*h1))
                        acc[q]+=χ*d0
                        acc1[q]+=χ*d1
                        acc2[q]+=χ*d2
                    end
                end
            end
            for q in 1:Mk
                Ds[q][a,b]=acc[q]
                D1s[q][a,b]=acc1[q]
                D2s[q][a,b]=acc2[q]
            end
        end
    end
    return nothing
end

####################################
######## FREDHOLM ASSEMBLY #########
####################################

"""
    construct_dlp_kress_matrices_chebyshev!(Fs::Vector{<:AbstractMatrix{ComplexF64}},pts::BoundaryPoints{T},ws::Union{DLPKressH1J1ChebWorkspace{T},DLPKressReducedH1J1ChebWorkspace{T}};multithreaded::Bool=true) where {T<:Real} → Nothing

Assemble the full or symmetry-reduced DLP-Kress Fredholm matrices

    F(k)=I-D(k)

for all wavenumbers stored in `ws`.
"""
function construct_dlp_kress_matrices_chebyshev!(Fs::Vector{<:AbstractMatrix{ComplexF64}},pts::BoundaryPoints{T},ws::Union{DLPKressH1J1ChebWorkspace{T},DLPKressReducedH1J1ChebWorkspace{T}};multithreaded::Bool=true) where {T<:Real}
    _construct_dlp_kress_matrices_chebyshev!(Fs,pts,ws;multithreaded=multithreaded)
    @inbounds for m in eachindex(Fs),j in axes(Fs[m],2),i in axes(Fs[m],1)
        Fs[m][i,j]*=-1
    end
    @inbounds for m in eachindex(Fs),i in axes(Fs[m],1)
        Fs[m][i,i]+=1.0+0im
    end
    return nothing
end

"""
    construct_dlp_kress_matrices_derivatives_chebyshev!(Fs::Vector{<:AbstractMatrix{ComplexF64}},F1s::Vector{<:AbstractMatrix{ComplexF64}},F2s::Vector{<:AbstractMatrix{ComplexF64}},pts::BoundaryPoints{T},ws::Union{DLPKressH0H1J0J1ChebWorkspace{T},DLPKressReducedH0H1J0J1ChebWorkspace{T}};multithreaded::Bool=true) where {T<:Real} → Nothing

Assemble

    F(k)=I-D(k),
    F'(k)=-D'(k),
    F''(k)=-D''(k)

for all wavenumbers stored in `ws`.
"""
function construct_dlp_kress_matrices_derivatives_chebyshev!(Fs::Vector{<:AbstractMatrix{ComplexF64}},F1s::Vector{<:AbstractMatrix{ComplexF64}},F2s::Vector{<:AbstractMatrix{ComplexF64}},pts::BoundaryPoints{T},ws::Union{DLPKressH0H1J0J1ChebWorkspace{T},DLPKressReducedH0H1J0J1ChebWorkspace{T}};multithreaded::Bool=true) where {T<:Real}
    @assert length(Fs)==length(F1s)==length(F2s)==_cheb_workspace_length(ws)
    _construct_dlp_kress_matrices_derivatives_chebyshev!(Fs,F1s,F2s,pts,ws;multithreaded=multithreaded)
    @inbounds for m in eachindex(Fs),j in axes(Fs[m],2),i in axes(Fs[m],1)
        Fs[m][i,j]*=-1
        F1s[m][i,j]*=-1
        F2s[m][i,j]*=-1
    end
    @inbounds for m in eachindex(Fs),i in axes(Fs[m],1)
        Fs[m][i,i]+=1.0+0im
    end
    return nothing
end

####################################
###### CHEBYSHEV BACKEND API #######
####################################

"""
    construct_matrices_chebyshev!(Tbufs::Vector{Matrix{ComplexF64}},::Val{:dlp_kress},solver::Union{DLP_kress,DLP_kress_global_corners},pts::BoundaryPoints{T},zj::AbstractVector{ComplexF64};multithreaded::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false) where {T<:Real} → Nothing

Construct the DLP-Kress Fredholm matrices

    F(k)=I-D(k)

for all complex wavenumbers in `zj` using Chebyshev interpolation of
`H₁⁽¹⁾(kr)` and `J₁(kr)`.

The output dimension is the full boundary dimension without symmetry and the
fundamental-orbit dimension when symmetry reduction is active.
"""
function construct_matrices_chebyshev!(Tbufs::Vector{Matrix{ComplexF64}},::Val{:dlp_kress},solver::Union{DLP_kress,DLP_kress_global_corners},pts::BoundaryPoints{T},zj::AbstractVector{ComplexF64};multithreaded::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false) where {T<:Real}
    @assert length(Tbufs)==length(zj)
    @blas_1 begin
        @benchit timeit=timeit "DLP-Kress workspace" directws=build_dlp_kress_workspace(solver,pts)
        @benchit timeit=timeit "DLP-Kress H1/J1 plans" chebws=build_dlp_kress_h1_j1_cheb_workspace(solver,pts,directws,ComplexF64.(zj);npanels_h=n_panels_h,npanels_j=n_panels_j,M_h=M_h,M_j=M_j,ntls=Threads.nthreads())
        n=_cheb_workspace_dim(chebws)
        @inbounds for q in eachindex(Tbufs)
            @assert size(Tbufs[q])==(n,n) "Tbufs[$q] has size $(size(Tbufs[q])), expected ($n,$n)"
            fill!(Tbufs[q],0.0+0.0im)
        end
        @benchit timeit=timeit "DLP-Kress Chebyshev" construct_dlp_kress_matrices_chebyshev!(Tbufs,pts,chebws;multithreaded=multithreaded)
    end
    return nothing
end

"""
    construct_matrices_chebyshev_with_derivatives!(Tbufs::Vector{Matrix{ComplexF64}},dTbufs::Vector{Matrix{ComplexF64}},ddTbufs::Vector{Matrix{ComplexF64}},::Val{:dlp_kress},solver::Union{DLP_kress,DLP_kress_global_corners},pts::BoundaryPoints{T},zj::AbstractVector{ComplexF64};multithreaded::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false) where {T<:Real} → Nothing

Construct the DLP-Kress Fredholm matrices and their first two wavenumber
derivatives,

    F(k)=I-D(k),
    F'(k)=-D'(k),
    F''(k)=-D''(k),

for all complex wavenumbers in `zj`.

Workspace construction and full/reduced dispatch are delegated to the common
derivative Chebyshev workspace API.
"""
function construct_matrices_chebyshev_with_derivatives!(Tbufs::Vector{Matrix{ComplexF64}},dTbufs::Vector{Matrix{ComplexF64}},ddTbufs::Vector{Matrix{ComplexF64}},::Val{:dlp_kress},solver::Union{DLP_kress,DLP_kress_global_corners},pts::BoundaryPoints{T},zj::AbstractVector{ComplexF64};multithreaded::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false) where {T<:Real}
    @assert length(Tbufs)==length(zj)
    @assert length(dTbufs)==length(zj)
    @assert length(ddTbufs)==length(zj)
    ws=build_derivative_chebyshev_workspace(solver,pts,zj;n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,timeit=timeit)
    n=_cheb_workspace_dim(ws)
    @inbounds for q in eachindex(Tbufs)
        @assert size(Tbufs[q])==(n,n) "Tbufs[$q] has size $(size(Tbufs[q])), expected ($n,$n)"
        @assert size(dTbufs[q])==(n,n) "dTbufs[$q] has size $(size(dTbufs[q])), expected ($n,$n)"
        @assert size(ddTbufs[q])==(n,n) "ddTbufs[$q] has size $(size(ddTbufs[q])), expected ($n,$n)"
    end
    @blas_1 @benchit timeit=timeit "DLP-Kress derivative Chebyshev" construct_matrices_chebyshev_with_derivatives!(Tbufs,dTbufs,ddTbufs,solver,pts,ws;multithreaded=multithreaded)
    return nothing
end

########################################
########### SOLVE-VECT BATCH ###########
########################################

"""
    adjoint_fredholm_matrix_from_dlp_chebyshev!(A::AbstractMatrix{ComplexF64},D::AbstractMatrix{ComplexF64},pts::BoundaryPoints{T},ws::DLPKressWorkspace{T}) where {T<:Real} → A

Construct the weighted formal-transpose Fredholm matrix from a full raw
Chebyshev-assembled DLP matrix,

    A=I-W⁻¹DᵀW,

with `W=diag(pts.ds)`.

No complex conjugation is applied.
"""
function adjoint_fredholm_matrix_from_dlp_chebyshev!(A::AbstractMatrix{ComplexF64},D::AbstractMatrix{ComplexF64},pts::BoundaryPoints{T},ws::DLPKressWorkspace{T}) where {T<:Real}
    N=ws.N
    ds=pts.ds
    @assert size(A)==(N,N)
    @assert size(D)==(N,N)
    fill!(A,0.0+0.0im)
    @inbounds for j in 1:N,i in 1:N
        A[i,j]=-D[j,i]*ds[j]/ds[i]
    end
    @inbounds for i in 1:N
        A[i,i]+=1.0+0.0im
    end
    return A
end

"""
    adjoint_fredholm_matrix_from_dlp_chebyshev!(A::AbstractMatrix{ComplexF64},D::AbstractMatrix{ComplexF64},pts::BoundaryPoints{T},ws::DLPKressReducedWorkspace{T}) where {T<:Real} → A

Construct the weighted formal transpose of the symmetry-reduced DLP-Kress
Fredholm operator.

For fundamental indices `a,b`, with

    i=Ifund[a],
    j=Ifund[b],

the discrete matrix is

    A[a,b]=-D[b,a]ds[j]/ds[i]+δ_ab.

No complex conjugation is applied.
"""
function adjoint_fredholm_matrix_from_dlp_chebyshev!(A::AbstractMatrix{ComplexF64},D::AbstractMatrix{ComplexF64},pts::BoundaryPoints{T},ws::DLPKressReducedWorkspace{T}) where {T<:Real}
    orbits=ws.orbits
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    ds=pts.ds
    @assert size(A)==(m,m)
    @assert size(D)==(m,m)
    fill!(A,0.0+0.0im)
    @inbounds for b in 1:m,a in 1:m
        i=Ifund[a]
        j=Ifund[b]
        A[a,b]=-D[b,a]*ds[j]/ds[i]
    end
    @inbounds for a in 1:m
        A[a,a]+=1.0+0.0im
    end
    return A
end

"""
    solve_vect(solver::Union{DLP_kress,DLP_kress_global_corners},billiard::Bi,basis::Ba,ks::Vector{T};batch_size::Int=40,multithreaded::Bool=true,use_chebyshev::Bool=true,cheb_tol::Real=1e-12,npanels_h_init::Int=15000,M_h_init::Int=5,npanels_j_init::Int=3000,M_j_init::Int=5,sampling_points::Int=50000,max_iter::Int=20,grow_panels::Real=1.5,grow_M::Int=2,cheb_verbose::Bool=false,tol=1e-12,maxiter::Int=2000,krylovdim::Int=40) where {T<:Real,Ba<:AbsBasis,Bi<:BilliardGeometry.AbsBilliard} → Tuple{Vector{Vector{ComplexF64}},Vector{BoundaryPoints{T}}}

Compute near-null DLP-Kress boundary vectors for multiple real wavenumbers.

States are processed in batches sharing one boundary discretization. When
`use_chebyshev=true`, the special-function interpolation parameters are first
validated for the complete batch and the raw DLP matrices are assembled
simultaneously.

The weighted formal-transpose Fredholm matrix is then formed from each raw DLP
matrix and its near-null vector is extracted with
[`smallest_nullvec_krylov!`](@ref).
"""
function solve_vect(solver::Union{DLP_kress,DLP_kress_global_corners},billiard::Bi,basis::Ba,ks::Vector{T};batch_size::Int=40,multithreaded::Bool=true,use_chebyshev::Bool=true,cheb_tol::Real=1e-12,npanels_h_init::Int=15_000,M_h_init::Int=5,npanels_j_init::Int=3_000,M_j_init::Int=5,sampling_points::Int=50_000,max_iter::Int=20,grow_panels::Real=1.5,grow_M::Int=2,cheb_verbose::Bool=false,tol=1e-12,maxiter::Int=2000,krylovdim::Int=40) where {T<:Real,Ba<:AbsBasis,Bi<:BilliardGeometry.AbsBilliard}
    Nk=length(ks)
    us_all=Vector{Vector{ComplexF64}}(undef,Nk)
    pts_all=Vector{BoundaryPoints{T}}(undef,Nk)
    nb=_nbatches(Nk,batch_size)
    @showprogress "solve_vect DLP Kress" for ibatch in 1:nb
        i1=_batch_first(ibatch,batch_size)
        i2=_batch_last(ibatch,batch_size,Nk)
        inds=i1:i2
        kbatch=@view ks[inds]
        pts=evaluate_points(solver,billiard,maximum(kbatch))
        if use_chebyshev
            zj=ComplexF64.(kbatch)
            nh,Mh,nj,Mj,plans0,plans1,plansj0,plansj1,errH0,errH1,errJ0,errJ1=chebyshev_params(solver,pts,zj;npanels_h_init=npanels_h_init,M_h_init=M_h_init,npanels_j_init=npanels_j_init,M_j_init=M_j_init,tol=cheb_tol,sampling_points=sampling_points,max_iter=max_iter,grow_panels=grow_panels,grow_M=grow_M,verbose=cheb_verbose)
            directws=build_dlp_kress_workspace(solver,pts)
            chebws=build_dlp_kress_h1_j1_cheb_workspace(solver,pts,directws,zj;npanels_h=nh,M_h=Mh,npanels_j=nj,M_j=Mj,ntls=Threads.nthreads())
            n=_cheb_workspace_dim(chebws)
            Ds=[Matrix{ComplexF64}(undef,n,n) for _ in eachindex(zj)]
            A=Matrix{ComplexF64}(undef,n,n)
            _construct_dlp_kress_matrices_chebyshev!(Ds,pts,chebws;multithreaded=multithreaded)
            ws_adj=chebws.direct
            for (jlocal,jglobal) in enumerate(inds)
                adjoint_fredholm_matrix_from_dlp_chebyshev!(A,Ds[jlocal],pts,ws_adj)
                _,u,_=smallest_nullvec_krylov!(A;nev=1,tol=tol,maxiter=maxiter,krylovdim=krylovdim)
                us_all[jglobal]=ComplexF64.(u)
                pts_all[jglobal]=pts
            end
        else
            ws=build_dlp_kress_workspace(solver,pts)
            n=_workspace_dim(ws)
            A=Matrix{Complex{T}}(undef,n,n)
            D=similar(A)
            for jglobal in inds
                adjoint_fredholm_matrix!(A,D,solver,pts,ws,ks[jglobal];multithreaded=multithreaded)
                _,u,_=smallest_nullvec_krylov!(A;nev=1,tol=tol,maxiter=maxiter,krylovdim=krylovdim)
                us_all[jglobal]=ComplexF64.(u)
                pts_all[jglobal]=pts
            end
        end
    end
    return us_all,pts_all
end