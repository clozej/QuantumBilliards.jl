################################################################################
# CFIE-Kress Chebyshev matrix assembly
#
# For G_k(x,y)=(i/4)H₀⁽¹⁾(kr), r=|x-y|, define the doubled DLP and SLP kernels
#
#     D(k;x,y)=(i k/2)H₁⁽¹⁾(kr)c(x,y),
#     S(k;x,y)=(i/2)H₀⁽¹⁾(kr),
#
# where
#
#     c(x,y)=((x-y)⋅n_y)/r.
#
# The CFIE Fredholm operator is
#
#     A(k)=I-[D(k)+i k S(k)],
#
# with derivatives
#
#     A'(k) =-[D'(k)+iS(k)+ikS'(k)],
#     A''(k)=-[D''(k)+2iS'(k)+ikS''(k)].
#
# Same-component interactions use the Kress logarithmic split. For boundary
# nodes x_i,x_j define
#
#     r_ij=|x_i-x_j|,   inner_ij=(x_i-x_j)⋅n_j,
#     L_ij=logterm_ij,   R_ij=Rkress_ij,   A_ij=R_ij-w_j L_ij.
#
# The Kress-corrected DLP derivatives are
#
#     D'_ij(k)=inner_ij k/(2π)
#              [-A_ij J₀(k r_ij)+iπw_j H₀⁽¹⁾(k r_ij)],
#
#     D''_ij(k)=inner_ij/(2π)
#               [-A_ij(J₀-k r_ij J₁)
#                +iπw_j(H₀⁽¹⁾-k r_ij H₁⁽¹⁾)].
#
# On the same-component diagonal,
#
#     D_ii(k)=w_i κ_i,   D'_ii(k)=D''_ii(k)=0,
#
# and
#
#     S_ii(k)=R_ii[-s_i/(2π)]
#             +w_i s_i[i/2-γ/π-(1/(2π))log((k²/4)s_i²)].
#
# Therefore
#
#     S'_ii(k) =-w_i s_i/(π k),
#     S''_ii(k)= w_i s_i/(π k²).
#
# Off-component interactions are smooth and are evaluated directly from the
# Hankel kernels.
#
# SYMMETRY
#
# When symmetry is active the complete physical boundary remains discretized.
# A SymmetryOrbitMap folds the integral-operator source sum onto the fundamental
# boundary:
#
#     B_red[a,b]=Σ_l χ_l B_full(i_a,j_l),
#
# where
#
#     B=D+i k S.
#
# The Fredholm identity is not part of the image sum:
#
#     A_red=I-B_red.
################################################################################
#
# Reference:
#   R. Kress, "Boundary Integral Equations in Time-Harmonic Acoustic
#   Scattering," Mathl. Comput. Modelling 15(3-5), 229-243 (1991).
################################################################################

const _TWO_PI=2*pi
const _INV_TWO_PI=inv(_TWO_PI)
const _EULER_OVER_PI=MathConstants.eulergamma/pi

################################################################################
############################ GEOMETRY CACHES ###################################
################################################################################

"""
    CFIEKressBlockCache{T} where {T<:Real}

Geometry and interpolation cache for one ordered target/source CFIE-Kress
component pair.
"""
struct CFIEKressBlockCache{T<:Real}
    same::Bool
    row_offset::Int
    col_offset::Int
    Ni::Int
    Nj::Int
    R::Matrix{T}
    invR::Matrix{T}
    inner::Matrix{T}
    speed_i::Vector{T}
    speed_j::Vector{T}
    wi::Vector{T}
    wj::Vector{T}
    pidx::Matrix{Int32}
    tloc::Matrix{Float64}
    pidxj::Matrix{Int32}
    tlocj::Matrix{Float64}
    logterm::Union{Nothing,Matrix{T}}
    kappa_i::Union{Nothing,Vector{T}}
    Rkress::Union{Nothing,Matrix{T}}
end

"""
    CFIEKressSystemCache{T} where {T<:Real}

Complete wavenumber-independent full-boundary CFIE-Kress geometry and
interpolation-index cache.
"""
struct CFIEKressSystemCache{T<:Real}
    blocks::Matrix{CFIEKressBlockCache{T}}
    offsets::Vector{Int}
    rmin::Float64
    rmax::Float64
end

"""
    CFIEKressReducedWorkspace{T,S} where {T<:Real}

Symmetry-reduced CFIE-Kress geometry workspace.

`S` is unrestricted because the symmetry interface may use either one symmetry
descriptor or a collection of rotation descriptors.
"""
struct CFIEKressReducedWorkspace{T<:Real,S}
    block_cache::CFIEKressSystemCache{T}
    sym::S
    orbits::SymmetryOrbitMap{T}
    global_to_block::Vector{Int}
    global_to_local::Vector{Int}
end

@inline _cfie_cheb_dim(cache::CFIEKressSystemCache)=cache.offsets[end]-1
@inline _cfie_cheb_dim(ws::CFIEKressReducedWorkspace)=fundamental_size(ws.orbits)

################################################################################
######################## SPECIAL-FUNCTION WORKSPACE ############################
################################################################################

"""
    CFIEKressH0H1J0J1BesselWorkspace

Thread-local buffers for simultaneous evaluation of `H₀⁽¹⁾`, `H₁⁽¹⁾`, `J₀`,
and `J₁`.
"""
struct CFIEKressH0H1J0J1BesselWorkspace
    h0_tls::Vector{Vector{ComplexF64}}
    h1_tls::Vector{Vector{ComplexF64}}
    j0_tls::Vector{Vector{ComplexF64}}
    j1_tls::Vector{Vector{ComplexF64}}
end

"""
    CFIEKressH0H1J0J1BesselWorkspace(Mk::Int;ntls::Int=Threads.nthreads()) → CFIEKressH0H1J0J1BesselWorkspace

Allocate thread-local CFIE-Kress interpolation buffers.
"""
function CFIEKressH0H1J0J1BesselWorkspace(Mk::Int;ntls::Int=Threads.nthreads())
    h0_tls=[Vector{ComplexF64}(undef,Mk) for _ in 1:ntls]
    h1_tls=[Vector{ComplexF64}(undef,Mk) for _ in 1:ntls]
    j0_tls=[Vector{ComplexF64}(undef,Mk) for _ in 1:ntls]
    j1_tls=[Vector{ComplexF64}(undef,Mk) for _ in 1:ntls]
    return CFIEKressH0H1J0J1BesselWorkspace(h0_tls,h1_tls,j0_tls,j1_tls)
end

################################################################################
############################ CHEBYSHEV PLANS ###################################
################################################################################

"""
    build_cfie_kress_plans(ks::AbstractVector{<:Number},rmin::Float64,rmax::Float64;npanels_h::Int=10000,M_h::Int=5,npanels_j::Int=3000,M_j::Int=5) → Tuple{Vector{ChebHankelPlanH},Vector{ChebHankelPlanH},Vector{ChebJPlan},Vector{ChebJPlan}}

Build `H₀⁽¹⁾`, `H₁⁽¹⁾`, `J₀`, and `J₁` Chebyshev plans.
"""
function build_cfie_kress_plans(ks::AbstractVector{<:Number},rmin::Float64,rmax::Float64;npanels_h::Int=10000,M_h::Int=5,npanels_j::Int=3000,M_j::Int=5)
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

################################################################################
############################ BLOCK-CACHE BUILDERS ##############################
################################################################################

"""
    build_cfie_kress_block_caches(solver::CFIE,comps::Vector{BoundaryPoints{T}};npanels_h::Int=10000,M_h::Int=5,npanels_j::Int=3000,M_j::Int=5,pad::Tuple{T,T}=(T(0.95),T(1.05)),rmin_cheb::Union{Nothing,Float64}=nothing) where {T<:Real} → CFIEKressSystemCache{T}

Build the wavenumber-independent ordered component-block cache used by
CFIE-Kress Chebyshev assembly.
"""
function build_cfie_kress_block_caches(solver::CFIE,comps::Vector{BoundaryPoints{T}};npanels_h::Int=10000,M_h::Int=5,npanels_j::Int=3000,M_j::Int=5,pad::Tuple{T,T}=(T(0.95),T(1.05)),rmin_cheb::Union{Nothing,Float64}=nothing) where {T<:Real}
    nc=length(comps)
    offs=component_offsets(comps)
    Gs=[boundary_geom_cache(p,_is_nontrivial_grading(p)) for p in comps]
    blocks=Matrix{CFIEKressBlockCache{T}}(undef,nc,nc)
    global_rmin=typemax(T)
    global_rmax=zero(T)
    for a in 1:nc,b in 1:nc
        pa=comps[a]
        pb=comps[b]
        Ga=Gs[a]
        Gb=Gs[b]
        Ni=length(pa.xy)
        Nj=length(pb.xy)
        same=a==b
        if same
            R=Ga.R
            invR=Ga.invR
            inner=Ga.inner
            speed_i=Ga.speed
            speed_j=Ga.speed
            wi=pa.ws
            wj=pa.ws
            rmin_blk=typemax(T)
            rmax_blk=zero(T)
            @inbounds for j in 1:Nj,i in 1:Ni
                i==j&&continue
                rij=R[i,j]
                if rij>eps(T)
                    rij<rmin_blk&&(rmin_blk=rij)
                    rij>rmax_blk&&(rmax_blk=rij)
                end
            end
            rmin_blk*=pad[1]
            rmax_blk*=pad[2]
            global_rmin=min(global_rmin,rmin_blk)
            global_rmax=max(global_rmax,rmax_blk)
            pidx=Matrix{Int32}(undef,Ni,Nj)
            tloc=Matrix{Float64}(undef,Ni,Nj)
            pidxj=Matrix{Int32}(undef,Ni,Nj)
            tlocj=Matrix{Float64}(undef,Ni,Nj)
            logterm=Ga.logterm
            kappa_i=Ga.kappa
            Rkress=zeros(T,Ni,Ni)
            kress_R!(Rkress)
            blocks[a,b]=CFIEKressBlockCache{T}(true,offs[a],offs[b],Ni,Nj,R,invR,inner,speed_i,speed_j,wi,wj,pidx,tloc,pidxj,tlocj,logterm,kappa_i,Rkress)
        else
            Xa=getindex.(pa.xy,1)
            Ya=getindex.(pa.xy,2)
            Xb=getindex.(pb.xy,1)
            Yb=getindex.(pb.xy,2)
            dXb=getindex.(pb.tangent,1)
            dYb=getindex.(pb.tangent,2)
            R=Matrix{T}(undef,Ni,Nj)
            invR=Matrix{T}(undef,Ni,Nj)
            inner=Matrix{T}(undef,Ni,Nj)
            @inbounds for j in 1:Nj,i in 1:Ni
                dx=Xa[i]-Xb[j]
                dy=Ya[i]-Yb[j]
                rij=hypot(dx,dy)
                R[i,j]=rij
                invR[i,j]=rij>eps(T) ? inv(rij) : zero(T)
                inner[i,j]=dYb[j]*dx-dXb[j]*dy
            end
            speed_i=Ga.speed
            speed_j=Gb.speed
            wi=pa.ws
            wj=pb.ws
            rmin_blk=typemax(T)
            rmax_blk=zero(T)
            @inbounds for j in 1:Nj,i in 1:Ni
                rij=R[i,j]
                if rij>eps(T)
                    rij<rmin_blk&&(rmin_blk=rij)
                    rij>rmax_blk&&(rmax_blk=rij)
                end
            end
            rmin_blk*=pad[1]
            rmax_blk*=pad[2]
            global_rmin=min(global_rmin,rmin_blk)
            global_rmax=max(global_rmax,rmax_blk)
            pidx=Matrix{Int32}(undef,Ni,Nj)
            tloc=Matrix{Float64}(undef,Ni,Nj)
            pidxj=Matrix{Int32}(undef,Ni,Nj)
            tlocj=Matrix{Float64}(undef,Ni,Nj)
            blocks[a,b]=CFIEKressBlockCache{T}(false,offs[a],offs[b],Ni,Nj,R,invR,inner,speed_i,speed_j,wi,wj,pidx,tloc,pidxj,tlocj,nothing,nothing,nothing)
        end
    end
    global_rmin_geom=Float64(global_rmin)
    global_rmax_geom=Float64(global_rmax)
    global_rmin_cheb=isnothing(rmin_cheb) ? global_rmin_geom : max(Float64(rmin_cheb),global_rmin_geom)
    pref_h=plan_h(0,1,1.0+0im,global_rmin_cheb,global_rmax_geom;npanels=npanels_h,M=M_h)
    pref_j=plan_j(0,1.0+0im,0.0,global_rmax_geom;npanels=npanels_j,M=M_j)
    pansh=pref_h.panels
    pansj=pref_j.panels
    for a in 1:nc,b in 1:nc
        blk=blocks[a,b]
        @inbounds for j in 1:blk.Nj,i in 1:blk.Ni
            if blk.same&&i==j
                blk.pidx[i,j]=Int32(1)
                blk.tloc[i,j]=0.0
                blk.pidxj[i,j]=Int32(1)
                blk.tlocj[i,j]=0.0
            else
                rij=Float64(blk.R[i,j])
                if rij<global_rmin_cheb
                    blk.pidx[i,j]=Int32(0)
                    blk.tloc[i,j]=0.0
                else
                    p=_find_panel(pref_h,rij)
                    P=pansh[p]
                    blk.pidx[i,j]=Int32(p)
                    blk.tloc[i,j]=(2rij-(P.b+P.a))/(P.b-P.a)
                end
                pj=_find_panel(pref_j,rij)
                Pj=pansj[pj]
                blk.pidxj[i,j]=Int32(pj)
                blk.tlocj[i,j]=(2rij-(Pj.b+Pj.a))/(Pj.b-Pj.a)
            end
        end
    end
    return CFIEKressSystemCache{T}(blocks,offs,global_rmin_cheb,global_rmax_geom)
end

"""
    build_cfie_kress_reduced_workspace(solver::CFIE,pts::Vector{BoundaryPoints{T}},sym::S;npanels_h::Int=10000,M_h::Int=5,npanels_j::Int=3000,M_j::Int=5,pad::Tuple{T,T}=(T(0.95),T(1.05)),rmin_cheb::Union{Nothing,Float64}=nothing) where {T<:Real,S} → CFIEKressReducedWorkspace{T,S}

Build the exact symmetry-reduced CFIE-Kress geometry workspace.
"""
function build_cfie_kress_reduced_workspace(solver::CFIE,pts::Vector{BoundaryPoints{T}},sym::S;npanels_h::Int=10000,M_h::Int=5,npanels_j::Int=3000,M_j::Int=5,pad::Tuple{T,T}=(T(0.95),T(1.05)),rmin_cheb::Union{Nothing,Float64}=nothing) where {T<:Real,S}
    block_cache=build_cfie_kress_block_caches(solver,pts;npanels_h=npanels_h,M_h=M_h,npanels_j=npanels_j,M_j=M_j,pad=pad,rmin_cheb=rmin_cheb)
    global_to_block,global_to_local=global_to_component_local(pts)
    orbits=symmetry_index_orbits(T,pts,sym)
    return CFIEKressReducedWorkspace(block_cache,sym,orbits,global_to_block,global_to_local)
end

################################################################################
########################### CHEBYSHEV WORKSPACE ################################
################################################################################

"""
    CFIEKressChebWorkspace{T,C} where {T<:Real}

Reusable full or symmetry-reduced CFIE-Kress Chebyshev workspace.

The workspace owns all geometry, exact symmetry mappings, special-function
plans, and thread-local buffers needed by both value-only and derivative-aware
matrix construction.
"""
struct CFIEKressChebWorkspace{T<:Real,C}
    cache::C
    plans0::Vector{ChebHankelPlanH}
    plans1::Vector{ChebHankelPlanH}
    plansj0::Vector{ChebJPlan}
    plansj1::Vector{ChebJPlan}
    bessel_ws::CFIEKressH0H1J0J1BesselWorkspace
    ks::Vector{ComplexF64}
    Mk::Int
end

@inline _cheb_workspace_dim(ws::CFIEKressChebWorkspace)=_cfie_cheb_dim(ws.cache)
@inline _cheb_workspace_length(ws::CFIEKressChebWorkspace)=ws.Mk

"""
    build_cfie_kress_cheb_workspace(solver::CFIE,pts::Vector{BoundaryPoints{T}},ks::AbstractVector{<:Number};n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,pad::Tuple{T,T}=(T(0.95),T(1.05)),rmin_cheb::Union{Nothing,Float64}=nothing,ntls::Int=Threads.nthreads(),timeit::Bool=false) where {T<:Real} → CFIEKressChebWorkspace

Build a reusable full or symmetry-reduced CFIE-Kress Chebyshev workspace.
"""
function build_cfie_kress_cheb_workspace(solver::CFIE,pts::Vector{BoundaryPoints{T}},ks::AbstractVector{<:Number};n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,pad::Tuple{T,T}=(T(0.95),T(1.05)),rmin_cheb::Union{Nothing,Float64}=nothing,ntls::Int=Threads.nthreads(),timeit::Bool=false) where {T<:Real}
    zks=ComplexF64.(ks)
    @benchit timeit=timeit "CFIE-Kress geometry cache" cache=isnothing(solver.symmetry) ? build_cfie_kress_block_caches(solver,pts;npanels_h=n_panels_h,M_h=M_h,npanels_j=n_panels_j,M_j=M_j,pad=pad,rmin_cheb=rmin_cheb) : build_cfie_kress_reduced_workspace(solver,pts,solver.symmetry;npanels_h=n_panels_h,M_h=M_h,npanels_j=n_panels_j,M_j=M_j,pad=pad,rmin_cheb=rmin_cheb)
    if cache isa CFIEKressSystemCache
        rmin=cache.rmin
        rmax=cache.rmax
    else
        rmin=cache.block_cache.rmin
        rmax=cache.block_cache.rmax
    end
    @benchit timeit=timeit "CFIE-Kress plans" plans0,plans1,plansj0,plansj1=build_cfie_kress_plans(zks,rmin,rmax;npanels_h=n_panels_h,M_h=M_h,npanels_j=n_panels_j,M_j=M_j)
    bessel_ws=CFIEKressH0H1J0J1BesselWorkspace(length(zks);ntls=ntls)
    return CFIEKressChebWorkspace{T,typeof(cache)}(cache,plans0,plans1,plansj0,plansj1,bessel_ws,zks,length(zks))
end

"""
    build_derivative_chebyshev_workspace(solver::CFIE,pts::Vector{BoundaryPoints{T}},ks::AbstractVector{<:Number};n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,pad::Tuple{T,T}=(T(0.95),T(1.05)),rmin_cheb::Union{Nothing,Float64}=nothing,timeit::Bool=false) where {T<:Real} → CFIEKressChebWorkspace

Build the reusable derivative-aware CFIE-Kress workspace used by common
higher-level algorithms such as EBIM.
"""
function build_derivative_chebyshev_workspace(solver::CFIE,pts::Vector{BoundaryPoints{T}},ks::AbstractVector{<:Number};n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,pad::Tuple{T,T}=(T(0.95),T(1.05)),rmin_cheb::Union{Nothing,Float64}=nothing,timeit::Bool=false) where {T<:Real}
    return build_cfie_kress_cheb_workspace(solver,pts,ks;n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,pad=pad,rmin_cheb=rmin_cheb,timeit=timeit)
end

################################################################################
############################ VALUE ASSEMBLY ####################################
################################################################################

function _all_k_nosymm_CFIE_chebyshev!(As::Vector{Matrix{ComplexF64}},pts::Vector{BoundaryPoints{T}},plans0::Vector{ChebHankelPlanH},plans1::Vector{ChebHankelPlanH},plansj0::Vector{ChebJPlan},plansj1::Vector{ChebJPlan},h0_tls::Vector{Vector{ComplexF64}},h1_tls::Vector{Vector{ComplexF64}},j0_tls::Vector{Vector{ComplexF64}},j1_tls::Vector{Vector{ComplexF64}},block_cache::CFIEKressSystemCache{T};multithreaded::Bool=true) where {T<:Real}
    Mk=length(plans0)
    ks=Vector{ComplexF64}(undef,Mk)
    αL1=Vector{ComplexF64}(undef,Mk)
    αL2=Vector{ComplexF64}(undef,Mk)
    iks=Vector{ComplexF64}(undef,Mk)
    @inbounds for m in 1:Mk
        km=ComplexF64(plans1[m].k)
        ks[m]=km
        αL1[m]=-km*_INV_TWO_PI
        αL2[m]=0.5im*km
        iks[m]=1im*km
        fill!(As[m],0.0+0.0im)
    end
    αM1=-_INV_TWO_PI
    αM2=0.5im
    blocks=block_cache.blocks
    nc=size(blocks,1)
    function same_block_col!(blk::CFIEKressBlockCache{T},j::Int,h0vals::Vector{ComplexF64},h1vals::Vector{ComplexF64},j0vals::Vector{ComplexF64},j1vals::Vector{ComplexF64}) where {T<:Real}
        ro=blk.row_offset
        co=blk.col_offset
        sj=blk.speed_j[j]
        wj=blk.wj[j]
        gj=co+j-1
        gi=ro+j-1
        κj=blk.kappa_i[j]
        rjj=blk.Rkress[j,j]
        @inbounds for m in 1:Mk
            km=ks[m]
            dval=ComplexF64(wj*κj,0.0)
            m1=αM1*sj
            m2=((0.5im-_EULER_OVER_PI)-_INV_TWO_PI*log((km^2/4)*(sj^2)))*sj
            sval=ComplexF64(rjj*m1,0.0)+wj*m2
            As[m][gi,gj]=1.0-(dval+iks[m]*sval)
        end
        @inbounds for i in (j+1):blk.Ni
            gi=ro+i-1
            r=blk.R[i,j]
            invr=blk.invR[i,j]
            lt=blk.logterm[i,j]
            Rij=blk.Rkress[i,j]
            inn_ij=blk.inner[i,j]
            inn_ji=blk.inner[j,i]
            si=blk.speed_i[i]
            wi=blk.wi[i]
            h0_h1_j0_j1_multi_ks_at_r!(h0vals,h1vals,j0vals,j1vals,plans0,plans1,plansj0,plansj1,blk.pidx[i,j],blk.tloc[i,j],blk.pidxj[i,j],blk.tlocj[i,j],Float64(r))
            cD1ij=Rij*inn_ij*invr
            cD2ij=wj*inn_ij*invr
            cD3ij=wj*lt*inn_ij*invr
            cD1ji=Rij*inn_ji*invr
            cD2ji=wi*inn_ji*invr
            cD3ji=wi*lt*inn_ji*invr
            cS1j=Rij*sj
            cS2j=wj*sj
            cS3j=wj*sj*lt
            cS1i=Rij*si
            cS2i=wi*si
            cS3i=wi*si*lt
            for m in 1:Mk
                h0=h0vals[m]
                h1=h1vals[m]
                j0=j0vals[m]
                j1=j1vals[m]
                L1=αL1[m]*j1
                L2=αL2[m]*h1
                M1=αM1*j0
                M2=αM2*h0
                dvalij=cD1ij*L1+cD2ij*L2-cD3ij*L1
                dvalji=cD1ji*L1+cD2ji*L2-cD3ji*L1
                svalij=cS1j*M1+cS2j*M2-cS3j*M1
                svalji=cS1i*M1+cS2i*M2-cS3i*M1
                As[m][gi,gj]=-(dvalij+iks[m]*svalij)
                As[m][gj,gi]=-(dvalji+iks[m]*svalji)
            end
        end
        return nothing
    end
    function off_block_col!(blk::CFIEKressBlockCache{T},j::Int,h0vals::Vector{ComplexF64},h1vals::Vector{ComplexF64}) where {T<:Real}
        ro=blk.row_offset
        co=blk.col_offset
        sj=blk.speed_j[j]
        wj=blk.wj[j]
        gj=co+j-1
        @inbounds for i in 1:blk.Ni
            gi=ro+i-1
            r=blk.R[i,j]
            invr=blk.invR[i,j]
            inn=blk.inner[i,j]
            h0_h1_multi_ks_at_r!(h0vals,h1vals,plans0,plans1,blk.pidx[i,j],blk.tloc[i,j],Float64(r))
            cD=wj*inn*invr
            cS=wj*sj
            for m in 1:Mk
                dval=cD*αL2[m]*h1vals[m]
                sval=cS*αM2*h0vals[m]
                As[m][gi,gj]=-(dval+iks[m]*sval)
            end
        end
        return nothing
    end
    for a in 1:nc
        blk=blocks[a,a]
        @use_threads multithreading=multithreaded for j in 1:blk.Nj
            tid=Threads.threadid()
            same_block_col!(blk,j,h0_tls[tid],h1_tls[tid],j0_tls[tid],j1_tls[tid])
        end
    end
    for a in 1:nc,b in 1:nc
        a==b&&continue
        blk=blocks[a,b]
        @use_threads multithreading=multithreaded for j in 1:blk.Nj
            tid=Threads.threadid()
            off_block_col!(blk,j,h0_tls[tid],h1_tls[tid])
        end
    end
    return nothing
end

function _one_k_nosymm_CFIE_chebyshev!(A::Matrix{ComplexF64},pts::Vector{BoundaryPoints{T}},plan0::ChebHankelPlanH,plan1::ChebHankelPlanH,planj0::ChebJPlan,planj1::ChebJPlan,h0_tls::Vector{Vector{ComplexF64}},h1_tls::Vector{Vector{ComplexF64}},j0_tls::Vector{Vector{ComplexF64}},j1_tls::Vector{Vector{ComplexF64}},block_cache::CFIEKressSystemCache{T};multithreaded::Bool=true) where {T<:Real}
    _all_k_nosymm_CFIE_chebyshev!([A],pts,[plan0],[plan1],[planj0],[planj1],h0_tls,h1_tls,j0_tls,j1_tls,block_cache;multithreaded=multithreaded)
    return nothing
end

function _all_k_symm_CFIE_chebyshev!(As::Vector{Matrix{ComplexF64}},plans0::Vector{ChebHankelPlanH},plans1::Vector{ChebHankelPlanH},plansj0::Vector{ChebJPlan},plansj1::Vector{ChebJPlan},h0_tls::Vector{Vector{ComplexF64}},h1_tls::Vector{Vector{ComplexF64}},j0_tls::Vector{Vector{ComplexF64}},j1_tls::Vector{Vector{ComplexF64}},ws::CFIEKressReducedWorkspace{T};multithreaded::Bool=true) where {T<:Real}
    Mk=length(plans0)
    orbits=ws.orbits
    Ifund=orbits.Ifund
    mred=fundamental_size(orbits)
    ng=orbit_size(orbits)
    ks=Vector{ComplexF64}(undef,Mk)
    αL1=Vector{ComplexF64}(undef,Mk)
    αL2=Vector{ComplexF64}(undef,Mk)
    iks=Vector{ComplexF64}(undef,Mk)
    @inbounds for q in 1:Mk
        k=ComplexF64(plans1[q].k)
        ks[q]=k
        αL1[q]=-k*_INV_TWO_PI
        αL2[q]=0.5im*k
        iks[q]=1im*k
        @assert size(As[q])==(mred,mred)
        fill!(As[q],0.0+0.0im)
    end
    αM1=-_INV_TWO_PI
    αM2=0.5im
    blocks=ws.block_cache.blocks
    g2b=ws.global_to_block
    g2l=ws.global_to_local
    ntls=length(h0_tls)
    acc_tls=[Vector{ComplexF64}(undef,Mk) for _ in 1:ntls]
    @use_threads multithreading=multithreaded for b in 1:mred
        tid=Threads.threadid()
        h0vals=h0_tls[tid]
        h1vals=h1_tls[tid]
        j0vals=j0_tls[tid]
        j1vals=j1_tls[tid]
        acc=acc_tls[tid]
        @inbounds for a in 1:mred
            fill!(acc,0.0+0.0im)
            ig=Ifund[a]
            ib=g2b[ig]
            i=g2l[ig]
            for l in 1:ng
                gj=orbits.fund_to_full[l,b]
                scale=ComplexF64(orbits.fund_to_scale[l,b])
                jb=g2b[gj]
                j=g2l[gj]
                blk=blocks[ib,jb]
                if blk.same
                    sj=blk.speed_j[j]
                    wj=blk.wj[j]
                    if i==j
                        κj=blk.kappa_i[j]
                        rjj=blk.Rkress[j,j]
                        for q in 1:Mk
                            k=ks[q]
                            dval=ComplexF64(wj*κj,0.0)
                            m1=αM1*sj
                            m2=((0.5im-_EULER_OVER_PI)-_INV_TWO_PI*log((k^2/4)*(sj^2)))*sj
                            sval=ComplexF64(rjj*m1,0.0)+wj*m2
                            acc[q]+=scale*(-(dval+iks[q]*sval))
                        end
                    else
                        r=blk.R[i,j]
                        invr=blk.invR[i,j]
                        lt=blk.logterm[i,j]
                        Rij=blk.Rkress[i,j]
                        inn=blk.inner[i,j]
                        h0_h1_j0_j1_multi_ks_at_r!(h0vals,h1vals,j0vals,j1vals,plans0,plans1,plansj0,plansj1,blk.pidx[i,j],blk.tloc[i,j],blk.pidxj[i,j],blk.tlocj[i,j],Float64(r))
                        cD1=scale*Rij*inn*invr
                        cD2=scale*wj*inn*invr
                        cD3=cD2*lt
                        cS1=scale*Rij*sj
                        cS2=scale*wj*sj
                        cS3=cS2*lt
                        for q in 1:Mk
                            h0=h0vals[q]
                            h1=h1vals[q]
                            j0=j0vals[q]
                            j1=j1vals[q]
                            L1=αL1[q]*j1
                            M1=αM1*j0
                            M2=αM2*h0
                            dval=cD1*L1+cD2*αL2[q]*h1-cD3*L1
                            sval=cS1*M1+cS2*M2-cS3*M1
                            acc[q]+=-(dval+iks[q]*sval)
                        end
                    end
                else
                    sj=blk.speed_j[j]
                    wj=blk.wj[j]
                    r=blk.R[i,j]
                    invr=blk.invR[i,j]
                    inn=blk.inner[i,j]
                    h0_h1_multi_ks_at_r!(h0vals,h1vals,plans0,plans1,blk.pidx[i,j],blk.tloc[i,j],Float64(r))
                    cD=scale*wj*inn*invr
                    cS=scale*wj*sj
                    for q in 1:Mk
                        dval=cD*αL2[q]*h1vals[q]
                        sval=cS*αM2*h0vals[q]
                        acc[q]+=-(dval+iks[q]*sval)
                    end
                end
            end
            if a==b
                for q in 1:Mk
                    acc[q]+=1.0+0.0im
                end
            end
            for q in 1:Mk
                As[q][a,b]=acc[q]
            end
        end
    end
    return nothing
end

function _one_k_symm_CFIE_chebyshev!(A::Matrix{ComplexF64},plan0::ChebHankelPlanH,plan1::ChebHankelPlanH,planj0::ChebJPlan,planj1::ChebJPlan,h0_tls::Vector{Vector{ComplexF64}},h1_tls::Vector{Vector{ComplexF64}},j0_tls::Vector{Vector{ComplexF64}},j1_tls::Vector{Vector{ComplexF64}},ws::CFIEKressReducedWorkspace{T};multithreaded::Bool=true) where {T<:Real}
    _all_k_symm_CFIE_chebyshev!([A],[plan0],[plan1],[planj0],[planj1],h0_tls,h1_tls,j0_tls,j1_tls,ws;multithreaded=multithreaded)
    return nothing
end

################################################################################
########################### DERIVATIVE ASSEMBLY ###############################
################################################################################

function _all_k_nosymm_CFIE_chebyshev_deriv!(As::Vector{Matrix{ComplexF64}},A1s::Vector{Matrix{ComplexF64}},A2s::Vector{Matrix{ComplexF64}},pts::Vector{BoundaryPoints{T}},plans0::Vector{ChebHankelPlanH},plans1::Vector{ChebHankelPlanH},plansj0::Vector{ChebJPlan},plansj1::Vector{ChebJPlan},h0_tls::Vector{Vector{ComplexF64}},h1_tls::Vector{Vector{ComplexF64}},j0_tls::Vector{Vector{ComplexF64}},j1_tls::Vector{Vector{ComplexF64}},block_cache::CFIEKressSystemCache{T};multithreaded::Bool=true) where {T<:Real}
    Mk=length(plans0)
    ks=Vector{ComplexF64}(undef,Mk)
    αL1=Vector{ComplexF64}(undef,Mk)
    αL2=Vector{ComplexF64}(undef,Mk)
    iks=Vector{ComplexF64}(undef,Mk)
    @inbounds for m in 1:Mk
        km=ComplexF64(plans1[m].k)
        ks[m]=km
        αL1[m]=-km*_INV_TWO_PI
        αL2[m]=0.5im*km
        iks[m]=1im*km
        fill!(As[m],0.0+0.0im)
        fill!(A1s[m],0.0+0.0im)
        fill!(A2s[m],0.0+0.0im)
    end
    αM1=-_INV_TWO_PI
    αM2=0.5im
    blocks=block_cache.blocks
    nc=size(blocks,1)
    function same_block_col_deriv!(blk::CFIEKressBlockCache{T},j::Int,h0vals::Vector{ComplexF64},h1vals::Vector{ComplexF64},j0vals::Vector{ComplexF64},j1vals::Vector{ComplexF64}) where {T<:Real}
        ro=blk.row_offset
        co=blk.col_offset
        sj=blk.speed_j[j]
        wj=blk.wj[j]
        gj=co+j-1
        gi=ro+j-1
        κj=blk.kappa_i[j]
        rjj=blk.Rkress[j,j]
        @inbounds for m in 1:Mk
            km=ks[m]
            dval=ComplexF64(wj*κj,0.0)
            m1=αM1*sj
            m2=((0.5im-_EULER_OVER_PI)-_INV_TWO_PI*log((km^2/4)*(sj^2)))*sj
            sval=ComplexF64(rjj*m1,0.0)+wj*m2
            sval1=-wj*sj/(pi*km)
            sval2=wj*sj/(pi*km^2)
            As[m][gi,gj]=1-(dval+iks[m]*sval)
            A1s[m][gi,gj]=-(1im*sval+iks[m]*sval1)
            A2s[m][gi,gj]=-(2im*sval1+iks[m]*sval2)
        end
        @inbounds for i in (j+1):blk.Ni
            gi=ro+i-1
            r=blk.R[i,j]
            invr=blk.invR[i,j]
            lt=blk.logterm[i,j]
            Rij=blk.Rkress[i,j]
            inn_ij=blk.inner[i,j]
            inn_ji=blk.inner[j,i]
            si=blk.speed_i[i]
            wi=blk.wi[i]
            h0_h1_j0_j1_multi_ks_at_r!(h0vals,h1vals,j0vals,j1vals,plans0,plans1,plansj0,plansj1,blk.pidx[i,j],blk.tloc[i,j],blk.pidxj[i,j],blk.tlocj[i,j],Float64(r))
            cDRij=Rij*inn_ij*invr
            cDWij=wj*inn_ij*invr
            cDLij=wj*lt*inn_ij*invr
            cDRji=Rij*inn_ji*invr
            cDWji=wi*inn_ji*invr
            cDLji=wi*lt*inn_ji*invr
            cRij=Rij*inn_ij*_INV_TWO_PI
            cWij=wj*inn_ij*_INV_TWO_PI
            cRji=Rij*inn_ji*_INV_TWO_PI
            cWji=wi*inn_ji*_INV_TWO_PI
            cSij=wj*sj
            cSji=wi*si
            cSRij=Rij*sj
            cSRji=Rij*si
            for m in 1:Mk
                km=ks[m]
                h0=h0vals[m]
                h1=h1vals[m]
                j0=j0vals[m]
                j1=j1vals[m]
                kr=km*r
                L1=αL1[m]*j1
                M1=αM1*j0
                M2=αM2*h0
                dval_ij=cDRij*L1+cDWij*αL2[m]*h1-cDLij*L1
                dval_ji=cDRji*L1+cDWji*αL2[m]*h1-cDLji*L1
                dval_ij_1=-cRij*km*j0+cWij*km*(lt*j0+1im*pi*h0)
                dval_ji_1=-cRji*km*j0+cWji*km*(lt*j0+1im*pi*h0)
                dval_ij_2=cRij*(kr*j1-j0)+cWij*(lt*(j0-kr*j1)+1im*pi*(h0-kr*h1))
                dval_ji_2=cRji*(kr*j1-j0)+cWji*(lt*(j0-kr*j1)+1im*pi*(h0-kr*h1))
                sval_ij=cSRij*M1+cSij*M2-cSij*lt*M1
                sval_ji=cSRji*M1+cSji*M2-cSji*lt*M1
                sval_ij_1=(r*sj*_INV_TWO_PI)*(Rij*j1-wj*(lt*j1+1im*pi*h1))
                sval_ji_1=(r*si*_INV_TWO_PI)*(Rij*j1-wi*(lt*j1+1im*pi*h1))
                sval_ij_2=(r*sj*_INV_TWO_PI/km)*(Rij*(kr*j0-j1)+wj*(lt*(j1-kr*j0)+1im*pi*(h1-kr*h0)))
                sval_ji_2=(r*si*_INV_TWO_PI/km)*(Rij*(kr*j0-j1)+wi*(lt*(j1-kr*j0)+1im*pi*(h1-kr*h0)))
                As[m][gi,gj]=-(dval_ij+iks[m]*sval_ij)
                A1s[m][gi,gj]=-(dval_ij_1+1im*sval_ij+iks[m]*sval_ij_1)
                A2s[m][gi,gj]=-(dval_ij_2+2im*sval_ij_1+iks[m]*sval_ij_2)
                As[m][gj,gi]=-(dval_ji+iks[m]*sval_ji)
                A1s[m][gj,gi]=-(dval_ji_1+1im*sval_ji+iks[m]*sval_ji_1)
                A2s[m][gj,gi]=-(dval_ji_2+2im*sval_ji_1+iks[m]*sval_ji_2)
            end
        end
        return nothing
    end
    function off_block_col_deriv!(blk::CFIEKressBlockCache{T},j::Int,h0vals::Vector{ComplexF64},h1vals::Vector{ComplexF64}) where {T<:Real}
        ro=blk.row_offset
        co=blk.col_offset
        sj=blk.speed_j[j]
        wj=blk.wj[j]
        gj=co+j-1
        @inbounds for i in 1:blk.Ni
            gi=ro+i-1
            r=blk.R[i,j]
            invr=blk.invR[i,j]
            inn=blk.inner[i,j]
            h0_h1_multi_ks_at_r!(h0vals,h1vals,plans0,plans1,blk.pidx[i,j],blk.tloc[i,j],Float64(r))
            cD=wj*inn*invr
            cD1=wj*inn
            cS=wj*sj
            for m in 1:Mk
                km=ks[m]
                h0=h0vals[m]
                h1=h1vals[m]
                dval=cD*αL2[m]*h1
                dval1=0.5im*cD1*km*h0
                dval2=0.5im*cD1*(h0-km*r*h1)
                sval=cS*αM2*h0
                sval1=-0.5im*cS*r*h1
                sval2=0.5im*cS*r*(h1-km*r*h0)/km
                As[m][gi,gj]=-(dval+iks[m]*sval)
                A1s[m][gi,gj]=-(dval1+1im*sval+iks[m]*sval1)
                A2s[m][gi,gj]=-(dval2+2im*sval1+iks[m]*sval2)
            end
        end
        return nothing
    end
    for a in 1:nc
        blk=blocks[a,a]
        @use_threads multithreading=multithreaded for j in 1:blk.Nj
            tid=Threads.threadid()
            same_block_col_deriv!(blk,j,h0_tls[tid],h1_tls[tid],j0_tls[tid],j1_tls[tid])
        end
    end
    for a in 1:nc,b in 1:nc
        a==b&&continue
        blk=blocks[a,b]
        @use_threads multithreading=multithreaded for j in 1:blk.Nj
            tid=Threads.threadid()
            off_block_col_deriv!(blk,j,h0_tls[tid],h1_tls[tid])
        end
    end
    return nothing
end

function _one_k_nosymm_CFIE_chebyshev_deriv!(A::Matrix{ComplexF64},A1::Matrix{ComplexF64},A2::Matrix{ComplexF64},pts::Vector{BoundaryPoints{T}},plan0::ChebHankelPlanH,plan1::ChebHankelPlanH,planj0::ChebJPlan,planj1::ChebJPlan,h0_tls::Vector{Vector{ComplexF64}},h1_tls::Vector{Vector{ComplexF64}},j0_tls::Vector{Vector{ComplexF64}},j1_tls::Vector{Vector{ComplexF64}},block_cache::CFIEKressSystemCache{T};multithreaded::Bool=true) where {T<:Real}
    _all_k_nosymm_CFIE_chebyshev_deriv!([A],[A1],[A2],pts,[plan0],[plan1],[planj0],[planj1],h0_tls,h1_tls,j0_tls,j1_tls,block_cache;multithreaded=multithreaded)
    return nothing
end

function _all_k_symm_CFIE_chebyshev_deriv!(As::Vector{Matrix{ComplexF64}},A1s::Vector{Matrix{ComplexF64}},A2s::Vector{Matrix{ComplexF64}},plans0::Vector{ChebHankelPlanH},plans1::Vector{ChebHankelPlanH},plansj0::Vector{ChebJPlan},plansj1::Vector{ChebJPlan},h0_tls::Vector{Vector{ComplexF64}},h1_tls::Vector{Vector{ComplexF64}},j0_tls::Vector{Vector{ComplexF64}},j1_tls::Vector{Vector{ComplexF64}},ws::CFIEKressReducedWorkspace{T};multithreaded::Bool=true) where {T<:Real}
    Mk=length(plans0)
    orbits=ws.orbits
    Ifund=orbits.Ifund
    mred=fundamental_size(orbits)
    ng=orbit_size(orbits)
    ks=Vector{ComplexF64}(undef,Mk)
    αL1=Vector{ComplexF64}(undef,Mk)
    αL2=Vector{ComplexF64}(undef,Mk)
    iks=Vector{ComplexF64}(undef,Mk)
    @inbounds for q in 1:Mk
        k=ComplexF64(plans1[q].k)
        ks[q]=k
        αL1[q]=-k*_INV_TWO_PI
        αL2[q]=0.5im*k
        iks[q]=1im*k
        @assert size(As[q])==(mred,mred)
        @assert size(A1s[q])==(mred,mred)
        @assert size(A2s[q])==(mred,mred)
        fill!(As[q],0.0+0.0im)
        fill!(A1s[q],0.0+0.0im)
        fill!(A2s[q],0.0+0.0im)
    end
    αM1=-_INV_TWO_PI
    αM2=0.5im
    blocks=ws.block_cache.blocks
    g2b=ws.global_to_block
    g2l=ws.global_to_local
    ntls=length(h0_tls)
    acc_tls=[Vector{ComplexF64}(undef,Mk) for _ in 1:ntls]
    acc1_tls=[Vector{ComplexF64}(undef,Mk) for _ in 1:ntls]
    acc2_tls=[Vector{ComplexF64}(undef,Mk) for _ in 1:ntls]
    @use_threads multithreading=multithreaded for b in 1:mred
        tid=Threads.threadid()
        h0vals=h0_tls[tid]
        h1vals=h1_tls[tid]
        j0vals=j0_tls[tid]
        j1vals=j1_tls[tid]
        acc=acc_tls[tid]
        acc1=acc1_tls[tid]
        acc2=acc2_tls[tid]
        @inbounds for a in 1:mred
            fill!(acc,0.0+0.0im)
            fill!(acc1,0.0+0.0im)
            fill!(acc2,0.0+0.0im)
            ig=Ifund[a]
            ib=g2b[ig]
            i=g2l[ig]
            for l in 1:ng
                gj=orbits.fund_to_full[l,b]
                scale=ComplexF64(orbits.fund_to_scale[l,b])
                jb=g2b[gj]
                j=g2l[gj]
                blk=blocks[ib,jb]
                if blk.same
                    sj=blk.speed_j[j]
                    wj=blk.wj[j]
                    if i==j
                        κj=blk.kappa_i[j]
                        rjj=blk.Rkress[j,j]
                        for q in 1:Mk
                            k=ks[q]
                            dval=ComplexF64(wj*κj,0.0)
                            m1=αM1*sj
                            m2=((0.5im-_EULER_OVER_PI)-_INV_TWO_PI*log((k^2/4)*(sj^2)))*sj
                            sval=ComplexF64(rjj*m1,0.0)+wj*m2
                            sval1=-wj*sj/(pi*k)
                            sval2=wj*sj/(pi*k^2)
                            acc[q]+=scale*(-(dval+iks[q]*sval))
                            acc1[q]+=scale*(-(1im*sval+iks[q]*sval1))
                            acc2[q]+=scale*(-(2im*sval1+iks[q]*sval2))
                        end
                    else
                        r=blk.R[i,j]
                        invr=blk.invR[i,j]
                        lt=blk.logterm[i,j]
                        Rij=blk.Rkress[i,j]
                        inn=blk.inner[i,j]
                        h0_h1_j0_j1_multi_ks_at_r!(h0vals,h1vals,j0vals,j1vals,plans0,plans1,plansj0,plansj1,blk.pidx[i,j],blk.tloc[i,j],blk.pidxj[i,j],blk.tlocj[i,j],Float64(r))
                        cDR=scale*Rij*inn*invr
                        cDW=scale*wj*inn*invr
                        cDL=scale*wj*lt*inn*invr
                        cR=scale*Rij*inn*_INV_TWO_PI
                        cW=scale*wj*inn*_INV_TWO_PI
                        cS=scale*wj*sj
                        cSR=scale*Rij*sj
                        for q in 1:Mk
                            k=ks[q]
                            h0=h0vals[q]
                            h1=h1vals[q]
                            j0=j0vals[q]
                            j1=j1vals[q]
                            kr=k*r
                            L1=αL1[q]*j1
                            M1=αM1*j0
                            M2=αM2*h0
                            dval=cDR*L1+cDW*αL2[q]*h1-cDL*L1
                            dval1=-cR*k*j0+cW*k*(lt*j0+1im*pi*h0)
                            dval2=cR*(kr*j1-j0)+cW*(lt*(j0-kr*j1)+1im*pi*(h0-kr*h1))
                            sval=cSR*M1+cS*M2-cS*lt*M1
                            sval1=(r*scale*sj*_INV_TWO_PI)*(Rij*j1-wj*(lt*j1+1im*pi*h1))
                            sval2=(r*scale*sj*_INV_TWO_PI/k)*(Rij*(kr*j0-j1)+wj*(lt*(j1-kr*j0)+1im*pi*(h1-kr*h0)))
                            acc[q]+=-(dval+iks[q]*sval)
                            acc1[q]+=-(dval1+1im*sval+iks[q]*sval1)
                            acc2[q]+=-(dval2+2im*sval1+iks[q]*sval2)
                        end
                    end
                else
                    sj=blk.speed_j[j]
                    wj=blk.wj[j]
                    r=blk.R[i,j]
                    invr=blk.invR[i,j]
                    inn=blk.inner[i,j]
                    h0_h1_multi_ks_at_r!(h0vals,h1vals,plans0,plans1,blk.pidx[i,j],blk.tloc[i,j],Float64(r))
                    cD=scale*wj*inn*invr
                    cD1=scale*wj*inn
                    cS=scale*wj*sj
                    for q in 1:Mk
                        k=ks[q]
                        h0=h0vals[q]
                        h1=h1vals[q]
                        dval=cD*αL2[q]*h1
                        dval1=0.5im*cD1*k*h0
                        dval2=0.5im*cD1*(h0-k*r*h1)
                        sval=cS*αM2*h0
                        sval1=-0.5im*cS*r*h1
                        sval2=0.5im*cS*r*(h1-k*r*h0)/k
                        acc[q]+=-(dval+iks[q]*sval)
                        acc1[q]+=-(dval1+1im*sval+iks[q]*sval1)
                        acc2[q]+=-(dval2+2im*sval1+iks[q]*sval2)
                    end
                end
            end
            if a==b
                for q in 1:Mk
                    acc[q]+=1.0+0.0im
                end
            end
            for q in 1:Mk
                As[q][a,b]=acc[q]
                A1s[q][a,b]=acc1[q]
                A2s[q][a,b]=acc2[q]
            end
        end
    end
    return nothing
end

function _one_k_symm_CFIE_chebyshev_deriv!(A::Matrix{ComplexF64},A1::Matrix{ComplexF64},A2::Matrix{ComplexF64},plan0::ChebHankelPlanH,plan1::ChebHankelPlanH,planj0::ChebJPlan,planj1::ChebJPlan,h0_tls::Vector{Vector{ComplexF64}},h1_tls::Vector{Vector{ComplexF64}},j0_tls::Vector{Vector{ComplexF64}},j1_tls::Vector{Vector{ComplexF64}},ws::CFIEKressReducedWorkspace{T};multithreaded::Bool=true) where {T<:Real}
    _all_k_symm_CFIE_chebyshev_deriv!([A],[A1],[A2],[plan0],[plan1],[planj0],[planj1],h0_tls,h1_tls,j0_tls,j1_tls,ws;multithreaded=multithreaded)
    return nothing
end

################################################################################
############################ PUBLIC ASSEMBLY ###################################
################################################################################

"""
    compute_kernel_matrices_CFIE_kress_chebyshev!(As::Vector{Matrix{ComplexF64}},pts::Vector{BoundaryPoints{T}},plans0::Vector{ChebHankelPlanH},plans1::Vector{ChebHankelPlanH},plansj0::Vector{ChebJPlan},plansj1::Vector{ChebJPlan},h0_tls::Vector{Vector{ComplexF64}},h1_tls::Vector{Vector{ComplexF64}},j0_tls::Vector{Vector{ComplexF64}},j1_tls::Vector{Vector{ComplexF64}},cache::Union{CFIEKressSystemCache{T},CFIEKressReducedWorkspace{T}};multithreaded::Bool=true) where {T<:Real} → Nothing

Assemble full or symmetry-reduced CFIE-Kress Fredholm matrices.
"""
function compute_kernel_matrices_CFIE_kress_chebyshev!(As::Vector{Matrix{ComplexF64}},pts::Vector{BoundaryPoints{T}},plans0::Vector{ChebHankelPlanH},plans1::Vector{ChebHankelPlanH},plansj0::Vector{ChebJPlan},plansj1::Vector{ChebJPlan},h0_tls::Vector{Vector{ComplexF64}},h1_tls::Vector{Vector{ComplexF64}},j0_tls::Vector{Vector{ComplexF64}},j1_tls::Vector{Vector{ComplexF64}},cache::Union{CFIEKressSystemCache{T},CFIEKressReducedWorkspace{T}};multithreaded::Bool=true) where {T<:Real}
    if cache isa CFIEKressSystemCache{T}
        _all_k_nosymm_CFIE_chebyshev!(As,pts,plans0,plans1,plansj0,plansj1,h0_tls,h1_tls,j0_tls,j1_tls,cache;multithreaded=multithreaded)
    else
        _all_k_symm_CFIE_chebyshev!(As,plans0,plans1,plansj0,plansj1,h0_tls,h1_tls,j0_tls,j1_tls,cache;multithreaded=multithreaded)
    end
    return nothing
end

"""
    compute_kernel_matrices_CFIE_kress_chebyshev!(A::Matrix{ComplexF64},pts::Vector{BoundaryPoints{T}},plan0::ChebHankelPlanH,plan1::ChebHankelPlanH,planj0::ChebJPlan,planj1::ChebJPlan,h0_tls::Vector{Vector{ComplexF64}},h1_tls::Vector{Vector{ComplexF64}},j0_tls::Vector{Vector{ComplexF64}},j1_tls::Vector{Vector{ComplexF64}},cache::Union{CFIEKressSystemCache{T},CFIEKressReducedWorkspace{T}};multithreaded::Bool=true) where {T<:Real} → Nothing

Assemble one full or symmetry-reduced CFIE-Kress Fredholm matrix.
"""
function compute_kernel_matrices_CFIE_kress_chebyshev!(A::Matrix{ComplexF64},pts::Vector{BoundaryPoints{T}},plan0::ChebHankelPlanH,plan1::ChebHankelPlanH,planj0::ChebJPlan,planj1::ChebJPlan,h0_tls::Vector{Vector{ComplexF64}},h1_tls::Vector{Vector{ComplexF64}},j0_tls::Vector{Vector{ComplexF64}},j1_tls::Vector{Vector{ComplexF64}},cache::Union{CFIEKressSystemCache{T},CFIEKressReducedWorkspace{T}};multithreaded::Bool=true) where {T<:Real}
    if cache isa CFIEKressSystemCache{T}
        _one_k_nosymm_CFIE_chebyshev!(A,pts,plan0,plan1,planj0,planj1,h0_tls,h1_tls,j0_tls,j1_tls,cache;multithreaded=multithreaded)
    else
        _one_k_symm_CFIE_chebyshev!(A,plan0,plan1,planj0,planj1,h0_tls,h1_tls,j0_tls,j1_tls,cache;multithreaded=multithreaded)
    end
    return nothing
end

"""
    compute_kernel_matrices_CFIE_kress_chebyshev!(As::Vector{Matrix{ComplexF64}},A1s::Vector{Matrix{ComplexF64}},A2s::Vector{Matrix{ComplexF64}},pts::Vector{BoundaryPoints{T}},plans0::Vector{ChebHankelPlanH},plans1::Vector{ChebHankelPlanH},plansj0::Vector{ChebJPlan},plansj1::Vector{ChebJPlan},h0_tls::Vector{Vector{ComplexF64}},h1_tls::Vector{Vector{ComplexF64}},j0_tls::Vector{Vector{ComplexF64}},j1_tls::Vector{Vector{ComplexF64}},cache::Union{CFIEKressSystemCache{T},CFIEKressReducedWorkspace{T}};multithreaded::Bool=true) where {T<:Real} → Nothing

Assemble full or symmetry-reduced CFIE-Kress `A`, `A'`, and `A''` matrices.
"""
function compute_kernel_matrices_CFIE_kress_chebyshev!(As::Vector{Matrix{ComplexF64}},A1s::Vector{Matrix{ComplexF64}},A2s::Vector{Matrix{ComplexF64}},pts::Vector{BoundaryPoints{T}},plans0::Vector{ChebHankelPlanH},plans1::Vector{ChebHankelPlanH},plansj0::Vector{ChebJPlan},plansj1::Vector{ChebJPlan},h0_tls::Vector{Vector{ComplexF64}},h1_tls::Vector{Vector{ComplexF64}},j0_tls::Vector{Vector{ComplexF64}},j1_tls::Vector{Vector{ComplexF64}},cache::Union{CFIEKressSystemCache{T},CFIEKressReducedWorkspace{T}};multithreaded::Bool=true) where {T<:Real}
    if cache isa CFIEKressSystemCache{T}
        _all_k_nosymm_CFIE_chebyshev_deriv!(As,A1s,A2s,pts,plans0,plans1,plansj0,plansj1,h0_tls,h1_tls,j0_tls,j1_tls,cache;multithreaded=multithreaded)
    else
        _all_k_symm_CFIE_chebyshev_deriv!(As,A1s,A2s,plans0,plans1,plansj0,plansj1,h0_tls,h1_tls,j0_tls,j1_tls,cache;multithreaded=multithreaded)
    end
    return nothing
end

"""
    compute_kernel_matrices_CFIE_kress_chebyshev!(A::Matrix{ComplexF64},A1::Matrix{ComplexF64},A2::Matrix{ComplexF64},pts::Vector{BoundaryPoints{T}},plan0::ChebHankelPlanH,plan1::ChebHankelPlanH,planj0::ChebJPlan,planj1::ChebJPlan,h0_tls::Vector{Vector{ComplexF64}},h1_tls::Vector{Vector{ComplexF64}},j0_tls::Vector{Vector{ComplexF64}},j1_tls::Vector{Vector{ComplexF64}},cache::Union{CFIEKressSystemCache{T},CFIEKressReducedWorkspace{T}};multithreaded::Bool=true) where {T<:Real} → Nothing

Assemble one full or symmetry-reduced CFIE-Kress `A`, `A'`, and `A''` triple.
"""
function compute_kernel_matrices_CFIE_kress_chebyshev!(A::Matrix{ComplexF64},A1::Matrix{ComplexF64},A2::Matrix{ComplexF64},pts::Vector{BoundaryPoints{T}},plan0::ChebHankelPlanH,plan1::ChebHankelPlanH,planj0::ChebJPlan,planj1::ChebJPlan,h0_tls::Vector{Vector{ComplexF64}},h1_tls::Vector{Vector{ComplexF64}},j0_tls::Vector{Vector{ComplexF64}},j1_tls::Vector{Vector{ComplexF64}},cache::Union{CFIEKressSystemCache{T},CFIEKressReducedWorkspace{T}};multithreaded::Bool=true) where {T<:Real}
    if cache isa CFIEKressSystemCache{T}
        _one_k_nosymm_CFIE_chebyshev_deriv!(A,A1,A2,pts,plan0,plan1,planj0,planj1,h0_tls,h1_tls,j0_tls,j1_tls,cache;multithreaded=multithreaded)
    else
        _one_k_symm_CFIE_chebyshev_deriv!(A,A1,A2,plan0,plan1,planj0,planj1,h0_tls,h1_tls,j0_tls,j1_tls,cache;multithreaded=multithreaded)
    end
    return nothing
end

################################################################################
###################### WORKSPACE-BASED COMMON API ##############################
################################################################################

"""
    construct_matrices_chebyshev!(As::Vector{Matrix{ComplexF64}},solver::CFIE,pts::Vector{BoundaryPoints{T}},ws::CFIEKressChebWorkspace;multithreaded::Bool=true) where {T<:Real} → Nothing

Construct all cached CFIE-Kress Fredholm matrices from a reusable workspace.
"""
function construct_matrices_chebyshev!(As::Vector{Matrix{ComplexF64}},solver::CFIE,pts::Vector{BoundaryPoints{T}},ws::CFIEKressChebWorkspace;multithreaded::Bool=true) where {T<:Real}
    Mk=_cheb_workspace_length(ws)
    n=_cheb_workspace_dim(ws)
    @assert length(As)==Mk
    @inbounds for q in 1:Mk
        @assert size(As[q])==(n,n) "As[$q] has size $(size(As[q])), expected ($n,$n)"
    end
    bfs=ws.bessel_ws
    compute_kernel_matrices_CFIE_kress_chebyshev!(As,pts,ws.plans0,ws.plans1,ws.plansj0,ws.plansj1,bfs.h0_tls,bfs.h1_tls,bfs.j0_tls,bfs.j1_tls,ws.cache;multithreaded=multithreaded)
    return nothing
end

"""
    construct_matrices_chebyshev_with_derivatives!(As::Vector{Matrix{ComplexF64}},A1s::Vector{Matrix{ComplexF64}},A2s::Vector{Matrix{ComplexF64}},solver::CFIE,pts::Vector{BoundaryPoints{T}},ws::CFIEKressChebWorkspace;multithreaded::Bool=true) where {T<:Real} → Nothing

Construct all cached CFIE-Kress Fredholm matrices and first two wavenumber
derivatives from a reusable workspace.
"""
function construct_matrices_chebyshev_with_derivatives!(As::Vector{Matrix{ComplexF64}},A1s::Vector{Matrix{ComplexF64}},A2s::Vector{Matrix{ComplexF64}},solver::CFIE,pts::Vector{BoundaryPoints{T}},ws::CFIEKressChebWorkspace;multithreaded::Bool=true) where {T<:Real}
    Mk=_cheb_workspace_length(ws)
    n=_cheb_workspace_dim(ws)
    @assert length(As)==Mk
    @assert length(A1s)==Mk
    @assert length(A2s)==Mk
    @inbounds for q in 1:Mk
        @assert size(As[q])==(n,n) "As[$q] has size $(size(As[q])), expected ($n,$n)"
        @assert size(A1s[q])==(n,n) "A1s[$q] has size $(size(A1s[q])), expected ($n,$n)"
        @assert size(A2s[q])==(n,n) "A2s[$q] has size $(size(A2s[q])), expected ($n,$n)"
    end
    bfs=ws.bessel_ws
    compute_kernel_matrices_CFIE_kress_chebyshev!(As,A1s,A2s,pts,ws.plans0,ws.plans1,ws.plansj0,ws.plansj1,bfs.h0_tls,bfs.h1_tls,bfs.j0_tls,bfs.j1_tls,ws.cache;multithreaded=multithreaded)
    return nothing
end

"""
    construct_matrix_chebyshev_at!(A::Matrix{ComplexF64},solver::CFIE,pts::Vector{BoundaryPoints{T}},ws::CFIEKressChebWorkspace,idx::Int;multithreaded::Bool=true) where {T<:Real} → Nothing

Construct the value-only CFIE-Kress Fredholm matrix corresponding to cached
wavenumber `idx`.
"""
function construct_matrix_chebyshev_at!(A::Matrix{ComplexF64},solver::CFIE,pts::Vector{BoundaryPoints{T}},ws::CFIEKressChebWorkspace,idx::Int;multithreaded::Bool=true) where {T<:Real}
    checkbounds(ws.ks,idx)
    n=_cheb_workspace_dim(ws)
    @assert size(A)==(n,n)
    bfs=ws.bessel_ws
    compute_kernel_matrices_CFIE_kress_chebyshev!(A,pts,ws.plans0[idx],ws.plans1[idx],ws.plansj0[idx],ws.plansj1[idx],bfs.h0_tls,bfs.h1_tls,bfs.j0_tls,bfs.j1_tls,ws.cache;multithreaded=multithreaded)
    return nothing
end

"""
    construct_matrix_chebyshev_with_derivatives_at!(A::Matrix{ComplexF64},A1::Matrix{ComplexF64},A2::Matrix{ComplexF64},solver::CFIE,pts::Vector{BoundaryPoints{T}},ws::CFIEKressChebWorkspace,idx::Int;multithreaded::Bool=true) where {T<:Real} → Nothing

Construct one cached CFIE-Kress Fredholm matrix and its first two wavenumber
derivatives.
"""
function construct_matrix_chebyshev_with_derivatives_at!(A::Matrix{ComplexF64},A1::Matrix{ComplexF64},A2::Matrix{ComplexF64},solver::CFIE,pts::Vector{BoundaryPoints{T}},ws::CFIEKressChebWorkspace,idx::Int;multithreaded::Bool=true) where {T<:Real}
    checkbounds(ws.ks,idx)
    n=_cheb_workspace_dim(ws)
    @assert size(A)==(n,n)
    @assert size(A1)==(n,n)
    @assert size(A2)==(n,n)
    bfs=ws.bessel_ws
    compute_kernel_matrices_CFIE_kress_chebyshev!(A,A1,A2,pts,ws.plans0[idx],ws.plans1[idx],ws.plansj0[idx],ws.plansj1[idx],bfs.h0_tls,bfs.h1_tls,bfs.j0_tls,bfs.j1_tls,ws.cache;multithreaded=multithreaded)
    return nothing
end

################################################################################
######################## CHEBYSHEV BACKEND API ################################
################################################################################

"""
    construct_matrices_chebyshev!(Tbufs::Vector{Matrix{ComplexF64}},::Val{:cfie_kress},solver::CFIE,pts::Vector{BoundaryPoints{T}},zj::AbstractVector{ComplexF64};multithreaded::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false) where {T<:Real} → Nothing

Construct the CFIE-Kress Fredholm matrices

    A(k)=I-[D(k)+ikS(k)]

for all complex wavenumbers in `zj`.
"""
function construct_matrices_chebyshev!(Tbufs::Vector{Matrix{ComplexF64}},::Val{:cfie_kress},solver::CFIE,pts::Vector{BoundaryPoints{T}},zj::AbstractVector{ComplexF64};multithreaded::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false) where {T<:Real}
    @assert length(Tbufs)==length(zj)
    ws=build_cfie_kress_cheb_workspace(solver,pts,zj;n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,ntls=Threads.nthreads(),timeit=timeit)
    n=_cheb_workspace_dim(ws)
    @inbounds for q in eachindex(Tbufs)
        @assert size(Tbufs[q])==(n,n) "Tbufs[$q] has size $(size(Tbufs[q])), expected ($n,$n)"
    end
    @blas_1 @benchit timeit=timeit "CFIE-Kress Chebyshev" construct_matrices_chebyshev!(Tbufs,solver,pts,ws;multithreaded=multithreaded)
    return nothing
end

"""
    construct_matrices_chebyshev_with_derivatives!(Tbufs::Vector{Matrix{ComplexF64}},dTbufs::Vector{Matrix{ComplexF64}},ddTbufs::Vector{Matrix{ComplexF64}},::Val{:cfie_kress},solver::CFIE,pts::Vector{BoundaryPoints{T}},zj::AbstractVector{ComplexF64};multithreaded::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false) where {T<:Real} → Nothing

Construct the CFIE-Kress Fredholm matrices and first two wavenumber derivatives

    A(k)=I-[D(k)+ikS(k)],
    A'(k)=-[D'(k)+iS(k)+ikS'(k)],
    A''(k)=-[D''(k)+2iS'(k)+ikS''(k)]

for all complex wavenumbers in `zj`.
"""
function construct_matrices_chebyshev_with_derivatives!(Tbufs::Vector{Matrix{ComplexF64}},dTbufs::Vector{Matrix{ComplexF64}},ddTbufs::Vector{Matrix{ComplexF64}},::Val{:cfie_kress},solver::CFIE,pts::Vector{BoundaryPoints{T}},zj::AbstractVector{ComplexF64};multithreaded::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false) where {T<:Real}
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
    @blas_1 @benchit timeit=timeit "CFIE-Kress derivative Chebyshev" construct_matrices_chebyshev_with_derivatives!(Tbufs,dTbufs,ddTbufs,solver,pts,ws;multithreaded=multithreaded)
    return nothing
end

################################################################################
############################ SOLVE-VECT BATCH ##################################
################################################################################

"""
    solve_vect(solver::CFIE,billiard::Bi,basis::Ba,ks::Vector{T};batch_size::Int=40,multithreaded::Bool=true,use_chebyshev::Bool=true,cheb_tol::Real=1e-12,npanels_h_init::Int=15000,M_h_init::Int=5,npanels_j_init::Int=3000,M_j_init::Int=5,sampling_points::Int=50000,max_iter::Int=20,grow_panels::Real=1.5,grow_M::Int=2,cheb_verbose::Bool=false,tol=1e-12,maxiter::Int=2000,krylovdim::Int=40) where {T<:Real,Ba<:AbsBasis,Bi<:BilliardGeometry.AbsBilliard} → Tuple{Vector{Vector{ComplexF64}},Vector{Vector{BoundaryPoints{T}}}}

Compute CFIE-Kress near-null boundary vectors for several real wavenumbers.

Each batch shares one boundary discretization. The Chebyshev pathway tunes the
special-function interpolation parameters once and then constructs a reusable
full or symmetry-reduced `CFIEKressChebWorkspace`.
"""
function solve_vect(solver::CFIE,billiard::Bi,basis::Ba,ks::Vector{T};batch_size::Int=40,multithreaded::Bool=true,use_chebyshev::Bool=true,cheb_tol::Real=1e-12,npanels_h_init::Int=15_000,M_h_init::Int=5,npanels_j_init::Int=3_000,M_j_init::Int=5,sampling_points::Int=50_000,max_iter::Int=20,grow_panels::Real=1.5,grow_M::Int=2,cheb_verbose::Bool=false,tol=1e-12,maxiter::Int=2000,krylovdim::Int=40) where {T<:Real,Ba<:AbsBasis,Bi<:BilliardGeometry.AbsBilliard}
    Nk=length(ks)
    us_all=Vector{Vector{ComplexF64}}(undef,Nk)
    pts_all=Vector{Vector{BoundaryPoints{T}}}(undef,Nk)
    nb=_nbatches(Nk,batch_size)
    @showprogress "solve_vect CFIE Kress" for ibatch in 1:nb
        i1=_batch_first(ibatch,batch_size)
        i2=_batch_last(ibatch,batch_size,Nk)
        inds=i1:i2
        kbatch=@view ks[inds]
        kmax=maximum(kbatch)
        pts=evaluate_points(solver,billiard,kmax)
        if use_chebyshev
            zj=ComplexF64.(kbatch)
            nh,Mh,nj,Mj,plans0,plans1,plansj0,plansj1,errH0,errH1,errJ0,errJ1=chebyshev_params(solver,pts,zj;npanels_h_init=npanels_h_init,M_h_init=M_h_init,npanels_j_init=npanels_j_init,M_j_init=M_j_init,tol=cheb_tol,sampling_points=sampling_points,max_iter=max_iter,grow_panels=grow_panels,grow_M=grow_M,verbose=cheb_verbose)
            ws=build_cfie_kress_cheb_workspace(solver,pts,zj;n_panels_h=nh,M_h=Mh,n_panels_j=nj,M_j=Mj,ntls=Threads.nthreads())
            Mk=_cheb_workspace_length(ws)
            Nmat=_cheb_workspace_dim(ws)
            As=[Matrix{ComplexF64}(undef,Nmat,Nmat) for _ in 1:Mk]
            construct_matrices_chebyshev!(As,solver,pts,ws;multithreaded=multithreaded)
            for (jlocal,jglobal) in enumerate(inds)
                _,u,_=smallest_nullvec_krylov!(As[jlocal];nev=1,tol=tol,maxiter=maxiter,krylovdim=krylovdim)
                us_all[jglobal]=ComplexF64.(u)
                pts_all[jglobal]=pts
            end
        else
            ws=build_cfie_kress_workspace(solver,pts)
            Nmat=_cfie_workspace_dim(ws)
            A=Matrix{Complex{T}}(undef,Nmat,Nmat)
            for jglobal in inds
                _,u=solve_vect(solver,basis,A,pts,ws,ks[jglobal];multithreaded=multithreaded,tol=tol,maxiter=maxiter,krylovdim=krylovdim)
                us_all[jglobal]=ComplexF64.(u)
                pts_all[jglobal]=pts
            end
        end
    end
    return us_all,pts_all
end