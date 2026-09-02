# Useful reading:
# - Kress, R., Boundary integral equations in time-harmonic acoustic scattering. Math. Comput. Modelling 15 (1991), 229-243.
# - Barnett, A. H. / Betcke, T., mpspack DLP implementation.
# - Zhao, L. / Barnett, A., Robust and efficient solution of the drum problem via Nyström approximation of the Fredholm determinant.

const two_pi=2*pi
const inv_two_pi=inv(two_pi)
abstract type DLP<:SweepSolver end
"""
    DLPKressWorkspace{T,M}

Workspace cache for repeated DLP-Kress matrix assembly on a fixed boundary
discretization.

## Description
This object stores all geometry-dependent data that can be reused across many
wavenumbers `k`. For repeated assembly on the same boundary nodes, the expensive
pairwise geometric quantities are therefore computed only once, while the
k-dependent special-function evaluations are updated for each new wavenumber.

The workspace contains the Kress logarithmic correction matrix, the pairwise
geometry cache, flat panel arrays and the matrix dimension.

## Attributes
* `Rmat::M`: Dense real Kress correction matrix for the logarithmic singular part.
* `G::BoundaryGeomCache{T}`: Pairwise geometry cache.
* `parr::BoundaryPanelArrays{T}`: Flat boundary-array cache.
* `N::Int`: Number of full-boundary points and matrix dimension.
"""
struct DLPKressWorkspace{T<:Real,M<:AbstractMatrix{T}}
    Rmat::M
    G::BoundaryGeomCache{T}
    parr::BoundaryPanelArrays{T}
    N::Int
end

"""
    DLPKressReducedWorkspace{T,M}

Workspace for symmetry-reduced DLP-Kress assembly.

The full periodic Kress workspace is retained because the singular same-copy
interaction must still be evaluated on the complete periodic discretization.
The symmetry information is stored in `orbits`, which maps each fundamental
boundary node to all full-boundary symmetry images and their irrep factors.

## Attributes
* `full::DLPKressWorkspace{T,M}`: Geometry and Kress data for the complete periodic boundary.
* `orbits::SymmetryOrbitMap{T}`: Exact full-to-fundamental and fundamental-to-full symmetry map.
"""
struct DLPKressReducedWorkspace{T<:Real,M<:AbstractMatrix{T}}
    full::DLPKressWorkspace{T,M}
    orbits::SymmetryOrbitMap{T}
end

const DLPKressAnyWorkspace{T}=Union{DLPKressWorkspace{T},DLPKressReducedWorkspace{T}}

"""
    DLP_kress{T,Bi,Sym} <: DLP

Solver for the smooth periodic Kress-corrected double-layer Fredholm
formulation.

## Description
`DLP_kress` implements the periodic Kress treatment of the Helmholtz
double-layer operator on a single smooth closed boundary.

The associated Fredholm matrix is

    F(k)=I-D(k),

where `D(k)` is the Kress-corrected Nyström discretization of the interior
Helmholtz double-layer operator.

The singular same-boundary interaction is split into the universal periodic
logarithmic kernel

    log(4sin²((t-s)/2))

and a smooth remainder.

## Attributes
* `sampler::Vector{BilliardGeometry.LinearNodes}`: Sampling descriptor retained for the common solver API.
* `pts_scaling_factor::Vector{T}`: Boundary-resolution scaling factors.
* `dim_scaling_factor::T`: Compatibility field used by generic refinement code.
* `eps::T`: Numerical tolerance.
* `min_dim::Int64`: Minimum dimension compatibility field.
* `min_pts::Int64`: Minimum number of boundary points.
* `billiard::Bi`: Underlying billiard geometry.
* `symmetry::Sym`: Optional reflection or rotation symmetry descriptor.
"""
struct DLP_kress{T<:Real,Bi<:BilliardGeometry.AbsBilliard,Sym}<:DLP
    sampler::Vector{BilliardGeometry.LinearNodes}
    pts_scaling_factor::Vector{T}
    dim_scaling_factor::T
    eps::T
    min_dim::Int64
    min_pts::Int64
    billiard::Bi
    symmetry::Sym
end

"""
    DLP_kress_global_corners{T,Bi,Sym} <: DLP

Solver for the globally graded Kress-corrected double-layer Fredholm
formulation on a piecewise smooth outer boundary.

## Description
This is the corner-capable counterpart of [`DLP_kress`](@ref). True corners
are handled by a global Kress grading map that clusters nodes near the corner
locations.

If the composite boundary contains no true corners, the implementation falls
back to an ungraded smooth-composite periodic discretization.

## Attributes
* `sampler::Vector{BilliardGeometry.LinearNodes}`: Sampling descriptor retained for the common solver API.
* `pts_scaling_factor::Vector{T}`: Boundary-resolution scaling factors.
* `dim_scaling_factor::T`: Compatibility field used by generic refinement code.
* `eps::T`: Numerical tolerance.
* `min_dim::Int64`: Minimum dimension compatibility field.
* `min_pts::Int64`: Minimum number of boundary points.
* `billiard::Bi`: Underlying billiard geometry.
* `symmetry::Sym`: Optional reflection or rotation symmetry descriptor.
* `kressq::Int`: Order of the Kress grading map.
* `min_t_spacing::Real`: Minimum allowed spacing after grading.
"""
struct DLP_kress_global_corners{T<:Real,Bi<:BilliardGeometry.AbsBilliard,Sym}<:DLP
    sampler::Vector{BilliardGeometry.LinearNodes}
    pts_scaling_factor::Vector{T}
    dim_scaling_factor::T
    eps::T
    min_dim::Int64
    min_pts::Int64
    billiard::Bi
    symmetry::Sym
    kressq::Int
    min_t_spacing::Real
end

function DLP_kress(pts_scaling_factor::Union{T,Vector{T}},billiard::Bi;min_pts=200,eps=T(1e-15),symmetry::Union{Nothing,BilliardGeometry.AbsSymmetry}=nothing) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    bs=pts_scaling_factor isa T ? [pts_scaling_factor] : pts_scaling_factor
    sampler=[BilliardGeometry.LinearNodes()]
    Sym=typeof(symmetry)
    return DLP_kress{T,Bi,Sym}(sampler,bs,bs[1],eps,min_pts,min_pts,billiard,symmetry)
end

function DLP_kress_global_corners(pts_scaling_factor::Union{T,Vector{T}},billiard::Bi;min_pts=200,eps=T(1e-15),symmetry::Union{Nothing,BilliardGeometry.AbsSymmetry}=nothing,kressq=2,min_t_spacing=1e-12) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    bs=pts_scaling_factor isa T ? [pts_scaling_factor] : pts_scaling_factor
    sampler=[BilliardGeometry.LinearNodes()]
    Sym=typeof(symmetry)
    return DLP_kress_global_corners{T,Bi,Sym}(sampler,bs,bs[1],eps,min_pts,min_pts,billiard,symmetry,kressq,min_t_spacing)
end

@inline function _is_nontrivial_grading(pts::BoundaryPoints{T}) where {T<:Real}
    length(pts.ws_der)==length(pts)||return false
    return maximum(abs.(pts.ws_der.-one(T)))>sqrt(eps(T))
end

@inline _is_dlp_kress_graded(::DLP,pts::BoundaryPoints)=_is_nontrivial_grading(pts)
@inline _dlp_kress_use_reduced(solver::DLP)=!isnothing(solver.symmetry)

"""
    build_Rmat_dlp_kress(
        solver::DLP,
        pts::BoundaryPoints{T},
    ) where {T<:Real} → Matrix{T}

Build the Kress matrix for a DLP-Kress boundary
discretization.

The rule is selected from the actual boundary discretization. A
nontrivially graded boundary uses [`kress_R_corner!`](@ref), while an ungraded
periodic boundary uses [`kress_R!`](@ref).

## Arguments
* `solver::DLP`: DLP-Kress solver. Retained for the common workspace API.
* `pts::BoundaryPoints{T}`: Boundary discretization.

## Returns
* `Rmat::Matrix{T}`: Dense `N×N` Kress matrix.
"""
function build_Rmat_dlp_kress(solver::DLP,pts::BoundaryPoints{T}) where {T<:Real}
    N=length(pts)
    Rmat=zeros(T,N,N)
    kress_R!(Rmat)
    return Rmat
end

@inline function _composite_arclength(comp::Vector{C},t::T) where {T<:Real,C<:BilliardGeometry.AbsCurve}
    _,_,Ltot=component_lengths(comp)
    τ=mod(t,T(two_pi))
    target=Ltot*τ/T(two_pi)
    offset=zero(T)
    @inbounds for j in eachindex(comp)
        Lj=T(comp[j].length)
        if target<offset+Lj||j==lastindex(comp)
            u=clamp((target-offset)/Lj,zero(T),one(T))
            return offset+BilliardGeometry.arc_length(comp[j],u)
        end
        offset+=Lj
    end
    return Ltot
end

# since we will need future compatibility for holes, these need to be added (for CFIE)
@inline _boundary_components(boundary::AbstractVector{<:BilliardGeometry.AbsCurve})=[boundary]
@inline _boundary_components(boundary::AbstractVector{<:AbstractVector{<:BilliardGeometry.AbsCurve}})=boundary
@inline _is_single_composite_boundary(boundary)=length(_boundary_components(boundary))==1

"""
    _evaluate_points(
        solver::DLP_kress{T},
        crv::C,
        k::T,
        idx::Int,
    ) where {T<:Real,C<:BilliardGeometry.AbsCurve} → BoundaryPoints{T}

Construct the periodic DLP-Kress discretization of a single smooth closed BilliardGeometry.curve.

The computational Kress variable is

    σ_j=2π(j-1/2)/N,

and the physical curve parameter is `u=σ/(2π)`. The node count is adjusted so
that it is compatible with both periodic Kress quadrature and the active
symmetry group.

## Arguments
* `solver::DLP_kress{T}`: Smooth DLP-Kress solver.
* `crv::C`: Smooth closed curve.
* `k::T`: Wavenumber controlling the discretization density.
* `idx::Int`: Boundary-component label.

## Returns
* `pts::BoundaryPoints{T}`: Periodic boundary discretization.
"""
function _evaluate_points(solver::DLP_kress{T},crv::C,k::T,idx::Int) where {T<:Real,C<:BilliardGeometry.AbsCurve}
    L=T(crv.length)
    N=max(solver.min_pts,round(Int,k*L*solver.pts_scaling_factor[1]/two_pi))
    needed=isnothing(solver.symmetry) ? 2 : lcm(2,symmetry_node_multiple(solver.symmetry))
    N=cld(N,needed)*needed
    ts=T[s_mid(j,N) for j in 1:N] # computational Kress parameter σ ∈ [0,2π)
    tphys=ts./T(two_pi) # physical BilliardGeometry curve parameter u ∈ [0,1)
    xy=BilliardGeometry.curve(crv,tphys)
    tangent_1st=tangent(crv,tphys)./T(two_pi) # dγ/dσ
    tangent_2nd=tangent_2(crv,tphys)./T(two_pi)^2 # d²γ/dσ²
    s=BilliardGeometry.arc_length(crv,tphys) # physical arclength position
    h=T(two_pi)/T(N)
    ds=Vector{T}(undef,N) # physical quadrature weights ds=|γ_σ|dσ
    @inbounds for i in 1:N
        v=tangent_1st[i]
        ds[i]=hypot(v[1],v[2])*h
    end
    ws=fill(h,N) # computational periodic quadrature weight dσ
    ws_der=ones(T,N) # no Kress grading
    z=SVector{2,T}(zero(T),zero(T))
    return BoundaryPoints(xy,tangent_1st,tangent_2nd,ts,tphys,ws,ws_der,s,ds,idx,true,z,z,z,z)
end

"""
    _evaluate_points_smooth_composite(
        solver::DLP_kress_global_corners{T},
        comp::Vector{C},
        k::T,
        idx::Int,
    ) where {T<:Real,C<:BilliardGeometry.AbsCurve} → BoundaryPoints{T}

Construct an ungraded periodic Nyström discretization of a smooth composite
boundary.

This is the fallback used by [`DLP_kress_global_corners`](@ref) when the
boundary consists of several joined curve pieces but contains no true corners.

## Arguments
* `solver::DLP_kress_global_corners{T}`: Global-corner DLP-Kress solver.
* `comp::Vector{C}`: Smooth curve pieces forming one closed component.
* `k::T`: Wavenumber controlling the discretization density.
* `idx::Int`: Boundary-component label.

## Returns
* `pts::BoundaryPoints{T}`: Ungraded periodic boundary discretization.
"""
function _evaluate_points_smooth_composite(solver::DLP_kress_global_corners{T},comp::Vector{C},k::T,idx::Int) where {T<:Real,C<:BilliardGeometry.AbsCurve}
    _,_,Ltot=component_lengths(comp)
    N=max(solver.min_pts,round(Int,k*Ltot*solver.pts_scaling_factor[1]/two_pi))
    needed=isnothing(solver.symmetry) ? 2 : lcm(2,symmetry_node_multiple(solver.symmetry))
    N=cld(N,needed)*needed
    ts=T[s_mid(j,N) for j in 1:N]
    tphys=copy(ts) # physical global composite parameter t ∈ [0,2π)
    xy=Vector{SVector{2,T}}(undef,N)
    tangent_1st=Vector{SVector{2,T}}(undef,N)
    tangent_2nd=Vector{SVector{2,T}}(undef,N)
    s=Vector{T}(undef,N)
    ds=Vector{T}(undef,N)
    h=T(two_pi)/T(N)
    @inbounds for i in 1:N
        q,γt,γtt=_eval_composite_geom_global_t(T,comp,tphys[i])
        xy[i]=q
        tangent_1st[i]=γt
        tangent_2nd[i]=γtt
        s[i]=_composite_arclength(comp,tphys[i])
        ds[i]=hypot(γt[1],γt[2])*h
    end
    ws=fill(h,N)
    ws_der=ones(T,N)
    z=SVector{2,T}(zero(T),zero(T))
    return BoundaryPoints(xy,tangent_1st,tangent_2nd,ts,tphys,ws,ws_der,s,ds,idx,true,z,z,z,z)
end

"""
    _evaluate_points(
        solver::DLP_kress_global_corners{T},
        comp::Vector{C},
        k::T,
        idx::Int,
    ) where {T<:Real,C<:BilliardGeometry.AbsCurve} → BoundaryPoints{T}

Construct a globally graded Nyström discretization of a piecewise smooth closed
boundary.

True corners are detected from tangent discontinuities between adjacent curve
segments. If no corners are present, the function delegates to
[`_evaluate_points_smooth_composite`](@ref).

For the grading map `t=t(σ)`,

    γ_σ=γ_t t_σ,
    γ_σσ=γ_tt(t_σ)²+γ_t t_σσ.

## Arguments
* `solver::DLP_kress_global_corners{T}`: Global-corner DLP-Kress solver.
* `comp::Vector{C}`: Curve segments forming one closed component.
* `k::T`: Wavenumber controlling the discretization density.
* `idx::Int`: Boundary-component label.

## Returns
* `pts::BoundaryPoints{T}`: Globally graded boundary discretization.
"""
function _evaluate_points(solver::DLP_kress_global_corners{T},comp::Vector{C},k::T,idx::Int) where {T<:Real,C<:BilliardGeometry.AbsCurve}
    corners=_component_corner_locations(T,comp)
    isempty(corners)&&return _evaluate_points_smooth_composite(solver,comp,k,idx)
    _,_,Ltot=component_lengths(comp)
    N=max(solver.min_pts,round(Int,k*Ltot*solver.pts_scaling_factor[1]/two_pi))
    needed=isnothing(solver.symmetry) ? 2 : symmetry_node_multiple(solver.symmetry)
    N=cld(N,needed)*needed
    σ,tmap,jac,jac2,_=multi_kress_graded_nodes_data(T,N,corners;q=solver.kressq,minsep_tol=solver.min_t_spacing)
    tphys=tmap
    xy=Vector{SVector{2,T}}(undef,N)
    tangent_1st=Vector{SVector{2,T}}(undef,N)
    tangent_2nd=Vector{SVector{2,T}}(undef,N)
    s=Vector{T}(undef,N)
    h=T(two_pi)/T(N)
    ds=Vector{T}(undef,N)
    @inbounds for i in 1:N
        q,γt,γtt=_eval_composite_geom_global_t(T,comp,tphys[i])
        tangent_1st[i]=γt*jac[i]
        tangent_2nd[i]=γtt*(jac[i]^2)+γt*jac2[i]
        xy[i]=q
        s[i]=_composite_arclength(comp,tphys[i])
        v=tangent_1st[i]
        ds[i]=hypot(v[1],v[2])*h
    end
    ws=fill(h,N)
    z=SVector{2,T}(zero(T),zero(T))
    return BoundaryPoints(xy,tangent_1st,tangent_2nd,σ,tphys,ws,jac,s,ds,idx,true,z,z,z,z)
end

"""
    evaluate_points(
        solver::DLP_kress{T},
        billiard::Bi,
        k::T,
    ) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard} → BoundaryPoints{T}

Construct the smooth periodic DLP-Kress discretization of `billiard`.

`DLP_kress` requires exactly one smooth closed outer boundary represented by a
single curve.

## Arguments
* `solver::DLP_kress{T}`: Smooth DLP-Kress solver.
* `billiard::Bi`: Billiard geometry.
* `k::T`: Wavenumber controlling the discretization density.

## Returns
* `pts::BoundaryPoints{T}`: Smooth periodic boundary discretization.
"""
function evaluate_points(solver::DLP_kress{T},billiard::Bi,k::T) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    boundary=billiard.full_boundary
    isempty(boundary)&&error("Boundary cannot be empty.")
    comps=_boundary_components(boundary)
    length(comps)==1||error("DLP_kress supports exactly one smooth outer boundary component. Geometries with holes or multiple closed components require a multiply-connected boundary-integral formulation.")
    comp=comps[1]
    length(comp)==1||error("DLP_kress requires a single smooth closed curve. This boundary is piecewise/composite; use DLP_kress_global_corners instead.")
    return _evaluate_points(solver,comp[1],k,1)
end

"""
    evaluate_points(
        solver::DLP_kress_global_corners{T},
        billiard::Bi,
        k::T,
    ) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard} → BoundaryPoints{T}

Construct the DLP-Kress discretization of a smooth or piecewise smooth
single-component boundary.

A single smooth curve is delegated to the ordinary periodic DLP-Kress
discretization. A composite outer boundary is globally graded whenever true
corners are present.

## Arguments
* `solver::DLP_kress_global_corners{T}`: Global-corner DLP-Kress solver.
* `billiard::Bi`: Billiard geometry.
* `k::T`: Wavenumber controlling the discretization density.

## Returns
* `pts::BoundaryPoints{T}`: Smooth or globally graded boundary discretization.
"""
function evaluate_points(solver::DLP_kress_global_corners{T},billiard::Bi,k::T) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    boundary=billiard.full_boundary
    isempty(boundary)&&error("Boundary cannot be empty.")
    comps=_boundary_components(boundary)
    length(comps)==1||error("DLP_kress_global_corners supports exactly one outer boundary component. Multiple closed components require a multiply-connected boundary-integral formulation.")
    comp=comps[1]
    if length(comp)==1
        base_solver=DLP_kress(solver.pts_scaling_factor,solver.billiard;min_pts=solver.min_pts,eps=solver.eps,symmetry=solver.symmetry)
        return _evaluate_points(base_solver,comp[1],k,1)
    end
    return _evaluate_points(solver,comp,k,1)
end

"""
    build_dlp_kress_workspace_full(
        solver::DLP,
        pts::BoundaryPoints{T},
    ) where {T<:Real} → DLPKressWorkspace

Build the complete geometry workspace for a fixed DLP-Kress discretization.

## Arguments
* `solver::DLP`: Smooth or global-corner DLP-Kress solver.
* `pts::BoundaryPoints{T}`: Boundary discretization.

## Returns
* `ws::DLPKressWorkspace`: Full DLP-Kress workspace.
"""
function build_dlp_kress_workspace_full(solver::DLP,pts::BoundaryPoints{T}) where {T<:Real}
    Rmat=build_Rmat_dlp_kress(solver,pts)
    G=boundary_geom_cache(pts,_is_dlp_kress_graded(solver,pts))
    parr=_boundary_panel_arrays_cache(pts)
    N=length(pts)
    return DLPKressWorkspace(Rmat,G,parr,N)
end

"""
    build_dlp_kress_reduced_workspace(
        solver::DLP,
        pts::BoundaryPoints{T},
    ) where {T<:Real} → DLPKressReducedWorkspace

Build the symmetry-reduced DLP-Kress workspace by combining the complete
periodic Kress geometry with the exact [`SymmetryOrbitMap`](@ref) of the active
symmetry.

## Arguments
* `solver::DLP`: DLP-Kress solver with an active symmetry.
* `pts::BoundaryPoints{T}`: Full periodic boundary discretization.

## Returns
* `rws::DLPKressReducedWorkspace`: Reduced workspace containing the full Kress workspace and symmetry-orbit map.
"""
function build_dlp_kress_reduced_workspace(solver::DLP,pts::BoundaryPoints{T}) where {T<:Real}
    isnothing(solver.symmetry)&&throw(ArgumentError("Cannot build a reduced DLP-Kress workspace without an active symmetry"))
    full=build_dlp_kress_workspace_full(solver,pts)
    orbits=symmetry_index_orbits(T,pts,solver.symmetry)
    return DLPKressReducedWorkspace(full,orbits)
end

"""
    build_dlp_kress_workspace(
        solver::DLP,
        pts::BoundaryPoints{T},
    ) where {T<:Real} → DLPKressAnyWorkspace{T}

Build the appropriate DLP-Kress workspace: full when no symmetry is active and
symmetry-reduced otherwise.

## Arguments
* `solver::DLP`: Smooth or global-corner DLP-Kress solver.
* `pts::BoundaryPoints{T}`: Full boundary discretization.

## Returns
* `ws::DLPKressWorkspace` when `solver.symmetry===nothing`.
* `rws::DLPKressReducedWorkspace` when symmetry reduction is active.
"""
function build_dlp_kress_workspace(solver::DLP,pts::BoundaryPoints{T}) where {T<:Real}
    return _dlp_kress_use_reduced(solver) ? build_dlp_kress_reduced_workspace(solver,pts) : build_dlp_kress_workspace_full(solver,pts)
end

@inline _workspace_dim(ws::DLPKressWorkspace)=ws.N
@inline _workspace_dim(ws::DLPKressReducedWorkspace)=fundamental_size(ws.orbits)
@inline function boundary_matrix_size(solver::DLP,pts::BoundaryPoints{T}) where {T<:Real}
    isnothing(solver.symmetry)&&return boundary_matrix_size(pts)
    return fundamental_size(symmetry_index_orbits(T,pts,solver.symmetry))
end

function _dlp_kress_workspace_from_Rmat(solver::DLP,pts::BoundaryPoints{T},Rmat::AbstractMatrix{T}) where {T<:Real}
    G=boundary_geom_cache(pts,_is_dlp_kress_graded(solver,pts))
    parr=_boundary_panel_arrays_cache(pts)
    full=DLPKressWorkspace(Rmat,G,parr,length(pts))
    isnothing(solver.symmetry)&&return full
    return DLPKressReducedWorkspace(full,symmetry_index_orbits(T,pts,solver.symmetry))
end

###############################################
############# NO SYMMETRY PATHWAY #############
###############################################

"""
    construct_dlp_matrix!(
        solver::DLP,
        D::AbstractMatrix{Complex{T}},
        pts::BoundaryPoints{T},
        Rmat::AbstractMatrix{T},
        G::BoundaryGeomCache{T},
        k::T;
        multithreaded::Bool=true,
    ) where {T<:Real} → D

Assemble the Kress-corrected Nyström matrix of the Helmholtz double-layer
operator.

For off-diagonal entries,

    D[i,j]=Rmat[i,j]*l1+pts.ws[j]*l2,

where

    l1=-(k/2π)*inner*J₁(kr)/r,

and

    l2=(ik/2)*inner*H₁^(1)(kr)/r-l1*logterm.

## Arguments
* `solver::DLP`: DLP-Kress solver.
* `D::AbstractMatrix{Complex{T}}`: Preallocated destination matrix.
* `pts::BoundaryPoints{T}`: Boundary discretization.
* `Rmat::AbstractMatrix{T}`: Kress logarithmic correction matrix.
* `G::BoundaryGeomCache{T}`: Pairwise geometry cache.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread the off-diagonal assembly.

## Returns
* `D::AbstractMatrix{Complex{T}}`: Assembled DLP matrix.
"""
function construct_dlp_matrix!(solver::DLP,D::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},Rmat::AbstractMatrix{T},G::BoundaryGeomCache{T},k::T;multithreaded::Bool=true) where {T<:Real}
    αL1=-k*inv_two_pi
    αL2=Complex{T}(0,k/2)
    fill!(D,zero(Complex{T}))
    N=length(pts)
    @inbounds for i in 1:N
        D[i,i]=Complex{T}(pts.ws[i]*G.kappa[i],zero(T))
    end
    @use_threads multithreading=(multithreaded&&N>=32) for j in 2:N
        @inbounds for i in 1:j-1
            r=G.R[i,j]
            invr=G.invR[i,j]
            lt=G.logterm[i,j]
            inn_ij=G.inner[i,j]
            inn_ji=G.inner[j,i]
            _,h1=hankel_pair01(k*r)
            j1=real(h1)
            l1_ij=αL1*inn_ij*j1*invr
            l2_ij=αL2*inn_ij*h1*invr-l1_ij*lt
            D[i,j]=Rmat[i,j]*l1_ij+pts.ws[j]*l2_ij
            l1_ji=αL1*inn_ji*j1*invr
            l2_ji=αL2*inn_ji*h1*invr-l1_ji*lt
            D[j,i]=Rmat[j,i]*l1_ji+pts.ws[i]*l2_ji
        end
    end
    return D
end

# INTERNAL - return the logarithmic and smooth Kress contributions separately.
function construct_dlp_split!(solver::DLP,Dlog::AbstractMatrix{Complex{T}},Dsmooth::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},Rmat::AbstractMatrix{T},G::BoundaryGeomCache{T},parr::BoundaryPanelArrays{T},k::T;multithreaded::Bool=true) where {T<:Real}
    αL1=-k*inv_two_pi
    αL2=Complex{T}(0,k/2)
    fill!(Dlog,zero(Complex{T}))
    fill!(Dsmooth,zero(Complex{T}))
    N=length(parr.X)
    @inbounds for i in 1:N
        Dsmooth[i,i]=Complex{T}(pts.ws[i]*G.kappa[i],zero(T))
    end
    @use_threads multithreading=(multithreaded&&N>=32) for j in 2:N
        @inbounds for i in 1:j-1
            r=G.R[i,j]
            invr=G.invR[i,j]
            lt=G.logterm[i,j]
            inn_ij=G.inner[i,j]
            inn_ji=G.inner[j,i]
            _,h1=hankel_pair01(k*r)
            j1=real(h1)
            l1_ij=αL1*inn_ij*j1*invr
            l2_ij=αL2*inn_ij*h1*invr-l1_ij*lt
            Dlog[i,j]=Rmat[i,j]*l1_ij
            Dsmooth[i,j]=pts.ws[j]*l2_ij
            l1_ji=αL1*inn_ji*j1*invr
            l2_ji=αL2*inn_ji*h1*invr-l1_ji*lt
            Dlog[j,i]=Rmat[j,i]*l1_ji
            Dsmooth[j,i]=pts.ws[i]*l2_ji
        end
    end
    return Dlog,Dsmooth
end

"""
    construct_fredholm_matrix!(
        solver::DLP,
        F::AbstractMatrix{Complex{T}},
        pts::BoundaryPoints{T},
        Rmat::AbstractMatrix{T},
        G::BoundaryGeomCache{T},
        parr::BoundaryPanelArrays{T},
        k::T;
        multithreaded::Bool=true,
    ) where {T<:Real} → F

Assemble the full DLP-Kress Fredholm matrix

    F(k)=I-D(k).

## Arguments
* `solver::DLP`: DLP-Kress solver.
* `F::AbstractMatrix{Complex{T}}`: Preallocated destination matrix.
* `pts::BoundaryPoints{T}`: Boundary discretization.
* `Rmat::AbstractMatrix{T}`: Kress logarithmic correction matrix.
* `G::BoundaryGeomCache{T}`: Pairwise geometry cache.
* `parr::BoundaryPanelArrays{T}`: Flat boundary-array cache.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread the off-diagonal assembly.

## Returns
* `F::AbstractMatrix{Complex{T}}`: Assembled Fredholm matrix.
"""
function construct_fredholm_matrix!(solver::DLP,F::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},Rmat::AbstractMatrix{T},G::BoundaryGeomCache{T},parr::BoundaryPanelArrays{T},k::T;multithreaded::Bool=true) where {T<:Real}
    αL1=-k*inv_two_pi
    αL2=Complex{T}(0,k/2)
    fill!(F,zero(Complex{T}))
    N=length(parr.X)
    @inbounds for i in 1:N
        F[i,i]=one(Complex{T})-Complex{T}(pts.ws[i]*G.kappa[i],zero(T))
    end
    @use_threads multithreading=(multithreaded&&N>=32) for j in 2:N
        @inbounds for i in 1:j-1
            r=G.R[i,j]
            invr=G.invR[i,j]
            lt=G.logterm[i,j]
            inn_ij=G.inner[i,j]
            inn_ji=G.inner[j,i]
            _,h1=hankel_pair01(k*r)
            j1=real(h1)
            l1_ij=αL1*inn_ij*j1*invr
            l2_ij=αL2*inn_ij*h1*invr-l1_ij*lt
            F[i,j]=-(Rmat[i,j]*l1_ij+pts.ws[j]*l2_ij)
            l1_ji=αL1*inn_ji*j1*invr
            l2_ji=αL2*inn_ji*h1*invr-l1_ji*lt
            F[j,i]=-(Rmat[j,i]*l1_ji+pts.ws[i]*l2_ji)
        end
    end
    return F
end

"""
    construct_dlp_matrix_derivatives!(
        solver::DLP,
        D::AbstractMatrix{Complex{T}},
        D1::AbstractMatrix{Complex{T}},
        D2::AbstractMatrix{Complex{T}},
        pts::BoundaryPoints{T},
        Rmat::AbstractMatrix{T},
        G::BoundaryGeomCache{T},
        parr::BoundaryPanelArrays{T},
        k::T;
        multithreaded::Bool=true,
    ) where {T<:Real} → D,D1,D2

Assemble the full Kress-corrected DLP matrix and its first two wavenumber
derivatives.

The diagonal DLP limit is independent of `k`, hence

    D1[i,i]=D2[i,i]=0.

## Arguments
* `solver::DLP`: DLP-Kress solver.
* `D::AbstractMatrix{Complex{T}}`: Destination matrix for `D(k)`.
* `D1::AbstractMatrix{Complex{T}}`: Destination matrix for `D'(k)`.
* `D2::AbstractMatrix{Complex{T}}`: Destination matrix for `D''(k)`.
* `pts::BoundaryPoints{T}`: Boundary discretization.
* `Rmat::AbstractMatrix{T}`: Kress logarithmic correction matrix.
* `G::BoundaryGeomCache{T}`: Pairwise geometry cache.
* `parr::BoundaryPanelArrays{T}`: Flat boundary-array cache.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread the off-diagonal assembly.

## Returns
* `D::AbstractMatrix{Complex{T}}`: DLP matrix.
* `D1::AbstractMatrix{Complex{T}}`: First wavenumber derivative.
* `D2::AbstractMatrix{Complex{T}}`: Second wavenumber derivative.
"""
function construct_dlp_matrix_derivatives!(solver::DLP,D::AbstractMatrix{Complex{T}},D1::AbstractMatrix{Complex{T}},D2::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},Rmat::AbstractMatrix{T},G::BoundaryGeomCache{T},parr::BoundaryPanelArrays{T},k::T;multithreaded::Bool=true) where {T<:Real}
    αL1=-k*inv_two_pi
    αL2=Complex{T}(0,k/2)
    fill!(D,zero(Complex{T}))
    fill!(D1,zero(Complex{T}))
    fill!(D2,zero(Complex{T}))
    N=length(parr.X)
    @inbounds for i in 1:N
        D[i,i]=Complex{T}(pts.ws[i]*G.kappa[i],zero(T))
        D1[i,i]=zero(Complex{T})
        D2[i,i]=zero(Complex{T})
    end
    @use_threads multithreading=(multithreaded&&N>=32) for j in 2:N
        @inbounds for i in 1:j-1
            r=G.R[i,j]
            invr=G.invR[i,j]
            lt=G.logterm[i,j]
            inn_ij=G.inner[i,j]
            inn_ji=G.inner[j,i]
            kr=k*r
            h0,h1=hankel_pair01(kr)
            j0=real(h0)
            j1=real(h1)
            l1_ij=αL1*inn_ij*j1*invr
            l2_ij=αL2*inn_ij*h1*invr-l1_ij*lt
            D[i,j]=Rmat[i,j]*l1_ij+pts.ws[j]*l2_ij
            l1_ij_1=-(inn_ij*k*j0)*inv_two_pi
            l1_ij_2=(inn_ij*(k*r*j1-j0))*inv_two_pi
            l2_ij_1=(inn_ij*k*(lt*j0+im*pi*h0))*inv_two_pi
            l2_ij_2=(inn_ij*(lt*(j0-k*r*j1)+im*pi*(h0-k*r*h1)))*inv_two_pi
            D1[i,j]=Rmat[i,j]*l1_ij_1+pts.ws[j]*l2_ij_1
            D2[i,j]=Rmat[i,j]*l1_ij_2+pts.ws[j]*l2_ij_2
            l1_ji=αL1*inn_ji*j1*invr
            l2_ji=αL2*inn_ji*h1*invr-l1_ji*lt
            D[j,i]=Rmat[j,i]*l1_ji+pts.ws[i]*l2_ji
            l1_ji_1=-(inn_ji*k*j0)*inv_two_pi
            l1_ji_2=(inn_ji*(k*r*j1-j0))*inv_two_pi
            l2_ji_1=(inn_ji*k*(lt*j0+im*pi*h0))*inv_two_pi
            l2_ji_2=(inn_ji*(lt*(j0-k*r*j1)+im*pi*(h0-k*r*h1)))*inv_two_pi
            D1[j,i]=Rmat[j,i]*l1_ji_1+pts.ws[i]*l2_ji_1
            D2[j,i]=Rmat[j,i]*l1_ji_2+pts.ws[i]*l2_ji_2
        end
    end
    return D,D1,D2
end

"""
    construct_fredholm_matrix_derivatives!(
        solver::DLP,
        F::AbstractMatrix{Complex{T}},
        F1::AbstractMatrix{Complex{T}},
        F2::AbstractMatrix{Complex{T}},
        pts::BoundaryPoints{T},
        Rmat::AbstractMatrix{T},
        G::BoundaryGeomCache{T},
        parr::BoundaryPanelArrays{T},
        k::T;
        multithreaded::Bool=true,
    ) where {T<:Real} → F,F1,F2

Assemble the full DLP-Kress Fredholm matrix and its first two derivatives,

    F=I-D,
    F1=-D1,
    F2=-D2.

## Arguments
* `solver::DLP`: DLP-Kress solver.
* `F::AbstractMatrix{Complex{T}}`: Destination matrix for `F(k)`.
* `F1::AbstractMatrix{Complex{T}}`: Destination matrix for `F'(k)`.
* `F2::AbstractMatrix{Complex{T}}`: Destination matrix for `F''(k)`.
* `pts::BoundaryPoints{T}`: Boundary discretization.
* `Rmat::AbstractMatrix{T}`: Kress logarithmic correction matrix.
* `G::BoundaryGeomCache{T}`: Pairwise geometry cache.
* `parr::BoundaryPanelArrays{T}`: Flat boundary-array cache.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread the DLP assembly.

## Returns
* `F::AbstractMatrix{Complex{T}}`: Fredholm matrix.
* `F1::AbstractMatrix{Complex{T}}`: First derivative.
* `F2::AbstractMatrix{Complex{T}}`: Second derivative.
"""
function construct_fredholm_matrix_derivatives!(solver::DLP,F::AbstractMatrix{Complex{T}},F1::AbstractMatrix{Complex{T}},F2::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},Rmat::AbstractMatrix{T},G::BoundaryGeomCache{T},parr::BoundaryPanelArrays{T},k::T;multithreaded::Bool=true) where {T<:Real}
    construct_dlp_matrix_derivatives!(solver,F,F1,F2,pts,Rmat,G,parr,k;multithreaded=multithreaded)
    @inbounds for j in axes(F,2),i in axes(F,1)
        F[i,j]*=-1
        F1[i,j]*=-1
        F2[i,j]*=-1
    end
    @inbounds for i in axes(F,1)
        F[i,i]+=one(Complex{T})
    end
    return F,F1,F2
end

#################################################
############# DESYMMETRIZED PATHWAY ############
#################################################

"""
    construct_dlp_matrix!(solver::DLP,D::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},rws::DLPKressReducedWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real} → D

Assemble the symmetry-reduced DLP-Kress matrix by folding the complete discrete
full-boundary Kress operator over each source symmetry orbit.

For fundamental target index `a` and reduced source index `b`,

    i=Ifund[a],

and each symmetry image of `b` is

    j_l=fund_to_full[l,b],
    χ_l=fund_to_scale[l,b].

The reduced matrix is

    Dred[a,b]=Σ_l χ_l Dfull[i,j_l].

Every image therefore uses exactly the same Kress product-integration entry as
the corresponding full matrix. In particular, nonidentity images on the same
periodic component must not be replaced by ordinary trapezoidal DLP entries.

The diagonal full-grid entry is

    Dfull[i,i]=ws[i]*κ[i].

## Arguments
* `solver::DLP`: DLP-Kress solver with active symmetry.
* `D::AbstractMatrix{Complex{T}}`: Preallocated reduced destination matrix.
* `pts::BoundaryPoints{T}`: Full periodic boundary discretization.
* `rws::DLPKressReducedWorkspace{T}`: Full Kress workspace and symmetry-orbit map.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread reduced-column assembly.

## Returns
* `D::AbstractMatrix{Complex{T}}`: Symmetry-reduced DLP-Kress matrix.
"""
function construct_dlp_matrix!(solver::DLP,D::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},rws::DLPKressReducedWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    full=rws.full
    orbits=rws.orbits
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    ng=orbit_size(orbits)
    @assert size(D)==(m,m)
    Rmat=full.Rmat
    G=full.G
    αL1=-k*inv_two_pi
    αL2=Complex{T}(0,k/2)
    fill!(D,zero(Complex{T}))
    @use_threads multithreading=(multithreaded&&m>=32) for b in 1:m
        @inbounds for a in 1:m
            i=Ifund[a]
            acc=zero(Complex{T})
            for l in 1:ng
                j=orbits.fund_to_full[l,b]
                χ=orbits.fund_to_scale[l,b]
                if i==j
                    dval=Complex{T}(pts.ws[i]*G.kappa[i],zero(T))
                else
                    r=G.R[i,j]
                    invr=G.invR[i,j]
                    lt=G.logterm[i,j]
                    inn=G.inner[i,j]
                    _,h1=hankel_pair01(k*r)
                    j1=real(h1)
                    l1=αL1*inn*j1*invr
                    l2=αL2*inn*h1*invr-l1*lt
                    dval=Rmat[i,j]*l1+pts.ws[j]*l2
                end
                acc+=χ*dval
            end
            D[a,b]=acc
        end
    end
    return D
end

"""
    construct_fredholm_matrix!(solver::DLP,F::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},rws::DLPKressReducedWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real} → F

Assemble the symmetry-reduced DLP-Kress Fredholm matrix

    F(k)=I-Dred(k),

where `Dred` is obtained by exact folding of the full discrete Kress operator
over the source symmetry orbits.

## Arguments
* `solver::DLP`: DLP-Kress solver with active symmetry.
* `F::AbstractMatrix{Complex{T}}`: Preallocated reduced destination matrix.
* `pts::BoundaryPoints{T}`: Full periodic boundary discretization.
* `rws::DLPKressReducedWorkspace{T}`: Full Kress workspace and symmetry-orbit map.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread reduced DLP assembly.

## Returns
* `F::AbstractMatrix{Complex{T}}`: Symmetry-reduced Fredholm matrix.
"""
function construct_fredholm_matrix!(solver::DLP,F::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},rws::DLPKressReducedWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    construct_dlp_matrix!(solver,F,pts,rws,k;multithreaded=multithreaded)
    @inbounds for j in axes(F,2),i in axes(F,1)
        F[i,j]*=-1
    end
    @inbounds for i in axes(F,1)
        F[i,i]+=one(Complex{T})
    end
    return F
end

"""
    construct_dlp_matrix_derivatives!(solver::DLP,D::AbstractMatrix{Complex{T}},D1::AbstractMatrix{Complex{T}},D2::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},rws::DLPKressReducedWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real} → D,D1,D2

Assemble the symmetry-reduced DLP-Kress matrix and its first two wavenumber
derivatives by folding the complete full-grid Kress discretization.

For every source symmetry image,

    Dred[a,b]  =Σ_l χ_l Dfull[i,j_l],
    Dred'[a,b] =Σ_l χ_l Dfull'[i,j_l],
    Dred''[a,b]=Σ_l χ_l Dfull''[i,j_l].

Thus all orbit images use the same Kress product-integration formula as their
corresponding full-matrix entries.

The analytic diagonal DLP limit is independent of `k`, hence

    Dfull'[i,i]=Dfull''[i,i]=0.

## Arguments
* `solver::DLP`: DLP-Kress solver with active symmetry.
* `D::AbstractMatrix{Complex{T}}`: Destination matrix for `Dred(k)`.
* `D1::AbstractMatrix{Complex{T}}`: Destination matrix for `Dred'(k)`.
* `D2::AbstractMatrix{Complex{T}}`: Destination matrix for `Dred''(k)`.
* `pts::BoundaryPoints{T}`: Full periodic boundary discretization.
* `rws::DLPKressReducedWorkspace{T}`: Full Kress workspace and symmetry-orbit map.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread reduced-column assembly.

## Returns
* `D::AbstractMatrix{Complex{T}}`: Reduced DLP matrix.
* `D1::AbstractMatrix{Complex{T}}`: First wavenumber derivative.
* `D2::AbstractMatrix{Complex{T}}`: Second wavenumber derivative.
"""
function construct_dlp_matrix_derivatives!(solver::DLP,D::AbstractMatrix{Complex{T}},D1::AbstractMatrix{Complex{T}},D2::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},rws::DLPKressReducedWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    full=rws.full
    orbits=rws.orbits
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    ng=orbit_size(orbits)
    @assert size(D)==(m,m)
    @assert size(D1)==(m,m)
    @assert size(D2)==(m,m)
    Rmat=full.Rmat
    G=full.G
    αL1=-k*inv_two_pi
    αL2=Complex{T}(0,k/2)
    fill!(D,zero(Complex{T}))
    fill!(D1,zero(Complex{T}))
    fill!(D2,zero(Complex{T}))
    @use_threads multithreading=(multithreaded&&m>=32) for b in 1:m
        @inbounds for a in 1:m
            i=Ifund[a]
            acc0=zero(Complex{T})
            acc1=zero(Complex{T})
            acc2=zero(Complex{T})
            for l in 1:ng
                j=orbits.fund_to_full[l,b]
                χ=orbits.fund_to_scale[l,b]
                if i==j
                    dval=Complex{T}(pts.ws[i]*G.kappa[i],zero(T))
                    dval1=zero(Complex{T})
                    dval2=zero(Complex{T})
                else
                    r=G.R[i,j]
                    invr=G.invR[i,j]
                    lt=G.logterm[i,j]
                    inn=G.inner[i,j]
                    kr=k*r
                    h0,h1=hankel_pair01(kr)
                    j0=real(h0)
                    j1=real(h1)
                    l1=αL1*inn*j1*invr
                    l2=αL2*inn*h1*invr-l1*lt
                    dval=Rmat[i,j]*l1+pts.ws[j]*l2
                    l1_1=-(inn*k*j0)*inv_two_pi
                    l1_2=(inn*(kr*j1-j0))*inv_two_pi
                    l2_1=(inn*k*(lt*j0+im*pi*h0))*inv_two_pi
                    l2_2=(inn*(lt*(j0-kr*j1)+im*pi*(h0-kr*h1)))*inv_two_pi
                    dval1=Rmat[i,j]*l1_1+pts.ws[j]*l2_1
                    dval2=Rmat[i,j]*l1_2+pts.ws[j]*l2_2
                end
                acc0+=χ*dval
                acc1+=χ*dval1
                acc2+=χ*dval2
            end
            D[a,b]=acc0
            D1[a,b]=acc1
            D2[a,b]=acc2
        end
    end
    return D,D1,D2
end

"""
    construct_fredholm_matrix_derivatives!(solver::DLP,F::AbstractMatrix{Complex{T}},F1::AbstractMatrix{Complex{T}},F2::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},rws::DLPKressReducedWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real} → F,F1,F2

Assemble the symmetry-reduced DLP-Kress Fredholm matrix and its first two
wavenumber derivatives,

    F=I-Dred,
    F1=-Dred',
    F2=-Dred''.

## Arguments
* `solver::DLP`: DLP-Kress solver with active symmetry.
* `F::AbstractMatrix{Complex{T}}`: Destination matrix for `F(k)`.
* `F1::AbstractMatrix{Complex{T}}`: Destination matrix for `F'(k)`.
* `F2::AbstractMatrix{Complex{T}}`: Destination matrix for `F''(k)`.
* `pts::BoundaryPoints{T}`: Full periodic boundary discretization.
* `rws::DLPKressReducedWorkspace{T}`: Full Kress workspace and symmetry-orbit map.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread reduced DLP assembly.

## Returns
* `F::AbstractMatrix{Complex{T}}`: Reduced Fredholm matrix.
* `F1::AbstractMatrix{Complex{T}}`: First derivative.
* `F2::AbstractMatrix{Complex{T}}`: Second derivative.
"""
function construct_fredholm_matrix_derivatives!(solver::DLP,F::AbstractMatrix{Complex{T}},F1::AbstractMatrix{Complex{T}},F2::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},rws::DLPKressReducedWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    construct_dlp_matrix_derivatives!(solver,F,F1,F2,pts,rws,k;multithreaded=multithreaded)
    @inbounds for j in axes(F,2),i in axes(F,1)
        F[i,j]*=-1
        F1[i,j]*=-1
        F2[i,j]*=-1
    end
    @inbounds for i in axes(F,1)
        F[i,i]+=one(Complex{T})
    end
    return F,F1,F2
end

########################################
######### NEEDED FOR HUSIMIS ###########
########################################

"""
    adjoint_fredholm_matrix!(
        A::AbstractMatrix{Complex{T}},
        D::AbstractMatrix{Complex{T}},
        solver::DLP,
        pts::BoundaryPoints{T},
        ws::DLPKressAnyWorkspace{T},
        k::T;
        multithreaded::Bool=true,
    ) where {T<:Real} → A

Assemble the adjoint DLP-Kress Fredholm matrix

    A=I-W⁻¹DᵀW,

with `W=diag(pts.ds)`.

For a symmetry-reduced workspace the weighted transpose is applied to the
fundamental-domain discretization.

## Arguments
* `A::AbstractMatrix{Complex{T}}`: Preallocated destination matrix for the adjoint Fredholm operator.
* `D::AbstractMatrix{Complex{T}}`: Preallocated DLP workspace matrix.
* `solver::DLP`: DLP-Kress solver.
* `pts::BoundaryPoints{T}`: Full boundary discretization.
* `ws::DLPKressAnyWorkspace{T}`: Full or symmetry-reduced workspace.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread DLP assembly.

## Returns
* `A::AbstractMatrix{Complex{T}}`: Adjoint Fredholm matrix.
"""
function adjoint_fredholm_matrix!(A::AbstractMatrix{Complex{T}},D::AbstractMatrix{Complex{T}},solver::DLP,pts::BoundaryPoints{T},ws::DLPKressWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    N=ws.N
    construct_dlp_matrix!(solver,D,pts,ws.Rmat,ws.G,k;multithreaded=multithreaded)
    ds=pts.ds
    @inbounds for i in 1:N,j in 1:N
        A[i,j]=-D[j,i]*ds[j]/ds[i]
    end
    @inbounds for i in 1:N
        A[i,i]+=one(Complex{T})
    end
    return A
end

function adjoint_fredholm_matrix!(A::AbstractMatrix{Complex{T}},D::AbstractMatrix{Complex{T}},solver::DLP,pts::BoundaryPoints{T},rws::DLPKressReducedWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    orbits=rws.orbits
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    construct_dlp_matrix!(solver,D,pts,rws,k;multithreaded=multithreaded)
    ds=pts.ds
    @inbounds for b in 1:m,a in 1:m
        i=Ifund[a]
        j=Ifund[b]
        A[a,b]=-D[b,a]*ds[j]/ds[i]
    end
    @inbounds for a in 1:m
        A[a,a]+=one(Complex{T})
    end
    return A
end

function adjoint_fredholm_matrix!(A::AbstractMatrix{Complex{T}},D::AbstractMatrix{Complex{T}},solver::DLP,pts::BoundaryPoints{T},k::T;multithreaded::Bool=true) where {T<:Real}
    ws=build_dlp_kress_workspace(solver,pts)
    return adjoint_fredholm_matrix!(A,D,solver,pts,ws,k;multithreaded=multithreaded)
end

function adjoint_fredholm_matrix!(A::AbstractMatrix{Complex{T}},D::AbstractMatrix{Complex{T}},solver::DLP,pts::BoundaryPoints{T},Rmat::AbstractMatrix{T},k::T;multithreaded::Bool=true) where {T<:Real}
    ws=_dlp_kress_workspace_from_Rmat(solver,pts,Rmat)
    return adjoint_fredholm_matrix!(A,D,solver,pts,ws,k;multithreaded=multithreaded)
end

########################################
######## MATRIX ASSEMBLY INTERFACE #####
########################################

"""
    construct_matrices!(
        solver::DLP,
        A::AbstractMatrix{Complex{T}},
        pts::BoundaryPoints{T},
        ws::DLPKressAnyWorkspace{T},
        k::T;
        multithreaded::Bool=true,
    ) where {T<:Real}

    construct_matrices!(
        solver::DLP,
        A::AbstractMatrix{Complex{T}},
        A1::AbstractMatrix{Complex{T}},
        A2::AbstractMatrix{Complex{T}},
        pts::BoundaryPoints{T},
        ws::DLPKressAnyWorkspace{T},
        k::T;
        multithreaded::Bool=true,
    ) where {T<:Real}

Assemble

    A(k)=I-D(k),

and, when requested,

    A1(k)=-D'(k),
    A2(k)=-D''(k).

## Arguments
* `solver::DLP`: DLP-Kress solver.
* `A::AbstractMatrix{Complex{T}}`: Destination matrix for the Fredholm operator.
* `A1::AbstractMatrix{Complex{T}}`: Destination matrix for the first derivative.
* `A2::AbstractMatrix{Complex{T}}`: Destination matrix for the second derivative.
* `pts::BoundaryPoints{T}`: Boundary discretization.
* `ws::DLPKressAnyWorkspace{T}`: Full or symmetry-reduced workspace.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread low-level matrix assembly.

## Returns
* `A::AbstractMatrix{Complex{T}}` for the single-matrix overload.
* `(A,A1,A2)` for the derivative overload.
"""
function construct_matrices!(solver::DLP,A::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},ws::DLPKressWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    return construct_fredholm_matrix!(solver,A,pts,ws.Rmat,ws.G,ws.parr,k;multithreaded=multithreaded)
end

function construct_matrices!(solver::DLP,A::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},rws::DLPKressReducedWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    return construct_fredholm_matrix!(solver,A,pts,rws,k;multithreaded=multithreaded)
end

function construct_matrices!(solver::DLP,A::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},Rmat::AbstractMatrix{T},k::T;multithreaded::Bool=true) where {T<:Real}
    ws=_dlp_kress_workspace_from_Rmat(solver,pts,Rmat)
    return construct_matrices!(solver,A,pts,ws,k;multithreaded=multithreaded)
end

function construct_matrices!(solver::DLP,A::AbstractMatrix{Complex{T}},A1::AbstractMatrix{Complex{T}},A2::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},ws::DLPKressWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    return construct_fredholm_matrix_derivatives!(solver,A,A1,A2,pts,ws.Rmat,ws.G,ws.parr,k;multithreaded=multithreaded)
end

function construct_matrices!(solver::DLP,A::AbstractMatrix{Complex{T}},A1::AbstractMatrix{Complex{T}},A2::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},rws::DLPKressReducedWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    return construct_fredholm_matrix_derivatives!(solver,A,A1,A2,pts,rws,k;multithreaded=multithreaded)
end

function construct_matrices!(solver::DLP,A::AbstractMatrix{Complex{T}},A1::AbstractMatrix{Complex{T}},A2::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},Rmat::AbstractMatrix{T},k::T;multithreaded::Bool=true) where {T<:Real}
    ws=_dlp_kress_workspace_from_Rmat(solver,pts,Rmat)
    return construct_matrices!(solver,A,A1,A2,pts,ws,k;multithreaded=multithreaded)
end

function construct_matrices!(solver::DLP,basis::AbstractHankelBasis,A::AbstractMatrix{Complex{T}},dA::AbstractMatrix{Complex{T}},ddA::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},k::T;multithreaded::Bool=true) where {T<:Real}
    ws=build_dlp_kress_workspace(solver,pts)
    construct_matrices!(solver,A,dA,ddA,pts,ws,k;multithreaded=multithreaded)
    return A,dA,ddA
end

function construct_matrices!(solver::DLP,basis::AbstractHankelBasis,A::AbstractMatrix{Complex{T}},dA::AbstractMatrix{Complex{T}},ddA::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},ws::DLPKressAnyWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    construct_matrices!(solver,A,dA,ddA,pts,ws,k;multithreaded=multithreaded)
    return A,dA,ddA
end

function construct_matrices!(solver::DLP,basis::AbstractHankelBasis,A::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},rws::DLPKressReducedWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    construct_matrices!(solver,A,pts,rws,k;multithreaded=multithreaded)
    return A
end

function construct_matrices(solver::DLP,basis::AbstractHankelBasis,pts::BoundaryPoints{T},rws::DLPKressReducedWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    n=_workspace_dim(rws)
    A=Matrix{Complex{T}}(undef,n,n)
    construct_matrices!(solver,basis,A,pts,rws,k;multithreaded=multithreaded)
    return A
end

########################################
############### SOLVE ##################
########################################

"""
    solve(
        solver::DLP,
        basis::Ba,
        pts::BoundaryPoints{T},
        k;
        multithreaded::Bool=true,
        use_krylov::Bool=true,
        which::Symbol=:det,
    ) where {T<:Real,Ba<:AbsBasis}

    solve(
        solver::DLP,
        basis::Ba,
        pts::BoundaryPoints{T},
        ws::DLPKressAnyWorkspace{T},
        k;
        multithreaded::Bool=true,
        use_krylov::Bool=true,
        which::Symbol=:det_argmin,
    ) where {T<:Real,Ba<:AbsBasis}

Assemble the DLP-Kress Fredholm matrix and evaluate the scalar spectral
diagnostic selected by `which`.

## Arguments
* `solver::DLP`: DLP-Kress solver.
* `basis::Ba`: Basis argument retained for the common solver interface.
* `pts::BoundaryPoints{T}`: Boundary discretization.
* `ws::DLPKressAnyWorkspace{T}`: Optional full or symmetry-reduced workspace.
* `k`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread matrix assembly.
* `use_krylov::Bool`: Whether the common spectral backend may use its Krylov pathway.
* `which::Symbol`: Scalar spectral diagnostic to evaluate.

## Returns
* Scalar spectral diagnostic selected by `which`.
"""
function solve(solver::DLP,basis::Ba,pts::BoundaryPoints{T},k;multithreaded::Bool=true,use_krylov::Bool=true,which::Symbol=:det) where {T<:Real,Ba<:AbsBasis}
    ws=build_dlp_kress_workspace(solver,pts)
    n=_workspace_dim(ws)
    A=Matrix{Complex{T}}(undef,n,n)
    @blas_1 construct_matrices!(solver,A,pts,ws,k;multithreaded=multithreaded)
    @svd_or_det_solve A use_krylov which MAX_BLAS_THREADS
end

function solve(solver::DLP,basis::Ba,pts::BoundaryPoints{T},ws::DLPKressAnyWorkspace{T},k;multithreaded::Bool=true,use_krylov::Bool=true,which::Symbol=:det_argmin) where {T<:Real,Ba<:AbsBasis}
    n=_workspace_dim(ws)
    A=Matrix{Complex{T}}(undef,n,n)
    @blas_1 construct_matrices!(solver,A,pts,ws,k;multithreaded=multithreaded)
    @svd_or_det_solve A use_krylov which MAX_BLAS_THREADS
end

function solve(solver::DLP,basis::Ba,A::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},k,Rmat::AbstractMatrix{T};multithreaded::Bool=true,use_krylov::Bool=true,which::Symbol=:det_argmin) where {T<:Real,Ba<:AbsBasis}
    @blas_1 construct_matrices!(solver,A,pts,Rmat,k;multithreaded=multithreaded)
    @svd_or_det_solve A use_krylov which MAX_BLAS_THREADS
end

function solve(solver::DLP,basis::Ba,A::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},ws::DLPKressAnyWorkspace{T},k;multithreaded::Bool=true,use_krylov::Bool=true,which::Symbol=:det_argmin) where {T<:Real,Ba<:AbsBasis}
    @blas_1 construct_matrices!(solver,A,pts,ws,k;multithreaded=multithreaded)
    @svd_or_det_solve A use_krylov which MAX_BLAS_THREADS
end

########################################
############## SOLVE VECT ##############
########################################

"""
    solve_vect(
        solver::DLP,
        basis::Ba,
        pts::BoundaryPoints{T},
        k;
        multithreaded::Bool=true,
        tol=1e-12,
        maxiter::Int=2000,
        krylovdim::Int=40,
    ) where {T<:Real,Ba<:AbsBasis}

    solve_vect(
        solver::DLP,
        basis::Ba,
        pts::BoundaryPoints{T},
        ws::DLPKressAnyWorkspace{T},
        k;
        multithreaded::Bool=true,
        tol=1e-12,
        maxiter::Int=2000,
        krylovdim::Int=40,
    ) where {T<:Real,Ba<:AbsBasis}

Compute the near-null vector of the adjoint DLP-Kress Fredholm matrix.

The returned vector represents the boundary normal derivative used for Husimi
and related boundary-function reconstruction.

## Arguments
* `solver::DLP`: DLP-Kress solver.
* `basis::Ba`: Basis argument retained for the common solver interface.
* `pts::BoundaryPoints{T}`: Boundary discretization.
* `ws::DLPKressAnyWorkspace{T}`: Optional full or symmetry-reduced workspace.
* `k`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread matrix assembly.
* `tol`: Krylov convergence tolerance.
* `maxiter::Int`: Maximum number of Krylov iterations.
* `krylovdim::Int`: Krylov subspace dimension.

## Returns
* `σ`: Smallest-eigenvalue or near-null residual proxy.
* `u`: Corresponding normalized near-null boundary vector.
"""
function solve_vect(solver::DLP,basis::Ba,A::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},k,Rmat::AbstractMatrix{T};multithreaded::Bool=true,tol=1e-12,maxiter::Int=2000,krylovdim::Int=40) where {T<:Real,Ba<:AbsBasis}
    ws=_dlp_kress_workspace_from_Rmat(solver,pts,Rmat)
    return solve_vect(solver,basis,A,pts,ws,k;multithreaded=multithreaded,tol=tol,maxiter=maxiter,krylovdim=krylovdim)
end

function solve_vect(solver::DLP,basis::Ba,A::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},ws::DLPKressAnyWorkspace{T},k;multithreaded::Bool=true,tol=1e-12,maxiter::Int=2000,krylovdim::Int=40) where {T<:Real,Ba<:AbsBasis}
    D=similar(A)
    @blas_1 adjoint_fredholm_matrix!(A,D,solver,pts,ws,k;multithreaded=multithreaded)
    σ,u,_=smallest_nullvec_krylov!(A;nev=1,tol=tol,maxiter=maxiter,krylovdim=krylovdim)
    return σ,u
end

function solve_vect(solver::DLP,basis::Ba,pts::BoundaryPoints{T},ws::DLPKressAnyWorkspace{T},k;multithreaded::Bool=true,tol=1e-12,maxiter::Int=2000,krylovdim::Int=40) where {T<:Real,Ba<:AbsBasis}
    n=_workspace_dim(ws)
    A=Matrix{Complex{T}}(undef,n,n)
    return solve_vect(solver,basis,A,pts,ws,k;multithreaded=multithreaded,tol=tol,maxiter=maxiter,krylovdim=krylovdim)
end

function solve_vect(solver::DLP,basis::Ba,pts::BoundaryPoints{T},k;multithreaded::Bool=true,tol=1e-12,maxiter::Int=2000,krylovdim::Int=40) where {T<:Real,Ba<:AbsBasis}
    ws=build_dlp_kress_workspace(solver,pts)
    return solve_vect(solver,basis,pts,ws,k;multithreaded=multithreaded,tol=tol,maxiter=maxiter,krylovdim=krylovdim)
end

########################################
############## SOLVE INFO ##############
########################################

# INTERNAL - for checking allocation patterns and execution time of the single-k solve variants.
function solve_INFO(solver::DLP,basis::Ba,pts::BoundaryPoints{T},ws::DLPKressAnyWorkspace{T},k;multithreaded::Bool=true,use_krylov::Bool=true,which::Symbol=:det_argmin) where {T<:Real,Ba<:AbsBasis}
    N=_workspace_dim(ws)
    A=Matrix{Complex{T}}(undef,N,N)
    t0=time()
    @info "Building boundary operator A from cached DLP-Kress workspace..."
    @blas_1 construct_matrices!(solver,A,pts,ws,k;multithreaded=multithreaded)
    t1=time()
    cA=cond(A)
    @info "Condition number of A: $(round(cA;sigdigits=4))"
    t2=time()
    s=@svd_or_det_solve A use_krylov which MAX_BLAS_THREADS
    t3=time()
    build_A=t1-t0
    svd_time=t3-t2
    total=build_A+svd_time
    println("────────── SOLVE_INFO SUMMARY ──────────")
    println("A-matrix build: ",100*build_A/total," %")
    println("SVD: ",100*svd_time/total," %")
    println("(total: ",total," s)")
    println("────────────────────────────────────────")
    return s
end