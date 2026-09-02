const TWO_PI=2*pi

# TODO Baoling Xie and Jun Lai,
# A Singularity Guided Nyström Method for Elastostatics on Two Dimensional
# Domains with Corners, arXiv:2512.18208, 2025

"""
    BoundaryIntegralMethod{T,Sym} <: SweepSolver

Configuration object for the standard boundary integral method (BIM) Fredholm
formulation based on the direct Helmholtz double-layer kernel.

## Description
The assembled Fredholm operator is

    A(k)=I-K(k),

where `K(k)` denotes the Nyström discretization of the doubled interior
Helmholtz double-layer operator.

When a symmetry is active, the complete boundary is still discretized. An exact
[`SymmetryOrbitMap`](@ref) is then used to fold the full source sum onto the
fundamental boundary indices.

## Attributes
* `dim_scaling_factor::T`: Compatibility field for the generic solver infrastructure.
* `pts_scaling_factor::Vector{T}`: Boundary-resolution scaling factors.
* `sampler::Vector`: Sampling rules used on the boundary curves.
* `eps::T`: Numerical tolerance placeholder.
* `min_dim::Int64`: Compatibility field mirroring the other solvers.
* `min_pts::Int64`: Minimum number of boundary points per component.
* `symmetry::Sym`: Optional symmetry descriptor.
"""
struct BoundaryIntegralMethod{T<:Real,Sym}<:SweepSolver
    dim_scaling_factor::T
    pts_scaling_factor::Vector{T}
    sampler::Vector
    eps::T
    min_dim::Int64
    min_pts::Int64
    symmetry::Sym
end

"""
    AbstractHankelBasis <: AbsBasis

Compatibility placeholder used by the direct boundary integral solver.
"""
struct AbstractHankelBasis <: AbsBasis end

"""
    resize_basis(
        basis::Ba,
        billiard::Bi,
        dim::Int,
        k,
    ) where {Ba<:AbstractHankelBasis,Bi<:BilliardGeometry.AbsBilliard} → AbstractHankelBasis

Return an empty basis compatibility object.

## Arguments
* `basis::Ba`: Existing basis placeholder.
* `billiard::Bi`: Billiard geometry.
* `dim::Int`: Requested basis dimension.
* `k`: Wavenumber.

## Returns
* `basis::AbstractHankelBasis`: New compatibility basis placeholder.
"""
function resize_basis(basis::Ba,billiard::Bi,dim::Int,k) where {Ba<:AbstractHankelBasis,Bi<:BilliardGeometry.AbsBilliard}
    return AbstractHankelBasis()
end

function BoundaryIntegralMethod(pts_scaling_factor::Union{T,Vector{T}},billiard::Bi;min_pts=20,symmetry::Union{Nothing,AbsSymmetry}=nothing) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    bs=pts_scaling_factor isa T ? [pts_scaling_factor] : pts_scaling_factor
    sampler=[BilliardGeometry.LinearNodes()]
    Sym=typeof(symmetry)
    return BoundaryIntegralMethod{T,Sym}(one(T),bs,sampler,eps(T),min_pts,min_pts,symmetry)
end

function BoundaryIntegralMethod(pts_scaling_factor::Union{T,Vector{T}},samplers::Vector,billiard::Bi;min_pts=20,symmetry::Union{Nothing,AbsSymmetry}=nothing) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    bs=pts_scaling_factor isa T ? [pts_scaling_factor] : pts_scaling_factor
    Sym=typeof(symmetry)
    return BoundaryIntegralMethod{T,Sym}(one(T),bs,samplers,eps(T),min_pts,min_pts,symmetry)
end

@inline function _dlp_symmetry_orbits(solver::BoundaryIntegralMethod,pts::BoundaryPoints{T}) where {T<:Real}
    isnothing(solver.symmetry)&&return nothing
    return symmetry_index_orbits(T,pts,solver.symmetry)
end

@inline _dlp_matrix_dim(pts::BoundaryPoints,::Nothing)=length(pts)
@inline _dlp_matrix_dim(pts::BoundaryPoints,orbits::SymmetryOrbitMap)=fundamental_size(orbits)
@inline function boundary_matrix_size(solver::BoundaryIntegralMethod,pts::BoundaryPoints{T}) where {T<:Real}
    orbits=_dlp_symmetry_orbits(solver,pts)
    return _dlp_matrix_dim(pts,orbits)
end

# Discretize the supplied physical boundary curves and concatenate them into one BoundaryPoints object.
function _evaluate_bim_curves(solver::BoundaryIntegralMethod,curves::AbstractVector,bs::AbstractVector,samplers::AbstractVector,k)
    T=eltype(solver.pts_scaling_factor)
    Ns=_determine_bp_sizes(curves,bs,k)
    M=length(Ns)
    xy_all=Vector{Vector{SVector{2,T}}}(undef,M)
    normal_all=Vector{Vector{SVector{2,T}}}(undef,M)
    s_all=Vector{Vector{T}}(undef,M)
    ds_all=Vector{Vector{T}}(undef,M)
    curvature_all=Vector{Vector{T}}(undef,M)
    L0=zero(T)
    @inbounds for i in eachindex(curves)
        crv=curves[i]
        t,dt=sample_points(samplers[i],Ns[i])
        xy,normal,s,ds=boundary_coords(crv,t,dt)
        xy_all[i]=xy
        normal_all[i]=normal
        s_all[i]=s.+L0
        ds_all[i]=ds
        curvature_all[i]=curvature(crv,t)
        L0+=crv.length
    end
    return BoundaryPoints(vcat(xy_all...);normal=vcat(normal_all...),s=vcat(s_all...),ds=vcat(ds_all...),curvature=vcat(curvature_all...))
end

# Expand a first-quadrant fundamental boundary into the complete D2-symmetric boundary in canonical CCW order.
function _expand_bim_boundary(pts::BoundaryPoints{T},::BilliardGeometry.XYAxisReflection,L::T) where {T<:Real}
    isempty(pts)&&throw(ArgumentError("Cannot expand an empty fundamental boundary"))
    sf=pts.s.-pts.s[1]
    sy=BilliardGeometry.YAxisReflection()
    sxy=BilliardGeometry.XYAxisReflection()
    sx=BilliardGeometry.XAxisReflection()
    xy=vcat(pts.xy,BilliardGeometry.apply_symmetry(sy,reverse(pts.xy)),BilliardGeometry.apply_symmetry(sxy,pts.xy),BilliardGeometry.apply_symmetry(sx,reverse(pts.xy)))
    normal=vcat(pts.normal,BilliardGeometry.apply_symmetry(sy,reverse(pts.normal)),BilliardGeometry.apply_symmetry(sxy,pts.normal),BilliardGeometry.apply_symmetry(sx,reverse(pts.normal)))
    s=vcat(sf,2*L.-reverse(sf),2*L.+sf,4*L.-reverse(sf))
    ds=vcat(pts.ds,reverse(pts.ds),pts.ds,reverse(pts.ds))
    curvature=vcat(pts.curvature,reverse(pts.curvature),pts.curvature,reverse(pts.curvature))
    return BoundaryPoints(xy;normal=normal,s=s,ds=ds,curvature=curvature)
end

# Expand a half-boundary fundamental domain across the x-axis into the complete physical boundary.
function _expand_bim_boundary(pts::BoundaryPoints{T},::BilliardGeometry.XAxisReflection,L::T) where {T<:Real}
    isempty(pts)&&throw(ArgumentError("Cannot expand an empty fundamental boundary"))
    sf=pts.s.-pts.s[1]
    sym=BilliardGeometry.XAxisReflection()
    xy=vcat(pts.xy,BilliardGeometry.apply_symmetry(sym,reverse(pts.xy)))
    normal=vcat(pts.normal,BilliardGeometry.apply_symmetry(sym,reverse(pts.normal)))
    s=vcat(sf,2L.-reverse(sf))
    ds=vcat(pts.ds,reverse(pts.ds))
    curvature=vcat(pts.curvature,reverse(pts.curvature))
    return BoundaryPoints(xy;normal=normal,s=s,ds=ds,curvature=curvature)
end

# Expand a half-boundary fundamental domain across the y-axis into the complete physical boundary.
function _expand_bim_boundary(pts::BoundaryPoints{T},::BilliardGeometry.YAxisReflection,L::T) where {T<:Real}
    isempty(pts)&&throw(ArgumentError("Cannot expand an empty fundamental boundary"))
    sf=pts.s.-pts.s[1]
    sym=BilliardGeometry.YAxisReflection()
    xy=vcat(pts.xy,BilliardGeometry.apply_symmetry(sym,reverse(pts.xy)))
    normal=vcat(pts.normal,BilliardGeometry.apply_symmetry(sym,reverse(pts.normal)))
    s=vcat(sf,2*L.-reverse(sf))
    ds=vcat(pts.ds,reverse(pts.ds))
    curvature=vcat(pts.curvature,reverse(pts.curvature))
    return BoundaryPoints(xy;normal=normal,s=s,ds=ds,curvature=curvature)
end

# Expand one rotational fundamental sector into the complete Cn-symmetric physical boundary.
function _expand_bim_boundary(pts::BoundaryPoints{T},sym::BilliardGeometry.NFoldRotation,L::T) where {T<:Real}
    isempty(pts)&&throw(ArgumentError("Cannot expand an empty fundamental boundary"))
    n=sym.order
    n>=2||throw(ArgumentError("NFoldRotation order must be at least 2; received $n"))
    sf=pts.s.-pts.s[1]
    M=length(pts)
    xy=Vector{SVector{2,T}}(undef,n*M)
    normal=Vector{SVector{2,T}}(undef,n*M)
    s=Vector{T}(undef,n*M)
    ds=Vector{T}(undef,n*M)
    curvature=Vector{T}(undef,n*M)
    @inbounds for l in 0:n-1
        off=l*M
        img=l==0 ? nothing : BilliardGeometry.NFoldRotation(n,l,sym.sector)
        for q in 1:M
            j=off+q
            xy[j]=l==0 ? pts.xy[q] : BilliardGeometry.apply_symmetry(img,pts.xy[q])
            normal[j]=l==0 ? pts.normal[q] : BilliardGeometry.apply_symmetry(img,pts.normal[q])
            s[j]=T(l)*L+sf[q]
            ds[j]=pts.ds[q]
            curvature[j]=pts.curvature[q]
        end
    end
    return BoundaryPoints(xy;normal=normal,s=s,ds=ds,curvature=curvature)
end

# Expand a rotational fundamental sector when the active Cn symmetry is stored as all nontrivial rotation images.
function _expand_bim_boundary(pts::BoundaryPoints{T},syms::AbstractVector{<:BilliardGeometry.NFoldRotation},L::T) where {T<:Real}
    isempty(syms)&&return pts
    n=syms[1].order
    sector=syms[1].sector
    all(sym.order==n&&sym.sector==sector for sym in syms)||throw(ArgumentError("All rotational symmetry images must have identical order and sector"))
    return _expand_bim_boundary(pts,BilliardGeometry.NFoldRotation(n,1,sector),L)
end

# Reject symmetry types for which no BIM fundamental-boundary expansion convention has been implemented.
function _expand_bim_boundary(pts::BoundaryPoints,sym::BilliardGeometry.AbsSymmetry,L)
    throw(ArgumentError("BoundaryIntegralMethod symmetry expansion is not implemented for $(typeof(sym))"))
end

"""
    evaluate_points(solver::BoundaryIntegralMethod{T,Sym},billiard::Bi,k) where {T<:Real,Sym,Bi<:BilliardGeometry.AbsBilliard} -> BoundaryPoints{T}

Construct the complete physical-boundary discretization used by the standard
boundary integral method.

Without symmetry, `billiard.full_boundary` is discretized directly.

When symmetry is active, the physical fundamental boundary returned by
`BilliardGeometry.get_boundary_curves(billiard)` is discretized first and then
expanded by geometric symmetry transformations into a complete,
canonically ordered physical boundary. Symmetry reduction of the resulting
operator is subsequently performed through `SymmetryOrbitMap`.

## Arguments
* `solver::BoundaryIntegralMethod{T}`: Boundary integral solver configuration.
* `billiard::Bi`: Billiard geometry.
* `k`: Wavenumber controlling the boundary discretization density.

## Returns
* `BoundaryPoints{T}`: Complete physical-boundary discretization.
"""
function evaluate_points(solver::BoundaryIntegralMethod{T},billiard::Bi,k) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    if isnothing(solver.symmetry)
        curves=billiard.full_boundary
        bs,samplers=_adjust_scaling_and_samplers(solver,length(curves))
        return _evaluate_bim_curves(solver,curves,bs,samplers,k)
    end
    curves=BilliardGeometry.get_boundary_curves(billiard)
    bs,samplers=adjust_scaling_and_samplers(solver,billiard)
    pts=_evaluate_bim_curves(solver,curves,bs,samplers,k)
    L=sum(T(crv.length) for crv in curves)
    return _expand_bim_boundary(pts,solver.symmetry,L)
end

################################################################################
########################## RAW DLP KERNEL ######################################
################################################################################

# Doubled 2D Helmholtz DLP scalar factors:
#
#   K   = (ik/2) H1(kr)
#   K'  = (ikr/2) H0(kr)
#   K'' = (i/2)[r H0(kr)-k r² H1(kr)].
#
# The geometric source-normal factor
#
#   c = n_y·(x-y)/r
#
# is applied separately.
@inline function _default_dlp_kernel_triplet(k::T,r::T) where {T<:Real}
    kr=k*r
    h0=Bessels.hankelh1(0,kr)
    h1=Bessels.hankelh1(1,kr)
    α=Complex{T}(zero(T),one(T)/2)
    K=α*k*h1
    dK=α*k*r*h0
    ddK=α*(r*h0-k*r*r*h1)
    return K,dK,ddK
end

@inline function _default_dlp_regular_entry(xi::T,yi::T,xj::T,yj::T,nxj::T,nyj::T,k::T) where {T<:Real}
    dx=xi-xj
    dy=yi-yj
    r=hypot(dx,dy)
    iszero(r)&&throw(ArgumentError("Regular DLP kernel received coincident target and source points"))
    c=(nxj*dx+nyj*dy)/r
    K,_,_=_default_dlp_kernel_triplet(k,r)
    return c*K
end

@inline function _default_dlp_regular_triplet(xi::T,yi::T,xj::T,yj::T,nxj::T,nyj::T,k::T) where {T<:Real}
    dx=xi-xj
    dy=yi-yj
    r=hypot(dx,dy)
    iszero(r)&&throw(ArgumentError("Regular DLP kernel received coincident target and source points"))
    c=(nxj*dx+nyj*dy)/r
    K,dK,ddK=_default_dlp_kernel_triplet(k,r)
    return c*K,c*dK,c*ddK
end

"""
    default_helmholtz_kernel_matrix(bp::BoundaryPoints{T},k::T;multithreaded::Bool=true) where {T<:Real} → Matrix{Complex{T}}

Assemble the raw doubled two-dimensional Helmholtz double-layer kernel matrix.

## Description
For distinct boundary points,

    K(x,y;k)=(ik/2)H₁⁽¹⁾(kr)
             [n_y⋅(x-y)/r].

The smooth diagonal limit is

    K(x,x;k)=-κ(x)/(2π).

No source quadrature weights or Fredholm identity shift are included.

## Arguments
* `bp::BoundaryPoints{T}`: Boundary discretization.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to use threaded kernel assembly.

## Returns
* `K::Matrix{Complex{T}}`: Raw DLP kernel matrix.
"""
function default_helmholtz_kernel_matrix(bp::BoundaryPoints{T},k::T;multithreaded::Bool=true) where {T<:Real}
    N=length(bp)
    M=Matrix{Complex{T}}(undef,N,N)
    xy=bp.xy
    nrm=bp.normal
    κ=bp.curvature
    @use_threads multithreading=multithreaded for i in 1:N
        pi=xy[i]
        ni=nrm[i]
        @inbounds for j in 1:i
            if i==j
                M[i,i]=-Complex{T}(κ[i]/TWO_PI,zero(T))
            else
                pj=xy[j]
                nj=nrm[j]
                dx=pi[1]-pj[1]
                dy=pi[2]-pj[2]
                r=hypot(dx,dy)
                K,_,_=_default_dlp_kernel_triplet(k,r)
                cij=(nj[1]*dx+nj[2]*dy)/r
                cji=(ni[1]*(-dx)+ni[2]*(-dy))/r
                M[i,j]=cij*K
                M[j,i]=cji*K
            end
        end
    end
    filter_matrix!(M)
    return M
end

"""
    default_helmholtz_kernel_derivative_matrix(bp::BoundaryPoints{T},k::T;multithreaded::Bool=true) where {T<:Real} → Matrix{Complex{T}}

Assemble the first wavenumber derivative of the raw DLP kernel.

For distinct points,

    ∂K/∂k=(ikr/2)H₀⁽¹⁾(kr)
           [n_y⋅(x-y)/r].

The diagonal derivative vanishes.

## Arguments
* `bp::BoundaryPoints{T}`: Boundary discretization.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to use threaded assembly.

## Returns
* `dK::Matrix{Complex{T}}`: First derivative of the raw DLP kernel.
"""
function default_helmholtz_kernel_derivative_matrix(bp::BoundaryPoints{T},k::T;multithreaded::Bool=true) where {T<:Real}
    N=length(bp)
    M=zeros(Complex{T},N,N)
    xy=bp.xy
    nrm=bp.normal
    @use_threads multithreading=multithreaded for i in 1:N
        pi=xy[i]
        ni=nrm[i]
        @inbounds for j in 1:i-1
            pj=xy[j]
            nj=nrm[j]
            dx=pi[1]-pj[1]
            dy=pi[2]-pj[2]
            r=hypot(dx,dy)
            _,dK,_=_default_dlp_kernel_triplet(k,r)
            cij=(nj[1]*dx+nj[2]*dy)/r
            cji=(ni[1]*(-dx)+ni[2]*(-dy))/r
            M[i,j]=cij*dK
            M[j,i]=cji*dK
        end
    end
    filter_matrix!(M)
    return M
end

"""
    default_helmholtz_kernel_second_derivative_matrix(bp::BoundaryPoints{T},k::T;multithreaded::Bool=true) where {T<:Real} → Matrix{Complex{T}}

Assemble the second wavenumber derivative of the raw DLP kernel.
For distinct points,

    ∂²K/∂k² =
        (i/2)[rH₀⁽¹⁾(kr)-kr²H₁⁽¹⁾(kr)]
        [n_y⋅(x-y)/r].

The diagonal second derivative vanishes.

## Arguments
* `bp::BoundaryPoints{T}`: Boundary discretization.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to use threaded assembly.

## Returns
* `ddK::Matrix{Complex{T}}`: Second derivative of the raw DLP kernel.
"""
function default_helmholtz_kernel_second_derivative_matrix(bp::BoundaryPoints{T},k::T;multithreaded::Bool=true) where {T<:Real}
    N=length(bp)
    M=zeros(Complex{T},N,N)
    xy=bp.xy
    nrm=bp.normal
    @use_threads multithreading=multithreaded for i in 1:N
        pi=xy[i]
        ni=nrm[i]
        @inbounds for j in 1:i-1
            pj=xy[j]
            nj=nrm[j]
            dx=pi[1]-pj[1]
            dy=pi[2]-pj[2]
            r=hypot(dx,dy)
            _,_,ddK=_default_dlp_kernel_triplet(k,r)
            cij=(nj[1]*dx+nj[2]*dy)/r
            cji=(ni[1]*(-dx)+ni[2]*(-dy))/r
            M[i,j]=cij*ddK
            M[j,i]=cji*ddK
        end
    end
    filter_matrix!(M)
    return M
end

################################################################################
############################ FULL-BOUNDARY ASSEMBLY ############################
################################################################################

"""
    compute_kernel_matrix!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},k::T;multithreaded::Bool=true) where {T<:Real} → K

Assemble the full raw DLP kernel matrix in place.

## Arguments
* `K::AbstractMatrix{Complex{T}}`: Preallocated full destination matrix.
* `bp::BoundaryPoints{T}`: Full boundary discretization.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to use threaded assembly.

## Returns
* `K::AbstractMatrix{Complex{T}}`: Raw full DLP kernel matrix.
"""
function compute_kernel_matrix!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},k::T;multithreaded::Bool=true) where {T<:Real}
    N=length(bp)
    @assert size(K)==(N,N)
    fill!(K,zero(Complex{T}))
    xy=bp.xy
    nrm=bp.normal
    κ=bp.curvature
    @use_threads multithreading=multithreaded for i in 1:N
        pi=xy[i]
        ni=nrm[i]
        @inbounds begin
            K[i,i]=-Complex{T}(κ[i]/TWO_PI,zero(T))
            for j in 1:i-1
                pj=xy[j]
                nj=nrm[j]
                dx=pi[1]-pj[1]
                dy=pi[2]-pj[2]
                r=hypot(dx,dy)
                hK,_,_=_default_dlp_kernel_triplet(k,r)
                cij=(nj[1]*dx+nj[2]*dy)/r
                cji=(ni[1]*(-dx)+ni[2]*(-dy))/r
                K[i,j]=cij*hK
                K[j,i]=cji*hK
            end
        end
    end
    return K
end

"""
    compute_kernel_matrix_with_derivatives!(K::AbstractMatrix{Complex{T}},dK::AbstractMatrix{Complex{T}},ddK::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},k::T;multithreaded::Bool=true) where {T<:Real} → K,dK,ddK

Assemble the full raw DLP kernel and its first two wavenumber derivatives.

## Arguments
* `K::AbstractMatrix{Complex{T}}`: Destination matrix for the raw DLP kernel.
* `dK::AbstractMatrix{Complex{T}}`: Destination matrix for the first derivative.
* `ddK::AbstractMatrix{Complex{T}}`: Destination matrix for the second derivative.
* `bp::BoundaryPoints{T}`: Full boundary discretization.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to use threaded assembly.

## Returns
* `K::AbstractMatrix{Complex{T}}`: Raw DLP kernel matrix.
* `dK::AbstractMatrix{Complex{T}}`: First derivative.
* `ddK::AbstractMatrix{Complex{T}}`: Second derivative.
"""
function compute_kernel_matrix_with_derivatives!(K::AbstractMatrix{Complex{T}},dK::AbstractMatrix{Complex{T}},ddK::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},k::T;multithreaded::Bool=true) where {T<:Real}
    N=length(bp)
    @assert size(K)==(N,N)
    @assert size(dK)==(N,N)
    @assert size(ddK)==(N,N)
    fill!(K,zero(Complex{T}))
    fill!(dK,zero(Complex{T}))
    fill!(ddK,zero(Complex{T}))
    xy=bp.xy
    nrm=bp.normal
    κ=bp.curvature
    @use_threads multithreading=multithreaded for i in 1:N
        pi=xy[i]
        ni=nrm[i]
        @inbounds begin
            K[i,i]=-Complex{T}(κ[i]/TWO_PI,zero(T))
            for j in 1:i-1
                pj=xy[j]
                nj=nrm[j]
                dx=pi[1]-pj[1]
                dy=pi[2]-pj[2]
                r=hypot(dx,dy)
                hK,hdK,hddK=_default_dlp_kernel_triplet(k,r)
                cij=(nj[1]*dx+nj[2]*dy)/r
                cji=(ni[1]*(-dx)+ni[2]*(-dy))/r
                K[i,j]=cij*hK
                dK[i,j]=cij*hdK
                ddK[i,j]=cij*hddK
                K[j,i]=cji*hK
                dK[j,i]=cji*hdK
                ddK[j,i]=cji*hddK
            end
        end
    end
    return K,dK,ddK
end

################################################################################
######################## SYMMETRY-REDUCED ASSEMBLY #############################
################################################################################

"""
    compute_kernel_matrix!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},k::T;multithreaded::Bool=true) where {T<:Real} → K

Assemble the symmetry-reduced raw DLP kernel in exact boundary-index space.

## Description
For fundamental source index `b`, let `j=Ifund[b]`. The reduced raw kernel is

    K_red[a,b]
      = Σ_g χ(g) K(i,g·j) ds[g·j]/ds[j],

where `i=Ifund[a]`.

The weight ratio is included so that subsequent multiplication of reduced
column `b` by `ds[j]` reproduces the complete full-boundary Nyström source sum,

    Σ_g χ(g)K(i,g·j)ds[g·j].

## Arguments
* `K::AbstractMatrix{Complex{T}}`: Preallocated reduced destination matrix.
* `bp::BoundaryPoints{T}`: Full boundary discretization.
* `orbits::SymmetryOrbitMap{T}`: Exact symmetry-orbit map.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to use threaded assembly.

## Returns
* `K::AbstractMatrix{Complex{T}}`: Symmetry-reduced raw DLP kernel.
"""
function compute_kernel_matrix!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},k::T;multithreaded::Bool=true) where {T<:Real}
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    ng=orbit_size(orbits)
    @assert size(K)==(m,m)
    fill!(K,zero(Complex{T}))
    xy=bp.xy
    nrm=bp.normal
    κ=bp.curvature
    ds=bp.ds
    @use_threads multithreading=multithreaded for b in 1:m
        j=Ifund[b]
        wj=ds[j]
        @inbounds for a in 1:m
            i=Ifund[a]
            pi=xy[i]
            val=zero(Complex{T})
            for l in 1:ng
                q=orbits.fund_to_full[l,b]
                scale=orbits.fund_to_scale[l,b]
                weight_ratio=ds[q]/wj
                if l==1&&i==j
                    val+=scale*(-Complex{T}(κ[i]/TWO_PI,zero(T)))*weight_ratio
                else
                    pq=xy[q]
                    nq=nrm[q]
                    dx=pi[1]-pq[1]
                    dy=pi[2]-pq[2]
                    r=hypot(dx,dy)
                    iszero(r)&&throw(ArgumentError("A nonidentity symmetry image coincides with a reduced target node"))
                    hK,_,_=_default_dlp_kernel_triplet(k,r)
                    c=(nq[1]*dx+nq[2]*dy)/r
                    val+=scale*c*hK*weight_ratio
                end
            end
            K[a,b]=val
        end
    end
    return K
end

"""
    compute_kernel_matrix_with_derivatives!(K::AbstractMatrix{Complex{T}},dK::AbstractMatrix{Complex{T}},ddK::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},k::T;multithreaded::Bool=true) where {T<:Real} → K,dK,ddK

Assemble the symmetry-reduced raw DLP kernel and its first two wavenumber
derivatives.

## Arguments
* `K::AbstractMatrix{Complex{T}}`: Destination matrix for the reduced raw kernel.
* `dK::AbstractMatrix{Complex{T}}`: Destination matrix for the first derivative.
* `ddK::AbstractMatrix{Complex{T}}`: Destination matrix for the second derivative.
* `bp::BoundaryPoints{T}`: Full boundary discretization.
* `orbits::SymmetryOrbitMap{T}`: Exact symmetry-orbit map.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to use threaded assembly.

## Returns
* `K::AbstractMatrix{Complex{T}}`: Reduced raw DLP kernel.
* `dK::AbstractMatrix{Complex{T}}`: First derivative.
* `ddK::AbstractMatrix{Complex{T}}`: Second derivative.
"""
function compute_kernel_matrix_with_derivatives!(K::AbstractMatrix{Complex{T}},dK::AbstractMatrix{Complex{T}},ddK::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},k::T;multithreaded::Bool=true) where {T<:Real}
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    ng=orbit_size(orbits)
    @assert size(K)==(m,m)
    @assert size(dK)==(m,m)
    @assert size(ddK)==(m,m)
    fill!(K,zero(Complex{T}))
    fill!(dK,zero(Complex{T}))
    fill!(ddK,zero(Complex{T}))
    xy=bp.xy
    nrm=bp.normal
    κ=bp.curvature
    ds=bp.ds
    @use_threads multithreading=multithreaded for b in 1:m
        j=Ifund[b]
        wj=ds[j]
        @inbounds for a in 1:m
            i=Ifund[a]
            pi=xy[i]
            v=zero(Complex{T})
            v1=zero(Complex{T})
            v2=zero(Complex{T})
            for l in 1:ng
                q=orbits.fund_to_full[l,b]
                scale=orbits.fund_to_scale[l,b]
                weight_ratio=ds[q]/wj
                if l==1&&i==j
                    v+=scale*(-Complex{T}(κ[i]/TWO_PI,zero(T)))*weight_ratio
                else
                    pq=xy[q]
                    nq=nrm[q]
                    dx=pi[1]-pq[1]
                    dy=pi[2]-pq[2]
                    r=hypot(dx,dy)
                    iszero(r)&&throw(ArgumentError("A nonidentity symmetry image coincides with a reduced target node"))
                    hK,hdK,hddK=_default_dlp_kernel_triplet(k,r)
                    c=(nq[1]*dx+nq[2]*dy)/r
                    s=scale*c*weight_ratio
                    v+=s*hK
                    v1+=s*hdK
                    v2+=s*hddK
                end
            end
            K[a,b]=v
            dK[a,b]=v1
            ddK[a,b]=v2
        end
    end
    return K,dK,ddK
end

################################################################################
############################ FREDHOLM MATRICES #################################
################################################################################

function fredholm_matrix!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},::Nothing,k::T;multithreaded::Bool=true) where {T<:Real}
    compute_kernel_matrix!(K,bp,k;multithreaded=multithreaded)
    ds=bp.ds
    @inbounds for j in eachindex(ds)
        @views K[:,j].*=ds[j]
    end
    K.*=-one(T)
    @inbounds for i in axes(K,1)
        K[i,i]+=one(Complex{T})
    end
    return K
end

function fredholm_matrix!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},k::T;multithreaded::Bool=true) where {T<:Real}
    compute_kernel_matrix!(K,bp,orbits,k;multithreaded=multithreaded)
    Ifund=orbits.Ifund
    @inbounds for b in 1:fundamental_size(orbits)
        @views K[:,b].*=bp.ds[Ifund[b]]
    end
    K.*=-one(T)
    @inbounds for i in axes(K,1)
        K[i,i]+=one(Complex{T})
    end
    return K
end

"""
    fredholm_matrix!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},symmetry,k::T;multithreaded::Bool=true) where {T<:Real} → K

Assemble the BIM Fredholm matrix

    A(k)=I-K(k)W.

When symmetry is active, the full source sum is first folded through the exact
symmetry-orbit map.

## Arguments
* `K::AbstractMatrix{Complex{T}}`: Preallocated full or reduced destination matrix.
* `bp::BoundaryPoints{T}`: Full boundary discretization.
* `symmetry`: Optional symmetry descriptor or `nothing`.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to use threaded assembly.

## Returns
* `K::AbstractMatrix{Complex{T}}`: Fredholm matrix.
"""
function fredholm_matrix!(K::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},symmetry,k::T;multithreaded::Bool=true) where {T<:Real}
    if isnothing(symmetry)
        return fredholm_matrix!(K,bp,nothing,k;multithreaded=multithreaded)
    end
    orbits=symmetry_index_orbits(T,bp,symmetry)
    return fredholm_matrix!(K,bp,orbits,k;multithreaded=multithreaded)
end

function fredholm_matrix_with_derivatives!(K::AbstractMatrix{Complex{T}},dK::AbstractMatrix{Complex{T}},ddK::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},::Nothing,k::T;multithreaded::Bool=true) where {T<:Real}
    compute_kernel_matrix_with_derivatives!(K,dK,ddK,bp,k;multithreaded=multithreaded)
    ds=bp.ds
    @inbounds for j in eachindex(ds)
        @views K[:,j].*=ds[j]
        @views dK[:,j].*=ds[j]
        @views ddK[:,j].*=ds[j]
    end
    K.*=-one(T)
    dK.*=-one(T)
    ddK.*=-one(T)
    @inbounds for i in axes(K,1)
        K[i,i]+=one(Complex{T})
    end
    return K,dK,ddK
end

function fredholm_matrix_with_derivatives!(K::AbstractMatrix{Complex{T}},dK::AbstractMatrix{Complex{T}},ddK::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},k::T;multithreaded::Bool=true) where {T<:Real}
    compute_kernel_matrix_with_derivatives!(K,dK,ddK,bp,orbits,k;multithreaded=multithreaded)
    Ifund=orbits.Ifund
    @inbounds for b in 1:fundamental_size(orbits)
        w=bp.ds[Ifund[b]]
        @views K[:,b].*=w
        @views dK[:,b].*=w
        @views ddK[:,b].*=w
    end
    K.*=-one(T)
    dK.*=-one(T)
    ddK.*=-one(T)
    @inbounds for i in axes(K,1)
        K[i,i]+=one(Complex{T})
    end
    return K,dK,ddK
end

"""
    fredholm_matrix_with_derivatives!(K::AbstractMatrix{Complex{T}},dK::AbstractMatrix{Complex{T}},ddK::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},symmetry,k::T;multithreaded::Bool=true) where {T<:Real} → K,dK,ddK

Assemble

    A(k)=I-K(k)W,
    A'(k)=-K'(k)W,
    A''(k)=-K''(k)W.

## Arguments
* `K::AbstractMatrix{Complex{T}}`: Destination matrix for the Fredholm operator.
* `dK::AbstractMatrix{Complex{T}}`: Destination matrix for the first derivative.
* `ddK::AbstractMatrix{Complex{T}}`: Destination matrix for the second derivative.
* `bp::BoundaryPoints{T}`: Full boundary discretization.
* `symmetry`: Optional symmetry descriptor or `nothing`.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to use threaded assembly.

## Returns
* `K::AbstractMatrix{Complex{T}}`: Fredholm matrix.
* `dK::AbstractMatrix{Complex{T}}`: First derivative.
* `ddK::AbstractMatrix{Complex{T}}`: Second derivative.
"""
function fredholm_matrix_with_derivatives!(K::AbstractMatrix{Complex{T}},dK::AbstractMatrix{Complex{T}},ddK::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},symmetry,k::T;multithreaded::Bool=true) where {T<:Real}
    if isnothing(symmetry)
        return fredholm_matrix_with_derivatives!(K,dK,ddK,bp,nothing,k;multithreaded=multithreaded)
    end
    orbits=symmetry_index_orbits(T,bp,symmetry)
    return fredholm_matrix_with_derivatives!(K,dK,ddK,bp,orbits,k;multithreaded=multithreaded)
end

function compute_kernel_matrix(bp::BoundaryPoints{T},k::T;multithreaded::Bool=true) where {T<:Real}
    N=length(bp)
    K=Matrix{Complex{T}}(undef,N,N)
    compute_kernel_matrix!(K,bp,k;multithreaded=multithreaded)
    return K
end

function compute_kernel_matrix(bp::BoundaryPoints{T},symmetry,k::T;multithreaded::Bool=true) where {T<:Real}
    if isnothing(symmetry)
        return compute_kernel_matrix(bp,k;multithreaded=multithreaded)
    end
    orbits=symmetry_index_orbits(T,bp,symmetry)
    m=fundamental_size(orbits)
    K=Matrix{Complex{T}}(undef,m,m)
    compute_kernel_matrix!(K,bp,orbits,k;multithreaded=multithreaded)
    return K
end

function fredholm_matrix(bp::BoundaryPoints{T},symmetry,k::T;multithreaded::Bool=true) where {T<:Real}
    orbits=isnothing(symmetry) ? nothing : symmetry_index_orbits(T,bp,symmetry)
    n=_dlp_matrix_dim(bp,orbits)
    K=Matrix{Complex{T}}(undef,n,n)
    fredholm_matrix!(K,bp,orbits,k;multithreaded=multithreaded)
    return K
end

function fredholm_matrix_with_derivatives(bp::BoundaryPoints{T},symmetry,k::T;multithreaded::Bool=true) where {T<:Real}
    orbits=isnothing(symmetry) ? nothing : symmetry_index_orbits(T,bp,symmetry)
    n=_dlp_matrix_dim(bp,orbits)
    K=Matrix{Complex{T}}(undef,n,n)
    dK=Matrix{Complex{T}}(undef,n,n)
    ddK=Matrix{Complex{T}}(undef,n,n)
    fredholm_matrix_with_derivatives!(K,dK,ddK,bp,orbits,k;multithreaded=multithreaded)
    return K,dK,ddK
end

################################################################################
############################ NEEDED FOR HUSIMIS ################################
################################################################################

function adjoint_fredholm_matrix!(A::AbstractMatrix{Complex{T}},D::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},::Nothing,k::T;multithreaded::Bool=true) where {T<:Real}
    N=length(bp)
    @assert size(A)==(N,N)
    @assert size(D)==(N,N)
    compute_kernel_matrix!(D,bp,k;multithreaded=multithreaded)
    ds=bp.ds
    @inbounds for j in 1:N
        @views D[:,j].*=ds[j]
    end
    @inbounds for i in 1:N,j in 1:N
        A[i,j]=-D[j,i]*ds[j]/ds[i]
    end
    @inbounds for i in 1:N
        A[i,i]+=one(Complex{T})
    end
    return A
end

function adjoint_fredholm_matrix!(A::AbstractMatrix{Complex{T}},D::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},orbits::SymmetryOrbitMap{T},k::T;multithreaded::Bool=true) where {T<:Real}
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    @assert size(A)==(m,m)
    @assert size(D)==(m,m)
    compute_kernel_matrix!(D,bp,orbits,k;multithreaded=multithreaded)
    @inbounds for b in 1:m
        @views D[:,b].*=bp.ds[Ifund[b]]
    end
    @inbounds for b in 1:m,a in 1:m
        wa=bp.ds[Ifund[a]]
        wb=bp.ds[Ifund[b]]
        A[a,b]=-D[b,a]*wb/wa
    end
    @inbounds for a in 1:m
        A[a,a]+=one(Complex{T})
    end
    return A
end

"""
    adjoint_fredholm_matrix!(
        A::AbstractMatrix{Complex{T}},
        D::AbstractMatrix{Complex{T}},
        bp::BoundaryPoints{T},
        symmetry,
        k::T;
        multithreaded::Bool=true,
    ) where {T<:Real} → A

Assemble the weighted-transpose DLP Fredholm matrix used for boundary-function
and Husimi postprocessing.

## Description
For the quadrature-weighted DLP matrix `D`, the discrete formal transpose is

    K'=W⁻¹DᵀW,

and

    A=I-K'.

This is a bilinear weighted transpose; no complex conjugation is applied.

## Arguments
* `A::AbstractMatrix{Complex{T}}`: Preallocated destination matrix.
* `D::AbstractMatrix{Complex{T}}`: Preallocated DLP workspace matrix.
* `bp::BoundaryPoints{T}`: Full boundary discretization.
* `symmetry`: Optional symmetry descriptor or `nothing`.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to use threaded assembly.

## Returns
* `A::AbstractMatrix{Complex{T}}`: Weighted-transpose Fredholm matrix.
"""
function adjoint_fredholm_matrix!(A::AbstractMatrix{Complex{T}},D::AbstractMatrix{Complex{T}},bp::BoundaryPoints{T},symmetry,k::T;multithreaded::Bool=true) where {T<:Real}
    if isnothing(symmetry)
        return adjoint_fredholm_matrix!(A,D,bp,nothing,k;multithreaded=multithreaded)
    end
    orbits=symmetry_index_orbits(T,bp,symmetry)
    return adjoint_fredholm_matrix!(A,D,bp,orbits,k;multithreaded=multithreaded)
end

"""
    smallest_nullvec_krylov!(
        A::AbstractMatrix{Complex{T}};
        nev::Int=1,
        tol=1e-12,
        maxiter::Int=2000,
        krylovdim::Int=40,
    ) where {T<:Real} → σ,u,info

Compute a near-null eigenvector by applying a Krylov eigensolver to `A⁻¹`.

## Description
The eigenvalue of `A⁻¹` with largest magnitude corresponds to an eigenvalue of
`A` nearest zero. If this inverse eigenvalue is `μ`, the returned scalar is

    σ=1/|μ|.

Thus `σ` is a near-zero eigenvalue-magnitude proxy; it is not, in general, a
singular value.

## Arguments
* `A::AbstractMatrix{Complex{T}}`: Matrix whose near-null vector is sought.

## Keyword Arguments
* `nev::Int`: Number of inverse eigenpairs requested.
* `tol`: Krylov convergence tolerance.
* `maxiter::Int`: Maximum number of Krylov iterations.
* `krylovdim::Int`: Krylov subspace dimension.

## Returns
* `σ`: Near-zero eigenvalue-magnitude proxy.
* `u`: Normalized approximate null vector.
* `info`: Krylov solver diagnostic information.
"""
function smallest_nullvec_krylov!(A::AbstractMatrix{Complex{T}};nev::Int=1,tol=1e-12,maxiter::Int=2000,krylovdim::Int=40) where {T<:Real}
    n=size(A,1)
    F=lu!(A)
    function op!(y,x)
        copyto!(y,x)
        ldiv!(F,y)
        return y
    end
    C=LinearMaps.LinearMap{eltype(A)}(op!,n,n;ismutating=true)
    μs,vecs,info=eigsolve(C,n,nev,:LM;tol=tol,maxiter=maxiter,krylovdim=krylovdim)
    μ=μs[1]
    u=vecs[1]./norm(vecs[1])
    σ=inv(abs(μ))
    return σ,u,info
end

################################################################################
############################ SOLVER INTERFACE ##################################
################################################################################

"""
    construct_matrices!(
        solver::BoundaryIntegralMethod,
        basis::Ba,
        A::AbstractMatrix{Complex{T}},
        pts::BoundaryPoints{T},
        k::T;
        multithreaded::Bool=true,
    ) where {Ba<:AbstractHankelBasis,T<:Real} → A

Assemble the BIM Fredholm matrix into a preallocated matrix.

## Arguments
* `solver::BoundaryIntegralMethod`: Boundary integral solver.
* `basis::Ba`: Compatibility basis placeholder.
* `A::AbstractMatrix{Complex{T}}`: Preallocated destination matrix.
* `pts::BoundaryPoints{T}`: Full boundary discretization.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to use threaded assembly.

## Returns
* `A::AbstractMatrix{Complex{T}}`: Fredholm matrix.
"""
function construct_matrices!(solver::BoundaryIntegralMethod,basis::Ba,A::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},k::T;multithreaded::Bool=true) where {Ba<:AbstractHankelBasis,T<:Real}
    @blas_1 fredholm_matrix!(A,pts,solver.symmetry,k;multithreaded=multithreaded)
    return A
end

"""
    construct_matrices!(
        solver::BoundaryIntegralMethod,
        basis::Ba,
        A::AbstractMatrix{Complex{T}},
        dA::AbstractMatrix{Complex{T}},
        ddA::AbstractMatrix{Complex{T}},
        pts::BoundaryPoints{T},
        k::T;
        multithreaded::Bool=true,
    ) where {Ba<:AbstractHankelBasis,T<:Real} → A,dA,ddA

Assemble the BIM Fredholm matrix and its first two wavenumber derivatives.

## Arguments
* `solver::BoundaryIntegralMethod`: Boundary integral solver.
* `basis::Ba`: Compatibility basis placeholder.
* `A::AbstractMatrix{Complex{T}}`: Destination matrix for the Fredholm operator.
* `dA::AbstractMatrix{Complex{T}}`: Destination matrix for the first derivative.
* `ddA::AbstractMatrix{Complex{T}}`: Destination matrix for the second derivative.
* `pts::BoundaryPoints{T}`: Full boundary discretization.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to use threaded assembly.

## Returns
* `A::AbstractMatrix{Complex{T}}`: Fredholm matrix.
* `dA::AbstractMatrix{Complex{T}}`: First derivative.
* `ddA::AbstractMatrix{Complex{T}}`: Second derivative.
"""
function construct_matrices!(solver::BoundaryIntegralMethod,basis::Ba,A::AbstractMatrix{Complex{T}},dA::AbstractMatrix{Complex{T}},ddA::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},k::T;multithreaded::Bool=true) where {Ba<:AbstractHankelBasis,T<:Real}
    @blas_1 fredholm_matrix_with_derivatives!(A,dA,ddA,pts,solver.symmetry,k;multithreaded=multithreaded)
    return A,dA,ddA
end

function construct_matrices(solver::BoundaryIntegralMethod,basis::Ba,pts::BoundaryPoints{T},k::T,A::AbstractMatrix{Complex{T}},dA::AbstractMatrix{Complex{T}},ddA::AbstractMatrix{Complex{T}};multithreaded::Bool=true) where {Ba<:AbstractHankelBasis,T<:Real}
    construct_matrices!(solver,basis,A,dA,ddA,pts,k;multithreaded=multithreaded)
    return A,dA,ddA
end

function construct_matrices(solver::BoundaryIntegralMethod,basis::Ba,pts::BoundaryPoints{T},k::T;multithreaded::Bool=true) where {Ba<:AbstractHankelBasis,T<:Real}
    orbits=_dlp_symmetry_orbits(solver,pts)
    N=_dlp_matrix_dim(pts,orbits)
    A=Matrix{Complex{T}}(undef,N,N)
    dA=Matrix{Complex{T}}(undef,N,N)
    ddA=Matrix{Complex{T}}(undef,N,N)
    if isnothing(orbits)
        @blas_1 fredholm_matrix_with_derivatives!(A,dA,ddA,pts,nothing,k;multithreaded=multithreaded)
    else
        @blas_1 fredholm_matrix_with_derivatives!(A,dA,ddA,pts,orbits,k;multithreaded=multithreaded)
    end
    return A,dA,ddA
end

"""
    solve(
        solver::BoundaryIntegralMethod,
        basis::Ba,
        pts::BoundaryPoints{T},
        k::T;
        multithreaded::Bool=true,
        use_krylov::Bool=true,
        which::Symbol=:det_argmin,
    ) where {Ba<:AbstractHankelBasis,T<:Real}

Evaluate the selected scalar spectral diagnostic of the BIM Fredholm matrix.

## Arguments
* `solver::BoundaryIntegralMethod`: Boundary integral solver.
* `basis::Ba`: Compatibility basis placeholder.
* `pts::BoundaryPoints{T}`: Full boundary discretization.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to use threaded assembly.
* `use_krylov::Bool`: Whether the scalar backend may use its Krylov pathway.
* `which::Symbol`: Scalar diagnostic to evaluate.

## Returns
* Scalar spectral diagnostic selected by `which`.
"""
function solve(solver::BoundaryIntegralMethod,basis::Ba,pts::BoundaryPoints{T},k::T;multithreaded::Bool=true,use_krylov::Bool=true,which::Symbol=:det_argmin) where {Ba<:AbstractHankelBasis,T<:Real}
    orbits=_dlp_symmetry_orbits(solver,pts)
    N=_dlp_matrix_dim(pts,orbits)
    A=Matrix{Complex{T}}(undef,N,N)
    if isnothing(orbits)
        @blas_1 fredholm_matrix!(A,pts,nothing,k;multithreaded=multithreaded)
    else
        @blas_1 fredholm_matrix!(A,pts,orbits,k;multithreaded=multithreaded)
    end
    @svd_or_det_solve A use_krylov which MAX_BLAS_THREADS
end

function solve(solver::BoundaryIntegralMethod,basis::Ba,A::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},k::T;multithreaded::Bool=true,use_krylov::Bool=true,which::Symbol=:det_argmin) where {Ba<:AbstractHankelBasis,T<:Real}
    @blas_1 construct_matrices!(solver,basis,A,pts,k;multithreaded=multithreaded)
    @svd_or_det_solve A use_krylov which MAX_BLAS_THREADS
end

"""
    solve_vect(solver::BoundaryIntegralMethod,basis::Ba,pts::BoundaryPoints{T},k::T;multithreaded::Bool=true,tol=1e-12,maxiter::Int=2000,krylovdim::Int=40) where {Ba<:AbstractHankelBasis,T<:Real} → σ,u

Compute a near-null vector of the weighted-transpose BIM Fredholm matrix.

## Arguments
* `solver::BoundaryIntegralMethod`: Boundary integral solver.
* `basis::Ba`: Compatibility basis placeholder.
* `pts::BoundaryPoints{T}`: Full boundary discretization.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to use threaded assembly.
* `tol`: Krylov convergence tolerance.
* `maxiter::Int`: Maximum number of Krylov iterations.
* `krylovdim::Int`: Krylov subspace dimension.

## Returns
* `σ`: Near-zero eigenvalue-magnitude proxy.
* `u`: Associated normalized near-null boundary vector.
"""
function solve_vect(solver::BoundaryIntegralMethod,basis::Ba,pts::BoundaryPoints{T},k::T;multithreaded::Bool=true,tol=1e-12,maxiter::Int=2000,krylovdim::Int=40) where {Ba<:AbstractHankelBasis,T<:Real}
    orbits=_dlp_symmetry_orbits(solver,pts)
    N=_dlp_matrix_dim(pts,orbits)
    A=Matrix{Complex{T}}(undef,N,N)
    D=similar(A)
    if isnothing(orbits)
        @blas_1 adjoint_fredholm_matrix!(A,D,pts,nothing,k;multithreaded=multithreaded)
    else
        @blas_1 adjoint_fredholm_matrix!(A,D,pts,orbits,k;multithreaded=multithreaded)
    end
    σ,u,_=smallest_nullvec_krylov!(A;nev=1,tol=tol,maxiter=maxiter,krylovdim=krylovdim)
    return σ,u
end

function solve_vect(solver::BoundaryIntegralMethod,basis::Ba,A::AbstractMatrix{Complex{T}},pts::BoundaryPoints{T},k::T;multithreaded::Bool=true,tol=1e-12,maxiter::Int=2000,krylovdim::Int=40) where {Ba<:AbstractHankelBasis,T<:Real}
    D=similar(A)
    @blas_1 adjoint_fredholm_matrix!(A,D,pts,solver.symmetry,k;multithreaded=multithreaded)
    σ,u,_=smallest_nullvec_krylov!(A;nev=1,tol=tol,maxiter=maxiter,krylovdim=krylovdim)
    return σ,u
end

# INTERNAL - only for testing performance of the solve workflow.
function solve_INFO(solver::BoundaryIntegralMethod,basis::Ba,pts::BoundaryPoints{T},k::T;multithreaded::Bool=true,use_krylov::Bool=true,which::Symbol=:det_argmin) where {Ba<:AbstractHankelBasis,T<:Real}
    orbits=_dlp_symmetry_orbits(solver,pts)
    N=_dlp_matrix_dim(pts,orbits)
    A=Matrix{Complex{T}}(undef,N,N)
    s_constr=time()
    @info "constructing Fredholm matrix A..."
    if isnothing(orbits)
        @blas_1 fredholm_matrix!(A,pts,nothing,k;multithreaded=multithreaded)
    else
        @blas_1 fredholm_matrix!(A,pts,orbits,k;multithreaded=multithreaded)
    end
    @info "Condition number of A for svd: $(cond(A))"
    e_constr=time()
    s_svd=time()
    mu=@svd_or_det_solve A use_krylov which MAX_BLAS_THREADS
    e_svd=time()
    total_time=(e_svd-s_svd)+(e_constr-s_constr)
    @info "Total solve time for test k: $(total_time)"
    println("%%%%% SUMMARY %%%%%")
    println("Percentage of total time (most relevant ones): ")
    println("Fredholm matrix A construction: $(100*(e_constr-s_constr)/total_time) %")
    println("SVD: $(100*(e_svd-s_svd)/total_time) %")
    println("%%%%%%%%%%%%%%%%%%%")
    return mu
end