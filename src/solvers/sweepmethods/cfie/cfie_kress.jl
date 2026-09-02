# Useful reading:
# - https://github.com/ahbarnett/mpspack - by Alex Barnett & Timo Betcke (MATLAB)
# - Kress, R., Boundary integral equations in time-harmonic acoustic scattering. Mathematics Comput. Modelling Vol 15, pp. 229-243). Pergamon Press, 1991, GB.
# - Barnett, A. H., & Betcke, T. (2007). Stability and convergence of the method of fundamental solutions for Helmholtz problems on analytic domains. Journal of Computational Physics, 227(14), 7003-7026.
# - Zhao, L., & Barnett, A. (2015). Robust and efficient solution of the drum problem via Nyström approximation of the Fredholm determinant. SIAM Journal on Numerical Analysis, Stable URL: https://www.jstor.org/stable/24512689

const euler_over_pi=MathConstants.eulergamma/pi
abstract type CFIE<:SweepSolver end
############################
#### CONSTRUCTOR KRESS ######
############################

"""
    CFIE_kress{T,Bi,Sym} <: CFIE

Combined-field integral equation solver using Kress periodic logarithmic
splitting on smooth closed boundary components.

## Description
`CFIE_kress` implements the smooth-boundary Kress version of the combined-field
Fredholm formulation. Each connected boundary component must be represented by
one smooth closed periodic curve.

The discretized operator has the form

    A(k)=I-(D(k)+ikS(k)),

where `D(k)` is the Helmholtz double-layer operator and `S(k)` is the
single-layer operator.

For same-component interactions, both kernels are decomposed into

    logarithmic coefficient × log(4sin²((t-s)/2))
    + smooth remainder.

The universal periodic logarithmic contribution is integrated using the Kress
correction matrix, while the smooth remainder is evaluated with periodic
trapezoidal quadrature.

The combined-field term removes the spurious interior-nullspace problem that can
occur for a pure double-layer formulation on multiply connected geometries.

## Attributes
* `sampler`: Placeholder periodic sampler retained for the common solver API.
* `pts_scaling_factor`: Boundary-resolution scaling factors.
* `dim_scaling_factor`: Compatibility field for the generic solver interface.
* `eps`: Numerical tolerance placeholder.
* `min_dim`: Compatibility field for the generic solver interface.
* `min_pts`: Minimum number of points on each boundary component.
* `billiard`: Underlying billiard geometry.
* `symmetry`: Optional symmetry descriptor.

## Notes
For a boundary component of length `L`, the nominal node count is

    N ≈ k*L*b/(2π),

where `b` is the corresponding boundary-resolution scaling factor.

Use this solver when every connected component is represented by one smooth
closed curve. Composite piecewise smooth components should instead use
[`CFIE_kress_global_corners`](@ref).
"""
struct CFIE_kress{T<:Real,Bi<:BilliardGeometry.AbsBilliard,Sym}<:CFIE
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
    CFIE_kress_corners{T,Bi,Sym} <: CFIE

Combined-field integral equation solver using a Kress grading transformation on
a single closed curve with corner-type parameter singularities.

## Description
This is the graded counterpart of [`CFIE_kress`](@ref) for a boundary component
that remains represented by one closed parameterized curve but requires
endpoint/corner clustering.

A uniform computational variable `σ` is transformed through a nonlinear map

    t=w(σ),

and the geometric derivatives are transformed by the chain rule. The resulting
operator remains

    A(k)=I-(D(k)+ikS(k)).

The Kress logarithmic splitting is performed in the computational periodic
variable.

## Attributes
* `sampler`: Placeholder periodic sampler retained for the common solver API.
* `pts_scaling_factor`: Boundary-resolution scaling factors.
* `dim_scaling_factor`: Compatibility field for the generic solver interface.
* `eps`: Numerical tolerance placeholder.
* `min_dim`: Compatibility field for the generic solver interface.
* `min_pts`: Minimum number of boundary points.
* `billiard`: Underlying billiard geometry.
* `symmetry`: Optional symmetry descriptor.
* `kressq`: Order of the Kress grading transformation.
* `min_t_spacing`: Minimum permitted physical-parameter spacing after grading.

## Notes
For `Float64`, excessively large grading orders may force neighboring mapped
nodes below machine-resolvable spacing. A grading order around `4` is therefore
the practical default.
"""
struct CFIE_kress_corners{T<:Real,Bi<:BilliardGeometry.AbsBilliard,Sym}<:CFIE
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

"""
    CFIE_kress_global_corners{T,Bi,Sym} <: CFIE

Combined-field integral equation solver using global periodic Kress grading on
closed composite boundary components.

## Description
This is the general Kress solver for closed components represented by several
joined curve segments.

For each connected component, all constituent segments are treated as one
global periodic boundary. True corner locations are detected from tangent jumps.
If corners are present, a global Kress grading map is constructed relative to
all corner locations simultaneously.

If the joins are smooth, the component is instead discretized on an ungraded
uniform periodic mesh.

The resulting Fredholm operator is

    A(k)=I-(D(k)+ikS(k)).

Treating each closed component globally is essential because the periodic
logarithmic singularity belongs to the complete closed boundary rather than to
the individual constituent segments.

## Attributes
* `sampler`: Placeholder periodic sampler retained for the common solver API.
* `pts_scaling_factor`: Boundary-resolution scaling factors.
* `dim_scaling_factor`: Compatibility field for the generic solver interface.
* `eps`: Numerical tolerance placeholder.
* `min_dim`: Compatibility field for the generic solver interface.
* `min_pts`: Minimum number of boundary points.
* `billiard`: Underlying billiard geometry.
* `symmetry`: Optional symmetry descriptor.
* `kressq`: Global Kress grading order.
* `min_t_spacing`: Minimum permitted physical-parameter spacing after grading.

## Notes
This solver supports composite outer boundaries and composite hole boundaries.
Each connected component is discretized as one periodic object.
"""
struct CFIE_kress_global_corners{T<:Real,Bi<:BilliardGeometry.AbsBilliard,Sym}<:CFIE
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

"""
    CFIE_kress(pts_scaling_factor,billiard;min_pts=20,eps=1e-15,symmetry=nothing) → solver::CFIE_kress

Constructs a smooth periodic Kress combined-field solver.

## Arguments
* `pts_scaling_factor`: Boundary-resolution scaling factor or vector of factors.
* `billiard`: Billiard geometry.

## Keyword arguments
* `min_pts`: Minimum number of points on each boundary component.
* `eps`: Numerical tolerance placeholder.
* `symmetry`: Optional reflection or rotation symmetry descriptor.

## Returns
* `solver`: Configured [`CFIE_kress`](@ref) instance.
"""
function CFIE_kress(pts_scaling_factor::Union{T,Vector{T}},billiard::Bi;min_pts=200,eps=T(1e-15),symmetry::Union{Nothing,AbsSymmetry}=nothing) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    bs=pts_scaling_factor isa T ? [pts_scaling_factor] : pts_scaling_factor
    sampler=[BilliardGeometry.LinearNodes()]
    Sym=typeof(symmetry)
    return CFIE_kress{T,Bi,Sym}(sampler,bs,bs[1],eps,min_pts,min_pts,billiard,symmetry)
end

"""
    CFIE_kress_corners(pts_scaling_factor,billiard;min_pts=200,eps=1e-15,symmetry=nothing,kressq=2,min_t_spacing=1e-12) → solver::CFIE_kress_corners

Constructs a single-curve corner-graded Kress combined-field solver.

## Arguments
* `pts_scaling_factor`: Boundary-resolution scaling factor or vector of factors.
* `billiard`: Billiard geometry.

## Keyword arguments
* `min_pts`: Minimum number of boundary points.
* `eps`: Numerical tolerance placeholder.
* `symmetry`: Optional reflection or rotation symmetry descriptor.
* `kressq`: Kress grading order.
* `min_t_spacing`: Minimum permitted mapped-parameter spacing.

## Returns
* `solver`: Configured [`CFIE_kress_corners`](@ref) instance.

## Notes
For `Float64`, values of `kressq` substantially larger than `4` should be used
with care because the graded nodes may become indistinguishable at machine
precision.
"""
function CFIE_kress_corners(pts_scaling_factor::Union{T,Vector{T}},billiard::Bi;min_pts=200,eps=T(1e-15),symmetry::Union{Nothing,AbsSymmetry}=nothing,kressq=2,min_t_spacing=1e-12) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    bs=pts_scaling_factor isa T ? [pts_scaling_factor] : pts_scaling_factor
    sampler=[BilliardGeometry.LinearNodes()]
    Sym=typeof(symmetry)
    return CFIE_kress_corners{T,Bi,Sym}(sampler,bs,bs[1],eps,min_pts,min_pts,billiard,symmetry,kressq,min_t_spacing)
end

"""
    CFIE_kress_global_corners(pts_scaling_factor,billiard;min_pts=200,eps=1e-15,symmetry=nothing,kressq=2,min_t_spacing=1e-12) → solver::CFIE_kress_global_corners

Constructs a globally graded Kress combined-field solver for composite boundary
components.

## Arguments
* `pts_scaling_factor`: Boundary-resolution scaling factor or vector of factors.
* `billiard`: Billiard geometry.

## Keyword arguments
* `min_pts`: Minimum number of boundary points.
* `eps`: Numerical tolerance placeholder.
* `symmetry`: Optional reflection or rotation symmetry descriptor.
* `kressq`: Global Kress grading order.
* `min_t_spacing`: Minimum permitted mapped-parameter spacing.

## Returns
* `solver`: Configured [`CFIE_kress_global_corners`](@ref) instance.
"""
function CFIE_kress_global_corners(pts_scaling_factor::Union{T,Vector{T}},billiard::Bi;min_pts=200,eps=T(1e-15),symmetry::Union{Nothing,AbsSymmetry}=nothing,kressq=2,min_t_spacing=1e-12) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    bs=pts_scaling_factor isa T ? [pts_scaling_factor] : pts_scaling_factor
    sampler=[BilliardGeometry.LinearNodes()]
    Sym=typeof(symmetry)
    return CFIE_kress_global_corners{T,Bi,Sym}(sampler,bs,bs[1],eps,min_pts,min_pts,billiard,symmetry,kressq,min_t_spacing)
end

"""
    CFIE_kress_composite_solver{T,Bi,Sym,CS} <: CFIE

Combined-field integral equation solver for multiply connected geometries whose
connected boundary components require different Kress discretizations.

## Description
`CFIE_kress_composite_solver` assigns one existing CFIE-Kress solver to each connected
physical boundary component. This allows, for example,

    outer boundary -> CFIE_kress
    inner polygon  -> CFIE_kress_global_corners

while assembling one globally coupled CFIE operator

    A(k)=I-(D(k)+ikS(k)).

The component solvers control only the discretization and same-component Kress
quadrature. Inter-component interactions remain smooth and are evaluated with
ordinary Nyström quadrature.

The first connected component is interpreted as the outer boundary. Components
`2:end` are interpreted as holes and are orientation reversed after
discretization.

## Attributes
* `component_solvers::CS`: Tuple containing one CFIE solver per connected boundary component.
* `sampler::Vector{BilliardGeometry.LinearNodes}`: Placeholder sampler retained for the common solver API.
* `pts_scaling_factor::Vector{T}`: Representative boundary-resolution factors.
* `dim_scaling_factor::T`: Compatibility field for the generic solver interface.
* `eps::T`: Numerical tolerance placeholder.
* `min_dim::Int64`: Compatibility field for the generic solver interface.
* `min_pts::Int64`: Minimum point count used for generic solver metadata.
* `billiard::Bi`: Underlying billiard geometry.
* `symmetry::Sym`: Optional active symmetry descriptor.

## Notes
The number of component solvers must equal the number of connected components in
`billiard.full_boundary`.

All component solvers must use the same active symmetry as the composite solver,
because symmetry-compatible node counts and the global symmetry-orbit map must
refer to one common representation.
"""
struct CFIE_kress_composite_solver{T<:Real,Bi<:BilliardGeometry.AbsBilliard,Sym,CS<:Tuple}<:CFIE
    component_solvers::CS
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
    CFIE_kress_composite_solver(component_solvers::Tuple,billiard::Bi;symmetry=nothing) where {Bi<:BilliardGeometry.AbsBilliard} → CFIE_kress_composite_solver

Construct a component-wise CFIE-Kress solver.

## Arguments
* `component_solvers::Tuple`: One `CFIE` solver for each connected physical boundary component.
* `billiard::Bi`: Multiply connected billiard geometry.

## Keyword Arguments
* `symmetry::Union{Nothing,AbsSymmetry}=nothing`: Active global symmetry descriptor.

## Returns
* `solver::CFIE_kress_composite_solver`: Component-wise CFIE solver.

## Throws
* `DimensionMismatch`: If the number of component solvers differs from the number of connected boundary components.
* `ArgumentError`: If a component solver uses a symmetry different from the global symmetry.
"""
function CFIE_kress_composite_solver(component_solvers::Tuple,billiard::Bi;symmetry::Union{Nothing,AbsSymmetry}=nothing) where {Bi<:BilliardGeometry.AbsBilliard}
    comps=_boundary_components(billiard.full_boundary)
    length(component_solvers)==length(comps)||throw(DimensionMismatch("Received $(length(component_solvers)) component solvers for $(length(comps)) connected boundary components"))
    all(s->s isa CFIE,component_solvers)||throw(ArgumentError("Every component solver must be a CFIE solver"))
    for (i,s) in enumerate(component_solvers) # there is still a global symmetry, so all component solvers must use the same symmetry
        isequal(s.symmetry,symmetry)||throw(ArgumentError("Component solver $i uses symmetry $(s.symmetry), but the composite solver uses symmetry $symmetry"))
    end
    T=promote_type(map(s->eltype(s.pts_scaling_factor),component_solvers)...)
    bs=T[s.pts_scaling_factor[1] for s in component_solvers] # each boundary has its own ppw
    eps=maximum(T(s.eps) for s in component_solvers)
    min_pts=maximum(s.min_pts for s in component_solvers)
    sampler=[BilliardGeometry.LinearNodes()]
    return CFIE_kress_composite_solver{T,Bi,typeof(symmetry),typeof(component_solvers)}(component_solvers,sampler,bs,bs[1],eps,min_pts,min_pts,billiard,symmetry)
end

#### Equispaced periodic parameters ####
@inline s(k::Int,N::Int)=two_pi*k/N
@inline s_mid(k::Int,N::Int)=two_pi*(k-0.5)/N

"""
    _reverse_component_orientation(solver::CFIE,pts::BoundaryPoints{T}) where {T<:Real} → pts_reversed::BoundaryPoints{T}

Reverses the orientation of one boundary component.

## Description
For multiply connected geometries, the outer boundary and hole boundaries must
carry opposite orientations.

The point ordering is reversed together with the first and second derivative
data. First derivatives change sign under orientation reversal, while second
derivatives retain their sign.

The endpoint metadata of open panels is also exchanged consistently.

The returned [`BoundaryPoints`](@ref) object automatically reconstructs its
outward-normal data from the reversed tangent orientation.

## Arguments
* `solver`: Combined-field solver.
* `pts`: Boundary component to reverse.

## Returns
* `pts_reversed`: New [`BoundaryPoints`](@ref) instance with reversed orientation.
"""
function _reverse_component_orientation(solver::S,pts::BoundaryPoints{T}) where {T<:Real,S<:CFIE}
    xy=reverse(pts.xy)
    tangent=reverse(-pts.tangent)
    tangent_2=reverse(pts.tangent_2)
    ts=reverse(pts.ts)
    tphys=reverse(pts.tphys)
    ws=reverse(pts.ws)
    ws_der=reverse(pts.ws_der)
    L=sum(pts.ds)
    s=L.-reverse(pts.s)
    ds=reverse(pts.ds)
    xL=pts.xR
    xR=pts.xL
    tL=-pts.tR
    tR=-pts.tL
    return BoundaryPoints(xy,tangent,tangent_2,ts,tphys,ws,ws_der,s,ds,pts.compid,pts.is_periodic,xL,xR,tL,tR)
end

###############
#### KRESS ####
###############

"""
    _evaluate_points(solver::CFIE_kress{T},crv::C,k::T,idx::Int) where {T<:Real,C<:BilliardGeometry.AbsCurve} → pts::BoundaryPoints{T}

Constructs the smooth periodic Kress discretization of one closed boundary
component.

## Description
For a component of length `L`, the nominal number of nodes is

    N ≈ k*L*b/(2π),

subject to the minimum node count and symmetry compatibility.

The Kress computational variable lies in `[0,2π)`, whereas the geometry curves
are parameterized on `[0,1]`. Therefore

    u=t/(2π),

and the geometric derivatives transform as

    γ_t=γ_u/(2π),

and

    γ_tt=γ_uu/(2π)².

The periodic quadrature weights are

    ws_j=2π/N.

Since no grading is applied,

    ws_der_j=1.

## Arguments
* `solver`: Smooth periodic Kress solver.
* `crv`: Smooth closed boundary curve.
* `k`: Real wavenumber controlling the node density.
* `idx`: Connected-component identifier.

## Returns
* `pts`: Periodic [`BoundaryPoints`](@ref) discretization.
"""
function _evaluate_points(solver::CFIE_kress{T},crv::C,k::T,idx::Int) where {T<:Real,C<:BilliardGeometry.AbsCurve}
    L=T(crv.length)
    N=max(solver.min_pts,round(Int,k*L*solver.pts_scaling_factor[1]/two_pi))
    needed=isnothing(solver.symmetry) ? 2 : lcm(2,symmetry_node_multiple(solver.symmetry))
    N=cld(N,needed)*needed
    ts=T[s_mid(j,N) for j in 1:N]
    tphys=ts./T(two_pi)
    xy=BilliardGeometry.curve(crv,tphys)
    tangent_1st=tangent(crv,tphys)./T(two_pi)
    tangent_2nd=tangent_2(crv,tphys)./T(two_pi)^2
    s=BilliardGeometry.arc_length(crv,tphys)
    h=T(two_pi)/T(N)
    ds=Vector{T}(undef,N)
    @inbounds for i in 1:N
        v=tangent_1st[i]
        ds[i]=hypot(v[1],v[2])*h
    end
    ws=fill(h,N)
    ws_der=ones(T,N)
    z=SVector{2,T}(zero(T),zero(T))
    return BoundaryPoints(xy,tangent_1st,tangent_2nd,ts,tphys,ws,ws_der,s,ds,idx,true,z,z,z,z)
end

"""
    _evaluate_points(solver::CFIE_kress_corners{T},crv::C,k::T,idx::Int) where {T<:Real,C<:BilliardGeometry.AbsCurve} → pts::BoundaryPoints{T}

Constructs a Kress-graded periodic discretization of one closed curve.

## Description
The uniform computational variable `σ` is mapped to a physical periodic
parameter

    t=t(σ),

with first and second derivatives `jac` and `jac2`.

Because the underlying geometry parameter is

    u=t/(2π),

the transformed derivatives are

    γ_σ=γ_u*jac/(2π),

and

    γ_σσ=γ_uu*(jac/(2π))²+γ_u*jac2/(2π).

The grading Jacobian is retained in `ws_der`.

## Arguments
* `solver`: Single-curve corner-graded Kress solver.
* `crv`: Closed boundary curve.
* `k`: Real wavenumber controlling the node density.
* `idx`: Connected-component identifier.

## Returns
* `pts`: Graded periodic [`BoundaryPoints`](@ref) discretization.
"""
function _evaluate_points(solver::CFIE_kress_corners{T},crv::C,k::T,idx::Int) where {T<:Real,C<:BilliardGeometry.AbsCurve}
    L=T(crv.length)
    N=max(solver.min_pts,round(Int,k*L*solver.pts_scaling_factor[1]/two_pi))
    needed=isnothing(solver.symmetry) ? 2 : lcm(2,symmetry_node_multiple(solver.symmetry))
    N=cld(N,needed)*needed
    σ,tmap,jac,jac2,_=kress_graded_nodes_data(T,N;q=solver.kressq,minsep_tol=solver.min_t_spacing)
    tphys=tmap./T(two_pi)
    xy=BilliardGeometry.curve(crv,tphys)
    γu=tangent(crv,tphys)
    γuu=tangent_2(crv,tphys)
    tangent_1st=Vector{SVector{2,T}}(undef,N)
    tangent_2nd=Vector{SVector{2,T}}(undef,N)
    @inbounds for i in 1:N
        a=jac[i]/T(two_pi)
        b=jac2[i]/T(two_pi)
        tangent_1st[i]=γu[i]*a
        tangent_2nd[i]=γuu[i]*a^2+γu[i]*b
    end
    s=BilliardGeometry.arc_length(crv,tphys)
    h=T(two_pi)/T(N)
    ds=Vector{T}(undef,N)
    @inbounds for i in 1:N
        v=tangent_1st[i]
        ds[i]=hypot(v[1],v[2])*h
    end
    ws=fill(h,N)
    z=SVector{2,T}(zero(T),zero(T))
    return BoundaryPoints(xy,tangent_1st,tangent_2nd,σ,tphys,ws,jac,s,ds,idx,true,z,z,z,z)
end

############################
#### KRESS MULTI CORNER ####
############################

"""
    _evaluate_points_smooth_composite(solver::CFIE_kress_global_corners{T},comp::Vector{C},k::T,idx::Int) where {T<:Real,C<:BilliardGeometry.AbsCurve} → pts::BoundaryPoints{T}

Constructs an ungraded periodic discretization of one smooth composite boundary
component.

## Description
This is the fallback used by [`CFIE_kress_global_corners`](@ref) when several
joined curve pieces form one closed component but all junctions are smooth.

A global parameter

    t∈[0,2π)

is used for the complete component and evaluated through
[`_eval_composite_geom_global_t`](@ref).

Since no grading is active, the quadrature is the ordinary periodic
trapezoidal rule.

## Arguments
* `solver`: Global-composite Kress solver.
* `comp`: Curve pieces forming one smooth closed component.
* `k`: Real wavenumber controlling the node density.
* `idx`: Connected-component identifier.

## Returns
* `pts`: Ungraded periodic [`BoundaryPoints`](@ref) discretization.
"""
function _evaluate_points_smooth_composite(solver::CFIE_kress_global_corners{T},comp::Vector{C},k::T,idx::Int) where {T<:Real,C<:BilliardGeometry.AbsCurve}
    _,_,Ltot=component_lengths(comp)
    N=max(solver.min_pts,round(Int,k*Ltot*solver.pts_scaling_factor[1]/two_pi))
    needed=isnothing(solver.symmetry) ? 2 : lcm(2,symmetry_node_multiple(solver.symmetry))
    N=cld(N,needed)*needed
    ts=T[s_mid(j,N) for j in 1:N]
    tphys=copy(ts)
    h=T(two_pi)/T(N)
    xy=Vector{SVector{2,T}}(undef,N)
    tangent_1st=Vector{SVector{2,T}}(undef,N)
    tangent_2nd=Vector{SVector{2,T}}(undef,N)
    s=Vector{T}(undef,N)
    ds=Vector{T}(undef,N)
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
    _evaluate_points(solver::CFIE_kress_global_corners{T},comp::Vector{C},k::T,idx::Int) where {T<:Real,C<:BilliardGeometry.AbsCurve} → pts::BoundaryPoints{T}

Constructs the global periodic Kress discretization of a composite closed
boundary component.

## Description
True corner locations are detected by
[`_component_corner_locations`](@ref).

If no true corners are present, the function delegates to
[`_evaluate_points_smooth_composite`](@ref).

Otherwise the global computational coordinate `σ` is transformed through a
multi-corner grading map

    t=t(σ).

For global composite derivatives,

    γ_σ=γ_t t_σ,

and

    γ_σσ=γ_tt(t_σ)²+γ_t t_σσ.

The physical arc-length quadrature element is

    ds_j=|γ_σ|*2π/N.

## Arguments
* `solver`: Global-composite Kress solver.
* `comp`: Curve segments forming one closed component.
* `k`: Real wavenumber controlling the node density.
* `idx`: Connected-component identifier.

## Returns
* `pts`: Globally graded [`BoundaryPoints`](@ref) discretization.
"""
function _evaluate_points(solver::CFIE_kress_global_corners{T},comp::Vector{C},k::T,idx::Int) where {T<:Real,C<:BilliardGeometry.AbsCurve}
    corners=_component_corner_locations(T,comp)
    isempty(corners)&&return _evaluate_points_smooth_composite(solver,comp,k,idx)
    _,_,Ltot=component_lengths(comp)
    N=max(solver.min_pts,round(Int,k*Ltot*solver.pts_scaling_factor[1]/two_pi))
    needed=isnothing(solver.symmetry) ? 2 : lcm(2,symmetry_node_multiple(solver.symmetry))
    N=cld(N,needed)*needed
    σ,tmap,jac,jac2,_=multi_kress_graded_nodes_data(T,N,corners;q=solver.kressq,minsep_tol=solver.min_t_spacing)
    tphys=tmap
    h=T(two_pi)/T(N)
    xy=Vector{SVector{2,T}}(undef,N)
    tangent_1st=Vector{SVector{2,T}}(undef,N)
    tangent_2nd=Vector{SVector{2,T}}(undef,N)
    s=Vector{T}(undef,N)
    ds=Vector{T}(undef,N)
    @inbounds for i in 1:N
        q,γt,γtt=_eval_composite_geom_global_t(T,comp,tphys[i])
        xy[i]=q
        tangent_1st[i]=γt*jac[i]
        tangent_2nd[i]=γtt*jac[i]^2+γt*jac2[i]
        s[i]=_composite_arclength(comp,tphys[i])
        v=tangent_1st[i]
        ds[i]=hypot(v[1],v[2])*h
    end
    ws=fill(h,N)
    z=SVector{2,T}(zero(T),zero(T))
    return BoundaryPoints(xy,tangent_1st,tangent_2nd,σ,tphys,ws,jac,s,ds,idx,true,z,z,z,z)
end

####################
#### HIGH LEVEL ####
####################

"""
    evaluate_points(solver::Union{CFIE_kress{T},CFIE_kress_corners{T}},billiard::Bi,k::T) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard} → pts::Vector{BoundaryPoints{T}}

Constructs the Kress boundary discretizations of all connected boundary
components.

## Description
Every connected component must consist of exactly one closed curve.

The first component is interpreted as the outer boundary. Every subsequent
component is interpreted as a hole and has its orientation reversed through
[`_reverse_component_orientation`](@ref).

## Arguments
* `solver`: Smooth or single-curve graded Kress solver.
* `billiard`: Billiard geometry.
* `k`: Real wavenumber controlling the boundary resolution.

## Returns
* `pts`: Vector of [`BoundaryPoints`](@ref), one per connected component.

## Notes
Composite multi-segment components require
[`CFIE_kress_global_corners`](@ref).
"""
function evaluate_points(solver::Union{CFIE_kress{T},CFIE_kress_corners{T}},billiard::Bi,k::T) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    comps=_boundary_components(billiard.full_boundary)
    isempty(comps)&&error("Boundary cannot be empty.")
    pts=Vector{BoundaryPoints{T}}(undef,length(comps))
    for (idx,comp) in enumerate(comps)
        isempty(comp)&&error("Boundary component cannot be empty.")
        length(comp)==1||error("Periodic Kress requires each boundary component to be represented by one closed curve. Use CFIE_kress_global_corners or CFIE_kress_composite_solver for composite components.")
        p=_evaluate_points(solver,comp[1],k,idx)
        pts[idx]=idx==1 ? p : _reverse_component_orientation(solver,p)
    end
    return pts
end

"""
    evaluate_points(solver::CFIE_kress_global_corners{T},billiard::Bi,k::T) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard} → pts::Vector{BoundaryPoints{T}}

Constructs globally periodic Kress discretizations for all connected boundary
components.

## Description
The function supports three geometry layouts:

1. one smooth closed curve,
2. one composite closed outer boundary,
3. multiple connected components, each smooth or composite.

For multiply connected geometries, all components after the first are
orientation-reversed to represent holes.

Smooth one-curve components are delegated to [`CFIE_kress`](@ref). Composite
components are treated globally and graded only when true corners are detected.

## Arguments
* `solver`: Global-composite Kress solver.
* `billiard`: Billiard geometry.
* `k`: Real wavenumber controlling the boundary resolution.

## Returns
* `pts`: Vector of [`BoundaryPoints`](@ref), one per connected component.
"""
function evaluate_points(solver::CFIE_kress_global_corners{T},billiard::Bi,k::T) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    comps=_boundary_components(billiard.full_boundary)
    isempty(comps)&&error("Boundary cannot be empty.")
    pts=Vector{BoundaryPoints{T}}(undef,length(comps))
    base=CFIE_kress(solver.pts_scaling_factor,solver.billiard;min_pts=solver.min_pts,eps=solver.eps,symmetry=solver.symmetry)
    for (idx,comp) in enumerate(comps)
        isempty(comp)&&error("Boundary component cannot be empty.")
        p=length(comp)==1 ? _evaluate_points(base,comp[1],k,idx) : _evaluate_points(solver,comp,k,idx)
        pts[idx]=idx==1 ? p : _reverse_component_orientation(solver,p)
    end
    return pts
end

"""
    evaluate_points(solver::CFIE_kress_composite_solver{T},billiard::Bi,k::T) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard} → pts::Vector{BoundaryPoints{T}}

Construct the boundary discretization of a multiply connected geometry using one
CFIE-Kress solver for each connected physical boundary component.

Each component is discretized independently by the corresponding entry in
`solver.component_solvers`. Smooth and single-curve graded solvers require one
closed curve, while `CFIE_kress_global_corners` handles composite components.
The first component retains its orientation and components `2:end` are reversed
to represent holes.

## Arguments
* `solver::CFIE_kress_composite_solver{T}`: Component-wise CFIE-Kress solver.
* `billiard::Bi`: Billiard geometry.
* `k::T`: Real wavenumber controlling boundary resolution.

## Returns
* `pts::Vector{BoundaryPoints{T}}`: One boundary discretization per connected component.

## Throws
* `DimensionMismatch`: If the number of component solvers does not match the number of connected physical boundary components.
* `ArgumentError`: If a component is empty or incompatible with its assigned solver.
"""
function evaluate_points(solver::CFIE_kress_composite_solver{T},billiard::Bi,k::T) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    comps=_boundary_components(billiard.full_boundary)
    length(comps)==length(solver.component_solvers)||throw(DimensionMismatch("Billiard has $(length(comps)) connected components but composite solver has $(length(solver.component_solvers)) component solvers"))
    pts=Vector{BoundaryPoints{T}}(undef,length(comps))
    @inbounds for i in eachindex(comps)
        comp=comps[i]
        isempty(comp)&&throw(ArgumentError("Boundary component $i cannot be empty"))
        s=solver.component_solvers[i]
        if s isa CFIE_kress
            length(comp)==1||throw(ArgumentError("CFIE_kress requires boundary component $i to contain exactly one closed smooth curve"))
            p=_evaluate_points(s,comp[1],k,i)
        elseif s isa CFIE_kress_corners
            length(comp)==1||throw(ArgumentError("CFIE_kress_corners requires boundary component $i to contain exactly one closed curve"))
            p=_evaluate_points(s,comp[1],k,i)
        elseif s isa CFIE_kress_global_corners
            if length(comp)==1
                base=CFIE_kress(s.pts_scaling_factor,s.billiard;min_pts=s.min_pts,eps=s.eps,symmetry=s.symmetry)
                p=_evaluate_points(base,comp[1],k,i)
            else
                p=_evaluate_points(s,comp,k,i)
            end
        else
            throw(ArgumentError("Unsupported component solver $(typeof(s)) for boundary component $i"))
        end
        pts[i]=i==1 ? p : _reverse_component_orientation(solver,p)
    end
    return pts
end

"""
    CFIEKressWorkspace{T,M,S}

Reusable geometry and symmetry cache for CFIE-Kress matrix assembly.

## Description
The workspace contains all quantities depending on the fixed boundary
discretization but not on the current wavenumber `k`.

For several connected components, each component has its own pairwise geometry
cache and flattened panel arrays, while `Rmat` contains the global
block-diagonal Kress logarithmic correction.

If symmetry is active, `orbits` stores the exact [`SymmetryOrbitMap`](@ref)
used to assemble the reduced operator directly without forming the full complex
matrix.

## Attributes
* `offs::Vector{Int}`: Component offsets into the full boundary matrix.
* `Rmat::M`: Global block-diagonal Kress logarithmic correction matrix.
* `Gs::Vector{BoundaryGeomCache{T}}`: Geometry cache for each connected component.
* `parr::Vector{BoundaryPanelArrays{T}}`: Flattened geometry arrays for each component.
* `Ntot::Int`: Full boundary matrix dimension.
* `symmetry::S`: Active symmetry descriptor or `nothing`.
* `orbits::Union{Nothing,SymmetryOrbitMap{T}}`: Exact symmetry-orbit map when reduction is active.
* `global_to_block::Vector{Int}`: Full global index to component index.
* `global_to_local::Vector{Int}`: Full global index to component-local index.
"""
struct CFIEKressWorkspace{T<:Real,M<:AbstractMatrix{T},S}
    offs::Vector{Int}
    Rmat::M
    Gs::Vector{BoundaryGeomCache{T}}
    parr::Vector{BoundaryPanelArrays{T}}
    Ntot::Int
    symmetry::S
    orbits::Union{Nothing,SymmetryOrbitMap{T}}
    global_to_block::Vector{Int}
    global_to_local::Vector{Int}
end

@inline _cfie_workspace_dim(ws::CFIEKressWorkspace)=isnothing(ws.orbits) ? ws.Ntot : fundamental_size(ws.orbits)
@inline function boundary_matrix_size(solver::CFIE,pts::Vector{BoundaryPoints{T}}) where {T<:Real}
    isnothing(solver.symmetry)&&return boundary_matrix_size(pts)
    return fundamental_size(symmetry_index_orbits(T,pts,solver.symmetry))
end

"""
    cfie_reduced_orbit_size(
        solver::CFIE,
        pts::Vector{BoundaryPoints{T}},
    ) where {T<:Real} → Int

Return the CFIE matrix dimension after applying the active symmetry reduction.

## Arguments
* `solver::CFIE`: CFIE-Kress solver.
* `pts::Vector{BoundaryPoints{T}}`: Full boundary discretizations.

## Returns
* `n::Int`: Full matrix dimension without symmetry or number of fundamental symmetry orbits with symmetry.
"""
function cfie_reduced_orbit_size(solver::CFIE,pts::Vector{BoundaryPoints{T}}) where {T<:Real}
    isnothing(solver.symmetry)&&return boundary_matrix_size(pts)
    return fundamental_size(symmetry_index_orbits(T,pts,solver.symmetry))
end

function global_to_component_local(pts::Vector{BoundaryPoints{T}}) where {T<:Real}
    offs=component_offsets(pts)
    Ntot=offs[end]-1
    global_to_block=Vector{Int}(undef,Ntot)
    global_to_local=Vector{Int}(undef,Ntot)
    @inbounds for a in eachindex(pts)
        Na=length(pts[a])
        off=offs[a]
        for j in 1:Na
            g=off+j-1
            global_to_block[g]=a
            global_to_local[g]=j
        end
    end
    return global_to_block,global_to_local
end

"""
    build_cfie_kress_workspace(
        solver::CFIE,
        pts::Vector{BoundaryPoints{T}},
    ) where {T<:Real} → CFIEKressWorkspace

Build a reusable CFIE-Kress workspace for a fixed boundary discretization.

## Description
The workspace caches component offsets, the global Kress correction matrix,
component geometry caches, flattened panel arrays and global-to-local index
maps.

If symmetry is active, the exact [`SymmetryOrbitMap`](@ref) is also constructed
and stored.

## Arguments
* `solver::CFIE`: CFIE-Kress solver.
* `pts::Vector{BoundaryPoints{T}}`: Full boundary discretizations.

## Returns
* `ws::CFIEKressWorkspace`: Reusable full or symmetry-reduced CFIE-Kress workspace.
"""
function build_cfie_kress_workspace(solver::CFIE,pts::Vector{BoundaryPoints{T}}) where {T<:Real}
    offs=component_offsets(pts)
    Rmat=build_Rmat_kress(solver,pts)
    Gs=[boundary_geom_cache(p,_is_nontrivial_grading(p)) for p in pts]
    parr=[_boundary_panel_arrays_cache(p) for p in pts]
    Ntot=offs[end]-1
    g2c,g2l=global_to_component_local(pts)
    orbits=isnothing(solver.symmetry) ? nothing : symmetry_index_orbits(T,pts,solver.symmetry)
    return CFIEKressWorkspace(offs,Rmat,Gs,parr,Ntot,solver.symmetry,orbits,g2c,g2l)
end

function _cfie_kress_workspace_from_Rmat(solver::CFIE,pts::Vector{BoundaryPoints{T}},Rmat::AbstractMatrix{T}) where {T<:Real}
    offs=component_offsets(pts)
    Gs=[boundary_geom_cache(p,_is_nontrivial_grading(p)) for p in pts]
    parr=[_boundary_panel_arrays_cache(p) for p in pts]
    Ntot=offs[end]-1
    g2c,g2l=global_to_component_local(pts)
    orbits=isnothing(solver.symmetry) ? nothing : symmetry_index_orbits(T,pts,solver.symmetry)
    return CFIEKressWorkspace(offs,Rmat,Gs,parr,Ntot,solver.symmetry,orbits,g2c,g2l)
end

"""
    build_Rmat_kress(solver::CFIE,pts::Vector{BoundaryPoints{T}}) where {T<:Real} → Rmat::Matrix{T}

Build the global block-diagonal Kress logarithmic correction matrix.

Each connected boundary component is treated independently.
Thus

    R=diag(R₁,R₂,...,R_nc),

with the correction rule selected from the actual discretization stored in each
[`BoundaryPoints`](@ref) object.

## Arguments
* `solver::CFIE`: CFIE solver. Retained for the common workspace API.
* `pts::Vector{BoundaryPoints{T}}`: Boundary discretizations.

## Returns
* `Rmat::Matrix{T}`: Global block-diagonal Kress correction matrix.
"""
function build_Rmat_kress(solver::CFIE,pts::Vector{BoundaryPoints{T}}) where {T<:Real}
    offs=component_offsets(pts)
    Ntot=offs[end]-1
    Rmat=zeros(T,Ntot,Ntot)
    @inbounds for a in eachindex(pts)
        ra=offs[a]:(offs[a+1]-1)
        kress_R!(@view Rmat[ra,ra])
    end
    return Rmat
end

################################################################################
############################ FULL MATRIX ASSEMBLY ##############################
################################################################################

"""
    construct_matrices!(
        solver::CFIE,
        A::AbstractMatrix{Complex{T}},
        pts::Vector{BoundaryPoints{T}},
        Rmat::AbstractMatrix{T},
        Gs::Vector{BoundaryGeomCache{T}},
        parr::Vector{BoundaryPanelArrays{T}},
        offs::Vector{Int},
        k::T;
        multithreaded::Bool=true,
    ) where {T<:Real} → A

Assemble the full CFIE-Kress Fredholm matrix

    A(k)=I-(D(k)+ikS(k)).

Same-component interactions use the Kress logarithmic split. Inter-component
interactions are smooth and use ordinary quadrature.

## Arguments
* `solver::CFIE`: CFIE-Kress solver.
* `A::AbstractMatrix{Complex{T}}`: Preallocated full destination matrix.
* `pts::Vector{BoundaryPoints{T}}`: Boundary discretizations.
* `Rmat::AbstractMatrix{T}`: Global Kress correction matrix.
* `Gs::Vector{BoundaryGeomCache{T}}`: Component geometry caches.
* `parr::Vector{BoundaryPanelArrays{T}}`: Component flattened geometry arrays.
* `offs::Vector{Int}`: Component offsets.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread sufficiently large assembly loops.

## Returns
* `A::AbstractMatrix{Complex{T}}`: Full CFIE-Kress Fredholm matrix.
"""
function construct_matrices!(solver::CFIE,A::AbstractMatrix{Complex{T}},pts::Vector{BoundaryPoints{T}},Rmat::AbstractMatrix{T},Gs::Vector{BoundaryGeomCache{T}},parr::Vector{BoundaryPanelArrays{T}},offs::Vector{Int},k::T;multithreaded::Bool=true) where {T<:Real}
    αL1=-k*inv_two_pi
    αL2=Complex{T}(0,k/2)
    αM1=-inv_two_pi
    αM2=Complex{T}(0,one(T)/2)
    ik=Complex{T}(0,k)
    fill!(A,zero(Complex{T}))
    nc=length(pts)
    for a in 1:nc
        pa=pts[a]
        Ga=Gs[a]
        Pa=parr[a]
        Na=length(Pa.X)
        ra=offs[a]:(offs[a+1]-1)
        @inbounds for i in 1:Na
            gi=ra[i]
            si=Ga.speed[i]
            κi=Ga.kappa[i]
            wi=pa.ws[i]
            dval=Complex{T}(wi*κi,zero(T))
            m1=αM1*si
            m2=((Complex{T}(0,one(T)/2)-euler_over_pi)-inv_two_pi*log((k^2/4)*si^2))*si
            sval=Complex{T}(Rmat[gi,gi]*m1,zero(T))+wi*m2
            A[gi,gi]=one(Complex{T})-(dval+ik*sval)
        end
        @use_threads multithreading=(multithreaded&&Na>=32) for j in 2:Na
            gj=ra[j]
            sj=Ga.speed[j]
            wj=pa.ws[j]
            @inbounds for i in 1:j-1
                gi=ra[i]
                si=Ga.speed[i]
                wi=pa.ws[i]
                r=Ga.R[i,j]
                invr=Ga.invR[i,j]
                lt=Ga.logterm[i,j]
                inn_ij=Ga.inner[i,j]
                inn_ji=Ga.inner[j,i]
                h0,h1=hankel_pair01(k*r)
                j0=real(h0)
                j1=real(h1)
                l1_ij=αL1*inn_ij*j1*invr
                l2_ij=αL2*inn_ij*h1*invr-l1_ij*lt
                dval_ij=Rmat[gi,gj]*l1_ij+wj*l2_ij
                m1_ij=αM1*j0*sj
                m2_ij=αM2*h0*sj-m1_ij*lt
                sval_ij=Rmat[gi,gj]*m1_ij+wj*m2_ij
                A[gi,gj]=-(dval_ij+ik*sval_ij)
                l1_ji=αL1*inn_ji*j1*invr
                l2_ji=αL2*inn_ji*h1*invr-l1_ji*lt
                dval_ji=Rmat[gj,gi]*l1_ji+wi*l2_ji
                m1_ji=αM1*j0*si
                m2_ji=αM2*h0*si-m1_ji*lt
                sval_ji=Rmat[gj,gi]*m1_ji+wi*m2_ji
                A[gj,gi]=-(dval_ji+ik*sval_ji)
            end
        end
    end
    for a in 1:nc,b in 1:nc
        a==b&&continue
        pb=pts[b]
        Pa=parr[a]
        Pb=parr[b]
        Na=length(Pa.X)
        Nb=length(Pb.X)
        ra=offs[a]:(offs[a+1]-1)
        rb=offs[b]:(offs[b+1]-1)
        Xa=Pa.X
        Ya=Pa.Y
        Xb=Pb.X
        Yb=Pb.Y
        dXb=Pb.dX
        dYb=Pb.dY
        sb=Pb.s
        @use_threads multithreading=(multithreaded&&Na>=16) for i in 1:Na
            gi=ra[i]
            xi=Xa[i]
            yi=Ya[i]
            @inbounds for j in 1:Nb
                gj=rb[j]
                dx=xi-Xb[j]
                dy=yi-Yb[j]
                r2=muladd(dx,dx,dy*dy)
                r2<=(eps(T))^2&&continue
                r=sqrt(r2)
                invr=inv(r)
                inn=dYb[j]*dx-dXb[j]*dy
                h0,h1=hankel_pair01(k*r)
                sj=sb[j]
                wj=pb.ws[j]
                dval=wj*(αL2*inn*h1*invr)
                sval=wj*(αM2*h0*sj)
                A[gi,gj]=-(dval+ik*sval)
            end
        end
    end
    return A
end

################################################################################
########################### REDUCED MATRIX ASSEMBLY ############################
################################################################################

"""
    construct_matrices_reduced!(
        solver::CFIE,
        A::AbstractMatrix{Complex{T}},
        pts::Vector{BoundaryPoints{T}},
        ws::CFIEKressWorkspace{T},
        k::T;
        multithreaded::Bool=true,
    ) where {T<:Real} → A

Assemble the CFIE-Kress matrix directly in the symmetry-reduced basis.

## Description
For a fundamental source orbit `b`, the reduced matrix element is

    A_red[a,b]=Σ_g χ(g)A_full[i_a,g·j_b].

All source quadrature weights are already included in the corresponding
full-space CFIE entries, so no additional orbit-weight ratio is required.

## Arguments
* `solver::CFIE`: CFIE-Kress solver.
* `A::AbstractMatrix{Complex{T}}`: Preallocated reduced destination matrix.
* `pts::Vector{BoundaryPoints{T}}`: Full boundary discretizations.
* `ws::CFIEKressWorkspace{T}`: Cached geometry and symmetry workspace.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread reduced column assembly.

## Returns
* `A::AbstractMatrix{Complex{T}}`: Symmetry-reduced CFIE-Kress matrix.
"""
function construct_matrices_reduced!(solver::CFIE,A::AbstractMatrix{Complex{T}},pts::Vector{BoundaryPoints{T}},ws::CFIEKressWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    orbits=ws.orbits
    isnothing(orbits)&&throw(ArgumentError("Reduced CFIE assembly requires an active symmetry"))
    Ifund=orbits.Ifund
    Nred=fundamental_size(orbits)
    ng=orbit_size(orbits)
    @assert size(A)==(Nred,Nred)
    αL1=-k*inv_two_pi
    αL2=Complex{T}(0,k/2)
    αM1=-inv_two_pi
    αM2=Complex{T}(0,one(T)/2)
    ik=Complex{T}(0,k)
    fill!(A,zero(Complex{T}))
    @use_threads multithreading=multithreaded for b in 1:Nred
        @inbounds for a in 1:Nred
            gi=Ifund[a]
            ib=ws.global_to_block[gi]
            i=ws.global_to_local[gi]
            acc=zero(Complex{T})
            for l in 1:ng
                gj=orbits.fund_to_full[l,b]
                scale=orbits.fund_to_scale[l,b]
                jb=ws.global_to_block[gj]
                j=ws.global_to_local[gj]
                if ib==jb
                    p=pts[ib]
                    G=ws.Gs[ib]
                    si=G.speed[i]
                    sj=G.speed[j]
                    wi=p.ws[i]
                    wj=p.ws[j]
                    if i==j
                        dval=Complex{T}(wi*G.kappa[i],zero(T))
                        m1=αM1*si
                        m2=((Complex{T}(0,one(T)/2)-euler_over_pi)-inv_two_pi*log((k^2/4)*si^2))*si
                        sval=Complex{T}(ws.Rmat[gi,gj]*m1,zero(T))+wi*m2
                        acc+=scale*(one(Complex{T})-(dval+ik*sval))
                    else
                        r=G.R[i,j]
                        invr=G.invR[i,j]
                        lt=G.logterm[i,j]
                        inn=G.inner[i,j]
                        h0,h1=hankel_pair01(k*r)
                        j0=real(h0)
                        j1=real(h1)
                        l1=αL1*inn*j1*invr
                        l2=αL2*inn*h1*invr-l1*lt
                        dval=ws.Rmat[gi,gj]*l1+wj*l2
                        m1=αM1*j0*sj
                        m2=αM2*h0*sj-m1*lt
                        sval=ws.Rmat[gi,gj]*m1+wj*m2
                        acc+=scale*-(dval+ik*sval)
                    end
                else
                    Pi=ws.parr[ib]
                    Pj=ws.parr[jb]
                    pb=pts[jb]
                    dx=Pi.X[i]-Pj.X[j]
                    dy=Pi.Y[i]-Pj.Y[j]
                    r2=muladd(dx,dx,dy*dy)
                    if r2>(eps(T))^2
                        r=sqrt(r2)
                        invr=inv(r)
                        inn=Pj.dY[j]*dx-Pj.dX[j]*dy
                        h0,h1=hankel_pair01(k*r)
                        sj=Pj.s[j]
                        wj=pb.ws[j]
                        dval=wj*(αL2*inn*h1*invr)
                        sval=wj*(αM2*h0*sj)
                        acc+=scale*-(dval+ik*sval)
                    end
                end
            end
            A[a,b]=acc
        end
    end
    return A
end

################################################################################
######################## FULL MATRIX WITH DERIVATIVES ##########################
################################################################################

"""
    construct_matrices!(
        solver::CFIE,
        A::AbstractMatrix{Complex{T}},
        A1::AbstractMatrix{Complex{T}},
        A2::AbstractMatrix{Complex{T}},
        pts::Vector{BoundaryPoints{T}},
        Rmat::AbstractMatrix{T},
        Gs::Vector{BoundaryGeomCache{T}},
        parr::Vector{BoundaryPanelArrays{T}},
        offs::Vector{Int},
        k::T;
        multithreaded::Bool=true,
    ) where {T<:Real} → A,A1,A2

Assemble the CFIE-Kress operator and its first two wavenumber derivatives.

## Description
For

    A(k)=I-(D(k)+ikS(k)),

the analytical derivatives are

    A'(k)=-(D'(k)+iS(k)+ikS'(k)),

and

    A''(k)=-(D''(k)+2iS'(k)+ikS''(k)).

## Arguments
* `solver::CFIE`: CFIE-Kress solver.
* `A::AbstractMatrix{Complex{T}}`: Destination matrix for `A(k)`.
* `A1::AbstractMatrix{Complex{T}}`: Destination matrix for `A'(k)`.
* `A2::AbstractMatrix{Complex{T}}`: Destination matrix for `A''(k)`.
* `pts::Vector{BoundaryPoints{T}}`: Boundary discretizations.
* `Rmat::AbstractMatrix{T}`: Global Kress correction matrix.
* `Gs::Vector{BoundaryGeomCache{T}}`: Component geometry caches.
* `parr::Vector{BoundaryPanelArrays{T}}`: Component flattened geometry arrays.
* `offs::Vector{Int}`: Component offsets.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread sufficiently large assembly loops.

## Returns
* `A::AbstractMatrix{Complex{T}}`: CFIE-Kress operator.
* `A1::AbstractMatrix{Complex{T}}`: First derivative.
* `A2::AbstractMatrix{Complex{T}}`: Second derivative.
"""
function construct_matrices!(solver::CFIE,A::AbstractMatrix{Complex{T}},A1::AbstractMatrix{Complex{T}},A2::AbstractMatrix{Complex{T}},pts::Vector{BoundaryPoints{T}},Rmat::AbstractMatrix{T},Gs::Vector{BoundaryGeomCache{T}},parr::Vector{BoundaryPanelArrays{T}},offs::Vector{Int},k::T;multithreaded::Bool=true) where {T<:Real}
    αL1=-k*inv_two_pi
    αL2=Complex{T}(0,k/2)
    αM1=-inv_two_pi
    αM2=Complex{T}(0,one(T)/2)
    ik=Complex{T}(0,k)
    fill!(A,zero(Complex{T}))
    fill!(A1,zero(Complex{T}))
    fill!(A2,zero(Complex{T}))
    nc=length(pts)
    for a in 1:nc
        pa=pts[a]
        Ga=Gs[a]
        Pa=parr[a]
        Na=length(Pa.X)
        ra=offs[a]:(offs[a+1]-1)
        @inbounds for i in 1:Na
            gi=ra[i]
            si=Ga.speed[i]
            κi=Ga.kappa[i]
            wi=pa.ws[i]
            dval=Complex{T}(wi*κi,zero(T))
            m1=αM1*si
            m2=((Complex{T}(0,one(T)/2)-euler_over_pi)-inv_two_pi*log((k^2/4)*si^2))*si
            sval=Complex{T}(Rmat[gi,gi]*m1,zero(T))+wi*m2
            sval1=wi*(-si/(pi*k))
            sval2=wi*(si/(pi*k^2))
            A[gi,gi]=one(Complex{T})-(dval+ik*sval)
            A1[gi,gi]=-(Complex{T}(0,1)*sval+ik*sval1)
            A2[gi,gi]=-(Complex{T}(0,2)*sval1+ik*sval2)
        end
        @use_threads multithreading=(multithreaded&&Na>=32) for j in 2:Na
            gj=ra[j]
            sj=Ga.speed[j]
            wj=pa.ws[j]
            @inbounds for i in 1:j-1
                gi=ra[i]
                si=Ga.speed[i]
                wi=pa.ws[i]
                r=Ga.R[i,j]
                invr=Ga.invR[i,j]
                lt=Ga.logterm[i,j]
                inn_ij=Ga.inner[i,j]
                inn_ji=Ga.inner[j,i]
                kr=k*r
                h0,h1=hankel_pair01(kr)
                j0=real(h0)
                j1=real(h1)
                l1_ij=αL1*inn_ij*j1*invr
                l2_ij=αL2*inn_ij*h1*invr-l1_ij*lt
                dval_ij=Rmat[gi,gj]*l1_ij+wj*l2_ij
                l1_ij_1=-(inn_ij*k*j0)*inv_two_pi
                l1_ij_2=(inn_ij*(k*r*j1-j0))*inv_two_pi
                l2_ij_1=(inn_ij*k*(lt*j0+im*pi*h0))*inv_two_pi
                l2_ij_2=(inn_ij*(lt*(j0-k*r*j1)+im*pi*(h0-k*r*h1)))*inv_two_pi
                dval_ij_1=Rmat[gi,gj]*l1_ij_1+wj*l2_ij_1
                dval_ij_2=Rmat[gi,gj]*l1_ij_2+wj*l2_ij_2
                m1_ij=αM1*j0*sj
                m2_ij=αM2*h0*sj-m1_ij*lt
                sval_ij=Rmat[gi,gj]*m1_ij+wj*m2_ij
                m1_ij_1=(r*sj*j1)*inv_two_pi
                m1_ij_2=(r*sj*(k*r*j0-j1))*inv_two_pi/k
                m2_ij_1=-(r*sj*(lt*j1+im*pi*h1))*inv_two_pi
                m2_ij_2=(r*sj*(lt*(j1-k*r*j0)-im*pi*k*r*h0+im*pi*h1))*inv_two_pi/k
                sval_ij_1=Rmat[gi,gj]*m1_ij_1+wj*m2_ij_1
                sval_ij_2=Rmat[gi,gj]*m1_ij_2+wj*m2_ij_2
                A[gi,gj]=-(dval_ij+ik*sval_ij)
                A1[gi,gj]=-(dval_ij_1+Complex{T}(0,1)*sval_ij+ik*sval_ij_1)
                A2[gi,gj]=-(dval_ij_2+Complex{T}(0,2)*sval_ij_1+ik*sval_ij_2)
                l1_ji=αL1*inn_ji*j1*invr
                l2_ji=αL2*inn_ji*h1*invr-l1_ji*lt
                dval_ji=Rmat[gj,gi]*l1_ji+wi*l2_ji
                l1_ji_1=-(inn_ji*k*j0)*inv_two_pi
                l1_ji_2=(inn_ji*(k*r*j1-j0))*inv_two_pi
                l2_ji_1=(inn_ji*k*(lt*j0+im*pi*h0))*inv_two_pi
                l2_ji_2=(inn_ji*(lt*(j0-k*r*j1)+im*pi*(h0-k*r*h1)))*inv_two_pi
                dval_ji_1=Rmat[gj,gi]*l1_ji_1+wi*l2_ji_1
                dval_ji_2=Rmat[gj,gi]*l1_ji_2+wi*l2_ji_2
                m1_ji=αM1*j0*si
                m2_ji=αM2*h0*si-m1_ji*lt
                sval_ji=Rmat[gj,gi]*m1_ji+wi*m2_ji
                m1_ji_1=(r*si*j1)*inv_two_pi
                m1_ji_2=(r*si*(k*r*j0-j1))*inv_two_pi/k
                m2_ji_1=-(r*si*(lt*j1+im*pi*h1))*inv_two_pi
                m2_ji_2=(r*si*(lt*(j1-k*r*j0)-im*pi*k*r*h0+im*pi*h1))*inv_two_pi/k
                sval_ji_1=Rmat[gj,gi]*m1_ji_1+wi*m2_ji_1
                sval_ji_2=Rmat[gj,gi]*m1_ji_2+wi*m2_ji_2
                A[gj,gi]=-(dval_ji+ik*sval_ji)
                A1[gj,gi]=-(dval_ji_1+Complex{T}(0,1)*sval_ji+ik*sval_ji_1)
                A2[gj,gi]=-(dval_ji_2+Complex{T}(0,2)*sval_ji_1+ik*sval_ji_2)
            end
        end
    end
    for a in 1:nc,b in 1:nc
        a==b&&continue
        pb=pts[b]
        Pa=parr[a]
        Pb=parr[b]
        Na=length(Pa.X)
        Nb=length(Pb.X)
        ra=offs[a]:(offs[a+1]-1)
        rb=offs[b]:(offs[b+1]-1)
        Xa=Pa.X
        Ya=Pa.Y
        Xb=Pb.X
        Yb=Pb.Y
        dXb=Pb.dX
        dYb=Pb.dY
        sb=Pb.s
        @use_threads multithreading=(multithreaded&&Na>=16) for i in 1:Na
            gi=ra[i]
            xi=Xa[i]
            yi=Ya[i]
            @inbounds for j in 1:Nb
                gj=rb[j]
                dx=xi-Xb[j]
                dy=yi-Yb[j]
                r2=muladd(dx,dx,dy*dy)
                r2<=(eps(T))^2&&continue
                r=sqrt(r2)
                invr=inv(r)
                inn=dYb[j]*dx-dXb[j]*dy
                h0,h1=hankel_pair01(k*r)
                sj=sb[j]
                wj=pb.ws[j]
                dval=wj*(Complex{T}(0,k/2)*inn*h1*invr)
                dval1=wj*(Complex{T}(0,one(T)/2)*inn*k*h0)
                dval2=wj*(Complex{T}(0,one(T)/2)*inn*(h0-k*r*h1))
                sval=wj*(Complex{T}(0,one(T)/2)*h0*sj)
                sval1=wj*(-Complex{T}(0,one(T)/2)*r*h1*sj)
                sval2=wj*(Complex{T}(0,one(T)/2)*r*(h1-k*r*h0)*sj/k)
                A[gi,gj]=-(dval+ik*sval)
                A1[gi,gj]=-(dval1+Complex{T}(0,1)*sval+ik*sval1)
                A2[gi,gj]=-(dval2+Complex{T}(0,2)*sval1+ik*sval2)
            end
        end
    end
    return A,A1,A2
end

################################################################################
####################### REDUCED MATRIX WITH DERIVATIVES ########################
################################################################################

"""
    construct_matrices_reduced_deriv!(
        solver::CFIE,
        A::AbstractMatrix{Complex{T}},
        A1::AbstractMatrix{Complex{T}},
        A2::AbstractMatrix{Complex{T}},
        pts::Vector{BoundaryPoints{T}},
        ws::CFIEKressWorkspace{T},
        k::T;
        multithreaded::Bool=true,
    ) where {T<:Real} → A,A1,A2

Assemble the symmetry-reduced CFIE-Kress operator and its first two analytical
wavenumber derivatives.

## Arguments
* `solver::CFIE`: CFIE-Kress solver.
* `A::AbstractMatrix{Complex{T}}`: Destination matrix for the reduced operator.
* `A1::AbstractMatrix{Complex{T}}`: Destination matrix for its first derivative.
* `A2::AbstractMatrix{Complex{T}}`: Destination matrix for its second derivative.
* `pts::Vector{BoundaryPoints{T}}`: Full boundary discretizations.
* `ws::CFIEKressWorkspace{T}`: Cached geometry and symmetry workspace.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread reduced column assembly.

## Returns
* `A::AbstractMatrix{Complex{T}}`: Reduced CFIE-Kress operator.
* `A1::AbstractMatrix{Complex{T}}`: First derivative.
* `A2::AbstractMatrix{Complex{T}}`: Second derivative.
"""
function construct_matrices_reduced_deriv!(solver::CFIE,A::AbstractMatrix{Complex{T}},A1::AbstractMatrix{Complex{T}},A2::AbstractMatrix{Complex{T}},pts::Vector{BoundaryPoints{T}},ws::CFIEKressWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    orbits=ws.orbits
    isnothing(orbits)&&throw(ArgumentError("Reduced CFIE assembly requires an active symmetry"))
    Ifund=orbits.Ifund
    Nred=fundamental_size(orbits)
    ng=orbit_size(orbits)
    @assert size(A)==(Nred,Nred)
    @assert size(A1)==(Nred,Nred)
    @assert size(A2)==(Nred,Nred)
    αL1=-k*inv_two_pi
    αL2=Complex{T}(0,k/2)
    αM1=-inv_two_pi
    αM2=Complex{T}(0,one(T)/2)
    ik=Complex{T}(0,k)
    fill!(A,zero(Complex{T}))
    fill!(A1,zero(Complex{T}))
    fill!(A2,zero(Complex{T}))
    @use_threads multithreading=multithreaded for b in 1:Nred
        @inbounds for a in 1:Nred
            gi=Ifund[a]
            ib=ws.global_to_block[gi]
            i=ws.global_to_local[gi]
            acc0=zero(Complex{T})
            acc1=zero(Complex{T})
            acc2=zero(Complex{T})
            for l in 1:ng
                gj=orbits.fund_to_full[l,b]
                scale=orbits.fund_to_scale[l,b]
                jb=ws.global_to_block[gj]
                j=ws.global_to_local[gj]
                if ib==jb
                    p=pts[ib]
                    G=ws.Gs[ib]
                    si=G.speed[i]
                    sj=G.speed[j]
                    wi=p.ws[i]
                    wj=p.ws[j]
                    if i==j
                        dval=Complex{T}(wi*G.kappa[i],zero(T))
                        m1=αM1*si
                        m2=((Complex{T}(0,one(T)/2)-euler_over_pi)-inv_two_pi*log((k^2/4)*si^2))*si
                        sval=Complex{T}(ws.Rmat[gi,gj]*m1,zero(T))+wi*m2
                        sval1=wi*(-si/(pi*k))
                        sval2=wi*(si/(pi*k^2))
                        acc0+=scale*(one(Complex{T})-(dval+ik*sval))
                        acc1+=scale*-(Complex{T}(0,1)*sval+ik*sval1)
                        acc2+=scale*-(Complex{T}(0,2)*sval1+ik*sval2)
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
                        dval=ws.Rmat[gi,gj]*l1+wj*l2
                        l1_1=-(inn*k*j0)*inv_two_pi
                        l1_2=(inn*(k*r*j1-j0))*inv_two_pi
                        l2_1=(inn*k*(lt*j0+im*pi*h0))*inv_two_pi
                        l2_2=(inn*(lt*(j0-k*r*j1)+im*pi*(h0-k*r*h1)))*inv_two_pi
                        dval1=ws.Rmat[gi,gj]*l1_1+wj*l2_1
                        dval2=ws.Rmat[gi,gj]*l1_2+wj*l2_2
                        m1=αM1*j0*sj
                        m2=αM2*h0*sj-m1*lt
                        sval=ws.Rmat[gi,gj]*m1+wj*m2
                        m1_1=(r*sj*j1)*inv_two_pi
                        m1_2=(r*sj*(k*r*j0-j1))*inv_two_pi/k
                        m2_1=-(r*sj*(lt*j1+im*pi*h1))*inv_two_pi
                        m2_2=(r*sj*(lt*(j1-k*r*j0)-im*pi*k*r*h0+im*pi*h1))*inv_two_pi/k
                        sval1=ws.Rmat[gi,gj]*m1_1+wj*m2_1
                        sval2=ws.Rmat[gi,gj]*m1_2+wj*m2_2
                        acc0+=scale*-(dval+ik*sval)
                        acc1+=scale*-(dval1+Complex{T}(0,1)*sval+ik*sval1)
                        acc2+=scale*-(dval2+Complex{T}(0,2)*sval1+ik*sval2)
                    end
                else
                    Pi=ws.parr[ib]
                    Pj=ws.parr[jb]
                    pb=pts[jb]
                    dx=Pi.X[i]-Pj.X[j]
                    dy=Pi.Y[i]-Pj.Y[j]
                    r2=muladd(dx,dx,dy*dy)
                    if r2>(eps(T))^2
                        r=sqrt(r2)
                        invr=inv(r)
                        inn=Pj.dY[j]*dx-Pj.dX[j]*dy
                        h0,h1=hankel_pair01(k*r)
                        sj=Pj.s[j]
                        wj=pb.ws[j]
                        dval=wj*(Complex{T}(0,k/2)*inn*h1*invr)
                        dval1=wj*(Complex{T}(0,one(T)/2)*inn*k*h0)
                        dval2=wj*(Complex{T}(0,one(T)/2)*inn*(h0-k*r*h1))
                        sval=wj*(Complex{T}(0,one(T)/2)*h0*sj)
                        sval1=wj*(-Complex{T}(0,one(T)/2)*r*h1*sj)
                        sval2=wj*(Complex{T}(0,one(T)/2)*r*(h1-k*r*h0)*sj/k)
                        acc0+=scale*-(dval+ik*sval)
                        acc1+=scale*-(dval1+Complex{T}(0,1)*sval+ik*sval1)
                        acc2+=scale*-(dval2+Complex{T}(0,2)*sval1+ik*sval2)
                    end
                end
            end
            A[a,b]=acc0
            A1[a,b]=acc1
            A2[a,b]=acc2
        end
    end
    return A,A1,A2
end

################################################################################
########################### HIGH-LEVEL INTERFACE ###############################
################################################################################

"""
    construct_matrices!(
        solver::CFIE,
        A::AbstractMatrix{Complex{T}},
        pts::Vector{BoundaryPoints{T}},
        ws::CFIEKressWorkspace{T},
        k::T;
        multithreaded::Bool=true,
    ) where {T<:Real} → A

    construct_matrices!(
        solver::CFIE,
        A::AbstractMatrix{Complex{T}},
        A1::AbstractMatrix{Complex{T}},
        A2::AbstractMatrix{Complex{T}},
        pts::Vector{BoundaryPoints{T}},
        ws::CFIEKressWorkspace{T},
        k::T;
        multithreaded::Bool=true,
    ) where {T<:Real} → A,A1,A2

Assemble the full or symmetry-reduced CFIE-Kress operator, selected
automatically from the cached workspace.

## Arguments
* `solver::CFIE`: CFIE-Kress solver.
* `A::AbstractMatrix{Complex{T}}`: Destination operator matrix.
* `A1::AbstractMatrix{Complex{T}}`: Destination first derivative.
* `A2::AbstractMatrix{Complex{T}}`: Destination second derivative.
* `pts::Vector{BoundaryPoints{T}}`: Full boundary discretizations.
* `ws::CFIEKressWorkspace{T}`: Cached CFIE-Kress workspace.
* `k::T`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread low-level assembly.

## Returns
* `A::AbstractMatrix{Complex{T}}` for the single-matrix overload.
* `(A,A1,A2)` for the derivative overload.
"""
function construct_matrices!(solver::CFIE,A::AbstractMatrix{Complex{T}},pts::Vector{BoundaryPoints{T}},ws::CFIEKressWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    if isnothing(ws.orbits)
        @blas_1 construct_matrices!(solver,A,pts,ws.Rmat,ws.Gs,ws.parr,ws.offs,k;multithreaded=multithreaded)
    else
        @blas_1 construct_matrices_reduced!(solver,A,pts,ws,k;multithreaded=multithreaded)
    end
    return A
end

function construct_matrices!(solver::CFIE,A::AbstractMatrix{Complex{T}},pts::Vector{BoundaryPoints{T}},Rmat::AbstractMatrix{T},k::T;multithreaded::Bool=true) where {T<:Real}
    ws=_cfie_kress_workspace_from_Rmat(solver,pts,Rmat)
    return construct_matrices!(solver,A,pts,ws,k;multithreaded=multithreaded)
end

function construct_matrices!(solver::CFIE,A::AbstractMatrix{Complex{T}},A1::AbstractMatrix{Complex{T}},A2::AbstractMatrix{Complex{T}},pts::Vector{BoundaryPoints{T}},ws::CFIEKressWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    if isnothing(ws.orbits)
        construct_matrices!(solver,A,A1,A2,pts,ws.Rmat,ws.Gs,ws.parr,ws.offs,k;multithreaded=multithreaded)
    else
        construct_matrices_reduced_deriv!(solver,A,A1,A2,pts,ws,k;multithreaded=multithreaded)
    end
    return A,A1,A2
end

function construct_matrices!(solver::CFIE,A::AbstractMatrix{Complex{T}},A1::AbstractMatrix{Complex{T}},A2::AbstractMatrix{Complex{T}},pts::Vector{BoundaryPoints{T}},Rmat::AbstractMatrix{T},k::T;multithreaded::Bool=true) where {T<:Real}
    ws=_cfie_kress_workspace_from_Rmat(solver,pts,Rmat)
    return construct_matrices!(solver,A,A1,A2,pts,ws,k;multithreaded=multithreaded)
end

function construct_matrices!(solver::CFIE,basis::AbstractHankelBasis,A::AbstractMatrix{Complex{T}},dA::AbstractMatrix{Complex{T}},ddA::AbstractMatrix{Complex{T}},pts::Vector{BoundaryPoints{T}},k::T;multithreaded::Bool=true) where {T<:Real}
    ws=build_cfie_kress_workspace(solver,pts)
    construct_matrices!(solver,A,dA,ddA,pts,ws,k;multithreaded=multithreaded)
    return A,dA,ddA
end

function construct_matrices!(solver::CFIE,basis::AbstractHankelBasis,A::AbstractMatrix{Complex{T}},dA::AbstractMatrix{Complex{T}},ddA::AbstractMatrix{Complex{T}},pts::Vector{BoundaryPoints{T}},ws::CFIEKressWorkspace{T},k::T;multithreaded::Bool=true) where {T<:Real}
    construct_matrices!(solver,A,dA,ddA,pts,ws,k;multithreaded=multithreaded)
    return A,dA,ddA
end

################################################################################
################################### SOLVE ######################################
################################################################################

"""
    solve(
        solver::CFIE,
        basis::Ba,
        pts::Vector{BoundaryPoints{T}},
        k;
        multithreaded::Bool=true,
        use_krylov::Bool=true,
        which::Symbol=:det,
    ) where {T<:Real,Ba<:AbsBasis}

Evaluate the selected scalar spectral diagnostic of the CFIE-Kress operator.

## Arguments
* `solver::CFIE`: CFIE-Kress solver.
* `basis::Ba`: Basis placeholder retained for the common solver API.
* `pts::Vector{BoundaryPoints{T}}`: Boundary discretizations.
* `k`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread matrix assembly.
* `use_krylov::Bool`: Whether the common spectral backend may use its Krylov pathway.
* `which::Symbol`: Spectral diagnostic to evaluate.

## Returns
* Scalar spectral diagnostic selected by `which`.
"""
function solve(solver::CFIE,basis::Ba,pts::Vector{BoundaryPoints{T}},k;multithreaded::Bool=true,use_krylov::Bool=true,which::Symbol=:det) where {T<:Real,Ba<:AbsBasis}
    ws=build_cfie_kress_workspace(solver,pts)
    N=_cfie_workspace_dim(ws)
    A=Matrix{Complex{T}}(undef,N,N)
    @blas_1 construct_matrices!(solver,A,pts,ws,k;multithreaded=multithreaded)
    @svd_or_det_solve A use_krylov which MAX_BLAS_THREADS
end

function solve(solver::CFIE,basis::Ba,pts::Vector{BoundaryPoints{T}},ws::CFIEKressWorkspace{T},k;multithreaded::Bool=true,use_krylov::Bool=true,which::Symbol=:det_argmin) where {T<:Real,Ba<:AbsBasis}
    N=_cfie_workspace_dim(ws)
    A=Matrix{Complex{T}}(undef,N,N)
    @blas_1 construct_matrices!(solver,A,pts,ws,k;multithreaded=multithreaded)
    @svd_or_det_solve A use_krylov which MAX_BLAS_THREADS
end

function solve(solver::CFIE,basis::Ba,A::AbstractMatrix{Complex{T}},pts::Vector{BoundaryPoints{T}},ws::CFIEKressWorkspace{T},k;multithreaded::Bool=true,use_krylov::Bool=true,which::Symbol=:det_argmin) where {T<:Real,Ba<:AbsBasis}
    @blas_1 construct_matrices!(solver,A,pts,ws,k;multithreaded=multithreaded)
    @svd_or_det_solve A use_krylov which MAX_BLAS_THREADS
end

################################################################################
################################ SOLVE VECT ####################################
################################################################################

"""
    solve_vect(
        solver::CFIE,
        basis::Ba,
        pts::Vector{BoundaryPoints{T}},
        ws::CFIEKressWorkspace{T},
        k;
        multithreaded::Bool=true,
        tol=1e-12,
        maxiter::Int=2000,
        krylovdim::Int=40,
    ) where {T<:Real,Ba<:AbsBasis} → σ,u

Compute a near-null vector of the CFIE-Kress matrix.

## Arguments
* `solver::CFIE`: CFIE-Kress solver.
* `basis::Ba`: Basis placeholder retained for API compatibility.
* `pts::Vector{BoundaryPoints{T}}`: Boundary discretizations.
* `ws::CFIEKressWorkspace{T}`: Cached CFIE-Kress workspace.
* `k`: Wavenumber.

## Keyword Arguments
* `multithreaded::Bool`: Whether to thread matrix assembly.
* `tol`: Krylov convergence tolerance.
* `maxiter::Int`: Maximum number of Krylov iterations.
* `krylovdim::Int`: Krylov subspace dimension.

## Returns
* `σ`: Near-zero eigenvalue-magnitude proxy returned by [`smallest_nullvec_krylov!`](@ref).
* `u`: Corresponding normalized near-null vector.
"""
function solve_vect(solver::CFIE,basis::Ba,A::AbstractMatrix{Complex{T}},pts::Vector{BoundaryPoints{T}},ws::CFIEKressWorkspace{T},k;multithreaded::Bool=true,tol=1e-12,maxiter::Int=2000,krylovdim::Int=40) where {T<:Real,Ba<:AbsBasis}
    @blas_1 construct_matrices!(solver,A,pts,ws,k;multithreaded=multithreaded)
    σ,u,_=smallest_nullvec_krylov!(A;nev=1,tol=tol,maxiter=maxiter,krylovdim=krylovdim)
    return σ,u
end

function solve_vect(solver::CFIE,basis::Ba,pts::Vector{BoundaryPoints{T}},ws::CFIEKressWorkspace{T},k;multithreaded::Bool=true,tol=1e-12,maxiter::Int=2000,krylovdim::Int=40) where {T<:Real,Ba<:AbsBasis}
    N=_cfie_workspace_dim(ws)
    A=Matrix{Complex{T}}(undef,N,N)
    return solve_vect(solver,basis,A,pts,ws,k;multithreaded=multithreaded,tol=tol,maxiter=maxiter,krylovdim=krylovdim)
end

function solve_vect(solver::CFIE,billiard::Bi,basis::Ba,pts::Vector{BoundaryPoints{T}},k;multithreaded::Bool=true,tol=1e-12,maxiter::Int=2000,krylovdim::Int=40) where {T<:Real,Ba<:AbsBasis,Bi<:BilliardGeometry.AbsBilliard}
    ws=build_cfie_kress_workspace(solver,pts)
    return solve_vect(solver,basis,pts,ws,k;multithreaded=multithreaded,tol=tol,maxiter=maxiter,krylovdim=krylovdim)
end

################################################################################
################################ SOLVE INFO ####################################
################################################################################

# INTERNAL - for benchmarking and diagnostics only.
function solve_INFO(solver::CFIE,basis::Ba,pts::Vector{BoundaryPoints{T}},ws::CFIEKressWorkspace{T},k;multithreaded::Bool=true,use_krylov::Bool=true,which::Symbol=:det_argmin) where {T<:Real,Ba<:AbsBasis}
    N=_cfie_workspace_dim(ws)
    A=Matrix{Complex{T}}(undef,N,N)
    t0=time()
    @info "Building boundary operator A from cached Kress workspace..."
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