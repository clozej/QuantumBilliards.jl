"""
    BoundaryPoints{T} <: AbsPoints

`BoundaryPoints` is the common boundary-discretization container used by the
quantum billiard solvers.

## Description
The type stores both the generic physical boundary data shared by all solvers and
additional parametrization data required by boundary-integral and Kress-type
discretizations.

The generic fields `xy`, `normal`, `s` and `ds` describe the sampled physical
boundary. Solver-specific fields such as `w`, `w_n`, `curvature` and `xy_int` are
populated only when required by the corresponding method. Parametric boundary
discretizations additionally populate `tangent`, `tangent_2`, `ts`, `tphys`,
`ws`, `ws_der` and the component metadata.

Fields that are not required by a given discretization default to empty vectors
or neutral scalar values. The inner constructor verifies that every non-empty
per-boundary-point vector has the same length as `xy`.

## Attributes
* `xy`: Boundary points in Cartesian coordinates.
* `normal`: Outward unit normal vectors at the boundary points.
* `s`: Arc-length coordinates of the boundary points.
* `ds`: Arc-length quadrature elements associated with the boundary points.
* `w`: Solver-specific boundary weights.
* `w_n`: Additional solver-specific normal weights.
* `curvature`: Curvature values at the boundary points.
* `xy_int`: Interior points used by methods requiring an interior sampling set.
* `shift_x`: Horizontal shift used when applying symmetry transformations.
* `shift_y`: Vertical shift used when applying symmetry transformations.
* `tangent`: First derivatives of the boundary parametrization at the nodes.
* `tangent_2`: Second derivatives of the boundary parametrization at the nodes.
* `ts`: Computational parameter values associated with the boundary nodes.
* `tphys`: Physical or global boundary parameter values associated with the nodes.
* `ws`: Quadrature weights in the computational parameter.
* `ws_der`: Derivatives of the computational quadrature weights.
* `compid`: Boundary-component index for multiply connected geometries.
* `is_periodic`: Whether the discretized boundary component is periodic.
* `xL`: Left endpoint of a non-periodic boundary component.
* `xR`: Right endpoint of a non-periodic boundary component.
* `tL`: Tangent vector at the left endpoint.
* `tR`: Tangent vector at the right endpoint.

## API
The following functions can be evaluated for this type:
- [`boundary_coords`](@ref)
- [`boundary_s`](@ref)
- [`boundary_matrix_size`](@ref)
- [`component_offsets`](@ref)
- `Base.length`
- `Base.isempty`
"""
struct BoundaryPoints{T<:Real}<:AbsPoints
    xy::Vector{SVector{2,T}}
    normal::Vector{SVector{2,T}}
    s::Vector{T}
    ds::Vector{T}
    w::Vector{T}
    w_n::Vector{T}
    curvature::Vector{T}
    xy_int::Vector{SVector{2,T}}
    shift_x::T
    shift_y::T
    tangent::Vector{SVector{2,T}}
    tangent_2::Vector{SVector{2,T}}
    ts::Vector{T}
    tphys::Vector{T}
    ws::Vector{T}
    ws_der::Vector{T}
    compid::Int
    is_periodic::Bool
    xL::SVector{2,T}
    xR::SVector{2,T}
    tL::SVector{2,T}
    tR::SVector{2,T}
    function BoundaryPoints{T}(xy,normal,s,ds,w,w_n,curvature,xy_int,shift_x,shift_y,tangent,tangent_2,ts,tphys,ws,ws_der,compid,is_periodic,xL,xR,tL,tR) where {T<:Real}
        n=length(xy)
        for (name,v) in ((:normal,normal),(:s,s),(:ds,ds),(:w,w),(:w_n,w_n),(:curvature,curvature),(:tangent,tangent),(:tangent_2,tangent_2),(:ts,ts),(:tphys,tphys),(:ws,ws),(:ws_der,ws_der))
            isempty(v)||length(v)==n||throw(DimensionMismatch("Length of $name ($(length(v))) must match xy ($n)"))
        end
        new{T}(xy,normal,s,ds,w,w_n,curvature,xy_int,shift_x,shift_y,tangent,tangent_2,ts,tphys,ws,ws_der,compid,is_periodic,xL,xR,tL,tR)
    end
end

"""
    BoundaryPoints(xy::Vector{SVector{2,T}}; kwargs...) where {T<:Real} → bp::BoundaryPoints{T}

Constructs a [`BoundaryPoints`](@ref) instance from the boundary points `xy`,
inferring the element type `T` and defaulting all optional data to empty vectors
or neutral values.

## Arguments
* `xy`: Boundary points in Cartesian coordinates.

## Keyword arguments
* `normal::Vector{SVector{2,T}} = SVector{2,T}[]`: Outward unit normal vectors.
* `s::Vector{T} = T[]`: Arc-length coordinates.
* `ds::Vector{T} = T[]`: Arc-length quadrature elements.
* `w::Vector{T} = T[]`: Solver-specific boundary weights.
* `w_n::Vector{T} = T[]`: Additional solver-specific normal weights.
* `curvature::Vector{T} = T[]`: Curvature values at the boundary points.
* `xy_int::Vector{SVector{2,T}} = SVector{2,T}[]`: Interior sampling points.
* `shift_x::T = zero(T)`: Horizontal symmetry shift.
* `shift_y::T = zero(T)`: Vertical symmetry shift.
* `tangent::Vector{SVector{2,T}} = SVector{2,T}[]`: First derivatives of the parametrization.
* `tangent_2::Vector{SVector{2,T}} = SVector{2,T}[]`: Second derivatives of the parametrization.
* `ts::Vector{T} = T[]`: Computational parameter nodes.
* `tphys::Vector{T} = T[]`: Physical or global boundary parameter values.
* `ws::Vector{T} = T[]`: Computational quadrature weights.
* `ws_der::Vector{T} = T[]`: Derivatives of the computational quadrature weights.
* `compid::Int = 1`: Boundary-component index.
* `is_periodic::Bool = true`: Whether the component is periodic.
* `xL::SVector{2,T}`: Left endpoint of a non-periodic component.
* `xR::SVector{2,T}`: Right endpoint of a non-periodic component.
* `tL::SVector{2,T}`: Tangent at the left endpoint.
* `tR::SVector{2,T}`: Tangent at the right endpoint.

## Returns
* `bp`: A [`BoundaryPoints{T}`](@ref) instance containing the supplied boundary data.
"""
function BoundaryPoints(xy::Vector{SVector{2,T}};normal=SVector{2,T}[],s=T[],ds=T[],w=T[],w_n=T[],curvature=T[],xy_int=SVector{2,T}[],shift_x=zero(T),shift_y=zero(T),tangent=SVector{2,T}[],tangent_2=SVector{2,T}[],ts=T[],tphys=T[],ws=T[],ws_der=T[],compid=1,is_periodic=true,xL=SVector{2,T}(zero(T),zero(T)),xR=SVector{2,T}(zero(T),zero(T)),tL=SVector{2,T}(zero(T),zero(T)),tR=SVector{2,T}(zero(T),zero(T))) where {T<:Real}
    return BoundaryPoints{T}(xy,normal,s,ds,w,w_n,curvature,xy_int,shift_x,shift_y,tangent,tangent_2,ts,tphys,ws,ws_der,compid,is_periodic,xL,xR,tL,tR)
end

"""
    BoundaryPoints(xy,normal,s,ds,w,w_n,curvature,xy_int,shift_x,shift_y) → bp::BoundaryPoints

Constructs a [`BoundaryPoints`](@ref) instance using the legacy positional
boundary-data layout.

## Description
This constructor is intended for solver discretizations that use the generic
boundary geometry together with method-specific weights but do not require the
additional parametric boundary data.

## Returns
* `bp`: A [`BoundaryPoints`](@ref) instance with the supplied generic and solver-specific fields populated.
"""
function BoundaryPoints(xy::Vector{SVector{2,T}},normal::Vector{SVector{2,T}},s::Vector{T},ds::Vector{T},w::Vector{T},w_n::Vector{T},curvature::Vector{T},xy_int::Vector{SVector{2,T}},shift_x::T,shift_y::T) where {T<:Real}
    return BoundaryPoints(xy;normal=normal,s=s,ds=ds,w=w,w_n=w_n,curvature=curvature,xy_int=xy_int,shift_x=shift_x,shift_y=shift_y)
end

"""
    BoundaryPoints(xy,tangent,tangent_2,ts,tphys,ws,ws_der,s,ds,compid,is_periodic,xL,xR,tL,tR) → bp::BoundaryPoints

Constructs a parametrized boundary discretization from sampled points, first and
second derivatives of the parametrization, computational nodes and quadrature
weights.

## Description
The outward unit normals are computed directly from `tangent` according to

    n = (t_y,-t_x)/|t|,

and the local arc-length coordinates are reconstructed from `ds`, with the first
boundary point assigned `s = 0`.

This constructor is used by boundary-integral discretizations that require the
underlying parametrization in addition to the physical boundary data.

## Returns
* `bp`: A [`BoundaryPoints`](@ref) instance with both physical and parametric boundary data populated.
"""
function BoundaryPoints(xy::Vector{SVector{2,T}},tangent::Vector{SVector{2,T}},tangent_2::Vector{SVector{2,T}},ts::Vector{T},tphys::Vector{T},ws::Vector{T},ws_der::Vector{T},s::Vector{T},ds::Vector{T},compid::Int,is_periodic::Bool,xL::SVector{2,T},xR::SVector{2,T},tL::SVector{2,T},tR::SVector{2,T}) where {T<:Real}
    n=length(xy)
    normal=Vector{SVector{2,T}}(undef,n)
    @inbounds for i in eachindex(tangent)
        tx,ty=tangent[i]
        sp=hypot(tx,ty)
        normal[i]=SVector{2,T}(ty/sp,-tx/sp)
    end
    return BoundaryPoints(xy;normal=normal,s=s,ds=ds,tangent=tangent,tangent_2=tangent_2,ts=ts,tphys=tphys,ws=ws,ws_der=ws_der,compid=compid,is_periodic=is_periodic,xL=xL,xR=xR,tL=tL,tR=tR)
end

"""
    length(pts::BoundaryPoints) → n::Int

Returns the number of sampled boundary points, `n = length(pts.xy)`.
"""
Base.length(pts::BoundaryPoints)=length(pts.xy)

"""
    isempty(pts::BoundaryPoints) → flag::Bool

Returns `true` if `pts` contains no sampled boundary points.
"""
Base.isempty(pts::BoundaryPoints)=isempty(pts.xy)

"""
    boundary_matrix_size(pts::BoundaryPoints) → N::Int

Returns the number of boundary degrees of freedom represented by `pts`.
"""
@inline boundary_matrix_size(pts::BoundaryPoints)=length(pts.xy)

"""
    boundary_matrix_size(pts::Vector{BoundaryPoints{T}}) where {T<:Real} → N::Int

Returns the total number of boundary degrees of freedom over all boundary
components in `pts`.
"""
function boundary_matrix_size(pts::Vector{BoundaryPoints{T}}) where {T<:Real}
    return sum(length,pts)
end

"""
    boundary_s(pts::BoundaryPoints) → s

Returns the stored arc-length coordinates of the boundary discretization.
"""
@inline boundary_s(pts::BoundaryPoints)=pts.s

"""
    boundary_s(pts::Vector{BoundaryPoints{T}}) where {T<:Real} → s::Vector{T}

Returns continuous arc-length coordinates for a vector of boundary components.

## Description
Each component stores its own local arc-length coordinates beginning at zero.
The components are shifted consecutively by their accumulated lengths so that the
returned vector defines a single continuous arc-length coordinate over the
concatenated boundary.
"""
function boundary_s(pts::Vector{BoundaryPoints{T}}) where {T<:Real}
    isempty(pts)&&return T[]
    s=T[]
    sizehint!(s,sum(length(p.s) for p in pts))
    soff=zero(T)
    for p in pts
        append!(s,p.s.+soff)
        soff+=sum(p.ds)
    end
    return s
end

"""
    component_offsets(pts::Vector{BoundaryPoints{T}}) where {T<:Real} → offs::Vector{Int}

Returns the starting indices of the boundary components in the flattened boundary
discretization.

## Description
For components containing `N₁,N₂,...,Nₘ` points, the returned offsets are

    [1,1+N₁,1+N₁+N₂,...,1+N₁+⋯+Nₘ].

Thus the points belonging to component `a` occupy

    offs[a]:offs[a+1]-1.
"""
function component_offsets(pts::Vector{BoundaryPoints{T}}) where {T<:Real}
    offs=Vector{Int}(undef,length(pts)+1)
    offs[1]=1
    @inbounds for i in eachindex(pts)
        offs[i+1]=offs[i]+length(pts[i])
    end
    return offs
end

"""
    component_offsets(pts::BoundaryPoints) → offs::Vector{Int}

Returns the component offsets for a single boundary component.
"""
@inline component_offsets(pts::BoundaryPoints)=[1,length(pts)+1]

"""
    _determine_bp_sizes(curves,bs,k) → Ns::Vector{Int64}

Determines the number of boundary points assigned to each curve.

## Description
For curve `i` with length `Lᵢ`, boundary sampling factor `bs[i]` and wavenumber
`k`, the number of points is estimated as

    Nᵢ = max(20,round(Int,k*Lᵢ*bs[i]/2π)).

The lower bound of `20` prevents very short boundary pieces from being sampled
with too few points.
"""
function _determine_bp_sizes(curves,bs,k)
    Ns=Vector{Int64}(undef,length(curves))
    @inbounds for i in eachindex(curves)
        Ns[i]=max(20,round(Int,k*curves[i].length*bs[i]/two_pi))
    end
    return Ns
end

"""
    boundary_coords(billiard::Bi,samplers::Vector{AbsSampler},Ns::Vector{Int64}) where {Bi<:BilliardGeometry.AbsBilliard} → bp::BoundaryPoints

Samples the complete physical boundary of `billiard` and returns the resulting
boundary discretization.

## Description
The physical boundary consists of curves carrying either `SpecularReflection` or
`QuantumSolverIgnore` boundary conditions. These curves are obtained from
`get_all_curves(billiard)` and filtered accordingly.

Each curve is sampled independently using the corresponding entry of `samplers`
and `Ns`. The local arc-length coordinate of each curve is shifted by the
cumulative length of all preceding curves so that `s` is continuous over the
concatenated physical boundary.

## Arguments
* `billiard`: The billiard whose physical boundary is sampled.
* `samplers`: One boundary sampler for each physical boundary curve.
* `Ns`: Number of sampling points assigned to each physical boundary curve.

## Returns
* `bp`: A [`BoundaryPoints`](@ref) instance with `xy`, `normal`, `s` and `ds` populated.
"""
function boundary_coords(billiard::Bi,samplers::Vector{<:BilliardGeometry.AbsSampler},Ns::Vector{Int64}) where {Bi<:BilliardGeometry.AbsBilliard}
    curves=filter(crv->crv.bc isa BilliardGeometry.SpecularReflection||crv.bc isa BilliardGeometry.QuantumSolverIgnore,BilliardGeometry.get_all_curves(billiard))
    T=typeof(curves[1].length)
    M=length(curves)
    length(samplers)==M||throw(DimensionMismatch("Expected $M samplers, received $(length(samplers))"))
    length(Ns)==M||throw(DimensionMismatch("Expected $M point counts, received $(length(Ns))"))
    xy_all=Vector{Vector{SVector{2,T}}}(undef,M)
    normal_all=Vector{Vector{SVector{2,T}}}(undef,M)
    s_all=Vector{Vector{T}}(undef,M)
    ds_all=Vector{Vector{T}}(undef,M)
    L0=zero(T)
    @inbounds for i in eachindex(curves)
        xy,normal,s,ds=boundary_coords(curves[i],samplers[i],Ns[i])
        xy_all[i]=xy
        normal_all[i]=normal
        s_all[i]=s.+L0
        ds_all[i]=ds
        L0+=curves[i].length
    end
    return BoundaryPoints(vcat(xy_all...);normal=vcat(normal_all...),s=vcat(s_all...),ds=vcat(ds_all...))
end

"""
    boundary_coords(billiard::Bi,sampler::BilliardGeometry.FourierNodes,N) where {Bi<:BilliardGeometry.AbsBilliard} → bp::BoundaryPoints

Samples the complete physical boundary of `billiard` using a global
[`FourierNodes`](@ref) discretization.

## Description
The physical boundary consists of curves carrying either `SpecularReflection` or
`QuantumSolverIgnore` boundary conditions. `sample_points` distributes the total
number of nodes `N` over these curves and returns one parameter grid and weight
vector for each component.

The local arc-length coordinate of every curve is shifted by the cumulative
length of the preceding curves so that the final `s` coordinate is continuous
over the complete physical boundary.

## Arguments
* `billiard`: The billiard whose physical boundary is sampled.
* `sampler`: Fourier-node sampling strategy.
* `N`: Total number of boundary sampling points.

## Returns
* `bp`: A [`BoundaryPoints`](@ref) instance with `xy`, `normal`, `s` and `ds` populated.
"""
function boundary_coords(billiard::Bi,sampler::BilliardGeometry.FourierNodes,N) where {Bi<:BilliardGeometry.AbsBilliard}
    curves=filter(crv->crv.bc isa BilliardGeometry.SpecularReflection||crv.bc isa BilliardGeometry.QuantumSolverIgnore,BilliardGeometry.get_all_curves(billiard))
    T=typeof(curves[1].length)
    M=length(curves)
    ts,dts=sample_points(sampler,N)
    length(ts)==M||throw(DimensionMismatch("Fourier sampler returned $(length(ts)) curve grids, expected $M"))
    xy_all=Vector{Vector{SVector{2,T}}}(undef,M)
    normal_all=Vector{Vector{SVector{2,T}}}(undef,M)
    s_all=Vector{Vector{T}}(undef,M)
    ds_all=Vector{Vector{T}}(undef,M)
    L0=zero(T)
    @inbounds for i in eachindex(curves)
        xy,normal,s,ds=boundary_coords(curves[i],ts[i],dts[i])
        xy_all[i]=xy
        normal_all[i]=normal
        s_all[i]=s.+L0
        ds_all[i]=ds
        L0+=curves[i].length
    end
    return BoundaryPoints(vcat(xy_all...);normal=vcat(normal_all...),s=vcat(s_all...),ds=vcat(ds_all...))
end

"""
    points_in_billiard(pts,billiard)

Returns the interior-membership mask of `pts` with respect to `billiard`.

## Description
Interior testing is delegated directly to `BilliardGeometry.is_inside`, avoiding
a separate polygonal approximation of the billiard boundary.
"""
@inline points_in_billiard(pts,billiard)=BilliardGeometry.is_inside(billiard,pts)

"""
    kress_R_even!(R0::AbstractMatrix{T}) where {T<:Real}

Constructs the periodic Kress logarithmic correction matrix for an even number of
nodes.

## Description
For `N = 2n`, the first column is computed spectrally using an inverse FFT of the
Fourier coefficients of the logarithmic kernel. The Nyquist correction term is
then added explicitly. The remaining columns are generated by circulant shifts.

The matrix corresponds to the periodic logarithmic kernel

    log(4sin²((t-τ)/2)).

## Arguments
* `R0`: Square matrix overwritten with the Kress correction matrix.

## Returns
* `nothing`.
"""
function kress_R_even!(R0::AbstractMatrix{T}) where {T<:Real}
    # Provides kress_R! to compute the circulant R matrix for the Kress method. kress_R! uses the FFT to compute the matrix efficiently, while kress_R! with ts computes it using a direct summation approach. Both functions modify the input matrix R0 in place.
    # Ref: Kress, R., Boundary integral equations in time-harmonic acoustic scattering. Mathematics Comput. Modelling Vol 15, pp. 229-243). Pergamon Press, 1991, GB.
    # Alex Barnett's code via ifft to get the circulant vector kernel and construct the circulant with circshift.
    N=size(R0,1)
    n=N÷2 # integer division
    a=zeros(Complex{T},N) #  build the spectral vector a (first col)
    @inbounds for m in 1:(n-1)
        a[m+1]=1/m     # positive freq
        a[N-m+1]=1/m     # negative freq
    end # leave a[n+1] == 0  (no 1/n term)
    rjn=real(FFTW.ifft(a)) # inverse FFT → rjn[j] = (2/N)*∑_{m=1..n-1} (1/m) cos(2π m (j-1)/N)
    ks=0:(N-1) # build the first column, adding the “alternating” correction
    alt=(-1).^ks # alt[j+1] = (-1)^j
    @. R0[:,1]=-two_pi*rjn-(2*two_pi/(N^2))*alt # R0[:,1] = -2π*rjn .- (4π/N^2)*alt, first col is ref
    @inbounds for j in 2:N # fill out the rest circulantly:
        @views R0[:,j].=circshift(R0[:,j-1],1) # shift by +1 wrt previous column
    end
    return nothing
end

"""
    kress_R_odd!(R0::AbstractMatrix{T}) where {T<:Real}

Constructs the periodic Kress logarithmic correction matrix for an odd number of
nodes.

## Description
Because an odd periodic grid has no Nyquist frequency, the first column is given
by the symmetric Fourier sum alone. The remaining columns are generated by
circulant shifts.

## Arguments
* `R0`: Square matrix overwritten with the Kress correction matrix.

## Returns
* `nothing`.
"""
function kress_R_odd!(R0::AbstractMatrix{T}) where {T<:Real}
    # This version of kress_R! computes the R matrix for the odd case (2n-1 points) where the Nyquist frequency is not included, so we only have m=1,...,n-1 positive and negative frequencies.
    # The first column is built using the same FFT approach, but with the appropriate range of m. The rest of the matrix is filled circulantly as before.
    N=size(R0,1)
    n=(N-1)÷2
    a=zeros(Complex{T},N)   # spectral first column
    for m in 1:n
        a[m+1]=1/m   # positive freq
        a[N-m+1]=1/m   # negative freq
    end
    rjn=real(FFTW.ifft(a)) # gives (1/N) * sum_m a_m exp(2πimj/N)
    @. R0[:,1]= -two_pi*rjn
    for j in 2:N
        @views R0[:,j].=circshift(R0[:,j-1],1)
    end
    return nothing
end

"""
    kress_R!(R0::AbstractMatrix{T}) where {T<:Real}

Constructs the periodic Kress logarithmic correction matrix in place.

## Description
The implementation automatically dispatches to [`kress_R_even!`](@ref) or
[`kress_R_odd!`](@ref) according to the parity of the matrix dimension.

## Returns
* `nothing`.
"""
function kress_R!(R0::AbstractMatrix{T}) where {T<:Real}
    iseven(size(R0,1)) ? kress_R_even!(R0) : kress_R_odd!(R0)
    return nothing
end

"""
    BoundaryPanelArrays{T}

Stores one-dimensional coordinate arrays extracted from a parametrized
[`BoundaryPoints`](@ref) discretization.

## Attributes
* `X`: x coordinates of the boundary nodes.
* `Y`: y coordinates of the boundary nodes.
* `dX`: x components of the parametrization derivative.
* `dY`: y components of the parametrization derivative.
* `speed`: Parametrization speed `|γ'(t)|` at each node.
"""
struct BoundaryPanelArrays{T<:Real}
    X::Vector{T}
    Y::Vector{T}
    dX::Vector{T}
    dY::Vector{T}
    speed::Vector{T}
end

"""
    _boundary_panel_arrays_cache(pts::BoundaryPoints{T}) where {T<:Real} → cache::BoundaryPanelArrays{T}

Extracts coordinate, tangent and speed arrays from a parametrized boundary
discretization.

## Returns
* `cache`: A [`BoundaryPanelArrays`](@ref) instance containing the flattened geometry arrays.
"""
@inline function _boundary_panel_arrays_cache(pts::BoundaryPoints{T}) where {T<:Real}
    X=getindex.(pts.xy,1)
    Y=getindex.(pts.xy,2)
    dX=getindex.(pts.tangent,1)
    dY=getindex.(pts.tangent,2)
    speed=@. hypot(dX,dY)
    return BoundaryPanelArrays(X,Y,dX,dY,speed)
end

"""
    _boundary_components(boundary)

Normalizes a boundary representation into a vector of connected components.

## Description
If `boundary` is already represented as a vector of component vectors, it is
returned unchanged. Otherwise each supplied curve is interpreted as an
independent component.

## Returns
* Vector of boundary components.
"""
function _boundary_components(boundary)
    isempty(boundary)&&throw(ArgumentError("Boundary cannot be empty"))
    return boundary[1] isa AbstractVector ? boundary : [[crv] for crv in boundary]
end

"""
    component_lengths(comp::Vector)

Computes segment lengths, cumulative segment lengths and the total length of one
composite boundary component.

## Returns
* `lens`: Length of each curve segment.
* `cum`: Cumulative lengths beginning at zero.
* `Ltot`: Total component length.
"""
function component_lengths(comp::Vector)
    lens=[crv.length for crv in comp]
    cum=Vector{eltype(lens)}(undef,length(lens)+1)
    cum[1]=zero(eltype(lens))
    @inbounds for j in eachindex(lens)
        cum[j+1]=cum[j]+lens[j]
    end
    return lens,cum,cum[end]
end

"""
    _unit_tangent_at_start(crv,::Type{T}) where {T<:Real}

Returns the unit tangent vector at the beginning of `crv`.
"""
@inline function _unit_tangent_at_start(crv,::Type{T}) where {T<:Real}
    v=SVector{2,T}(tangent(crv,zero(T)))
    return v/hypot(v[1],v[2])
end

"""
    _unit_tangent_at_end(crv,::Type{T}) where {T<:Real}

Returns the unit tangent vector at the end of `crv`.
"""
@inline function _unit_tangent_at_end(crv,::Type{T}) where {T<:Real}
    v=SVector{2,T}(tangent(crv,one(T)))
    return v/hypot(v[1],v[2])
end

"""
    _junction_angle(cleft,cright,::Type{T}) where {T<:Real} → angle::T

Returns the unsigned turning angle between the tangent leaving `cleft` and the
tangent entering `cright`.

The result lies in `[0,π]`.
"""
@inline function _junction_angle(cleft,cright,::Type{T}) where {T<:Real}
    tL=_unit_tangent_at_end(cleft,T)
    tR=_unit_tangent_at_start(cright,T)
    cr=tL[1]*tR[2]-tL[2]*tR[1]
    dt=clamp(tL[1]*tR[1]+tL[2]*tR[2],-one(T),one(T))
    return atan(abs(cr),dt)
end

"""
    _is_true_corner(cleft,cright,::Type{T};angle_tol=T(1e-8)) where {T<:Real} → flag::Bool

Returns `true` when the junction between `cleft` and `cright` has a tangent
discontinuity larger than `angle_tol`.
"""
@inline function _is_true_corner(cleft,cright,::Type{T};angle_tol=T(1e-8)) where {T<:Real}
    return _junction_angle(cleft,cright,T)>angle_tol
end

"""
    _component_corner_locations(::Type{T},comp::Vector;angle_tol=T(1e-8)) where {T<:Real} → corners::Vector{T}

Returns the locations of true corners of a composite boundary component in the
global periodic parameter `σ ∈ [0,2π)`.

## Description
Smooth joins are ignored. If the periodic seam between the final and first
segments is a true corner, it is represented by `σ = 0`.
"""
function _component_corner_locations(::Type{T},comp::Vector;angle_tol=T(1e-8)) where {T<:Real}
    _,cum,Ltot=component_lengths(comp)
    corners=T[]
    m=length(comp)
    _is_true_corner(comp[end],comp[1],T;angle_tol=angle_tol)&&push!(corners,zero(T))
    @inbounds for j in 1:m-1
        _is_true_corner(comp[j],comp[j+1],T;angle_tol=angle_tol)&&push!(corners,T(two_pi)*cum[j+1]/Ltot)
    end
    return corners
end

"""
    print_component_junctions(comp::Vector;T=Float64,angle_tol=1e-8)

Prints diagnostic information for all joins of a composite boundary component.

For every junction, the global periodic parameter, tangent discontinuity angle
and corner classification are displayed.
"""
function print_component_junctions(comp::Vector;T=Float64,angle_tol=1e-8)
    _,cum,Ltot=component_lengths(comp)
    m=length(comp)
    println("junction diagnostics:")
    @inbounds for j in 1:m
        jr=j==m ? 1 : j+1
        σ=j==m ? zero(T) : T(two_pi)*cum[j+1]/Ltot
        a=_junction_angle(comp[j],comp[jr],T)
        flag=a>T(angle_tol) ? "TRUE CORNER" : "smooth join"
        println("  $j -> $jr : σ = $σ, angle = $a, $flag")
    end
end

"""
    component_normals(pts::BoundaryPoints{T}) where {T<:Real}

Returns the Cartesian components of the stored outward normals together with the
parametrization speed.

## Returns
* `nx`: x components of the outward normals.
* `ny`: y components of the outward normals.
* `speed`: Parametrization speed `|γ'(t)|`; equal to one when tangent data are not stored.
"""
function component_normals(pts::BoundaryPoints{T}) where {T<:Real}
    length(pts.normal)==length(pts)||throw(ArgumentError("BoundaryPoints does not contain normal data"))
    nx=getindex.(pts.normal,1)
    ny=getindex.(pts.normal,2)
    if length(pts.tangent)==length(pts)
        tx=getindex.(pts.tangent,1)
        ty=getindex.(pts.tangent,2)
        speed=@. hypot(tx,ty)
    else
        speed=ones(T,length(pts))
    end
    return nx,ny,speed
end

"""
    flatten_boundary_components(comps::Vector{BoundaryPoints{T}}) where {T<:Real}

Flattens multiple boundary components into contiguous coordinate, normal and
quadrature arrays.

## Returns
A named tuple containing:
* `x`: x coordinates.
* `y`: y coordinates.
* `nx`: x components of outward normals.
* `ny`: y components of outward normals.
* `ds`: Arc-length quadrature elements.
* `offs`: Component offsets in the flattened arrays.
"""
function flatten_boundary_components(comps::Vector{BoundaryPoints{T}}) where {T<:Real}
    N=sum(length,comps)
    x=Vector{T}(undef,N)
    y=Vector{T}(undef,N)
    nx=Vector{T}(undef,N)
    ny=Vector{T}(undef,N)
    ds=Vector{T}(undef,N)
    offs=component_offsets(comps)
    p=1
    @inbounds for c in comps
        cnx,cny,_=component_normals(c)
        for j in eachindex(c.xy)
            q=c.xy[j]
            x[p]=q[1]
            y[p]=q[2]
            nx[p]=cnx[j]
            ny[p]=cny[j]
            ds[p]=c.ds[j]
            p+=1
        end
    end
    return (;x,y,nx,ny,ds,offs)
end

"""
    flatten_boundary_ds(comps::Vector{BoundaryPoints{T}}) where {T<:Real} → ds::Vector{T}

Concatenates the arc-length quadrature elements of all boundary components into a
single contiguous vector.
"""
function flatten_boundary_ds(comps::Vector{BoundaryPoints{T}}) where {T<:Real}
    ds=Vector{T}(undef,boundary_matrix_size(comps))
    p=1
    @inbounds for c in comps
        n=length(c.ds)
        ds[p:p+n-1].=c.ds
        p+=n
    end
    return ds
end

"""
    _global_t_to_segment_u(::Type{T},comp::Vector,t::T) where {T<:Real}

Maps a global periodic component parameter `t ∈ [0,2π)` to the corresponding
curve segment and local parameter.

## Returns
* `j`: Index of the active curve segment.
* `u`: Local curve parameter in `[0,1]`.
"""
function _global_t_to_segment_u(::Type{T},comp::Vector,t::T) where {T<:Real}
    lens,cum,Ltot=component_lengths(comp)
    s=(t/T(two_pi))*Ltot
    s>=Ltot&&return 1,zero(T)
    j=clamp(searchsortedlast(cum,s),1,length(comp))
    while j<length(comp)&&s>=cum[j+1]
        j+=1
    end
    slocal=s-cum[j]
    u=lens[j]==zero(T) ? zero(T) : slocal/lens[j]
    return j,clamp(u,zero(T),one(T))
end

"""
    _eval_composite_geom_global_t(::Type{T},comp::Vector,t::T) where {T<:Real}

Evaluates a composite boundary component at the global periodic parameter
`t ∈ [0,2π)`.

## Description
The global parameter is mapped to the appropriate curve segment and local
parameter. The first and second local curve derivatives are then transformed to
derivatives with respect to the global periodic parameter.

## Returns
* `xy`: Boundary point.
* `γt`: First derivative with respect to the global parameter.
* `γtt`: Second derivative with respect to the global parameter.
"""
function _eval_composite_geom_global_t(::Type{T},comp::Vector,t::T) where {T<:Real}
    lens,_,Ltot=component_lengths(comp)
    j,u=_global_t_to_segment_u(T,comp,t)
    crv=comp[j]
    xy=BilliardGeometry.curve(crv,u)
    du_dt=lens[j]==zero(T) ? zero(T) : Ltot/(T(two_pi)*lens[j])
    γu=tangent(crv,u)
    γuu=tangent_2(crv,u)
    γt=γu*du_dt
    γtt=γuu*du_dt^2
    return xy,γt,γtt
end

"""
    BoundaryGeomCache{T}

Stores pairwise and local geometric quantities reused during boundary-integral
matrix assembly.

## Attributes
* `R`: Pairwise distances between boundary nodes.
* `invR`: Pairwise inverse distances, with zero diagonal.
* `inner`: Tangential interaction term between source tangents and point differences.
* `logterm`: Periodic logarithmic kernel used in Kress splitting.
* `speed`: Parametrization speed at each boundary node.
* `kappa`: Scaled curvature term used in diagonal kernel limits.
* `original_ts`: Copy of the computational parameter nodes for corner-graded Kress discretizations.
"""
struct BoundaryGeomCache{T<:Real}
    R::Matrix{T}
    invR::Matrix{T}
    inner::Matrix{T}
    logterm::Matrix{T}
    speed::Vector{T}
    kappa::Vector{T}
    original_ts::Vector{T}
end

"""
    boundary_geom_cache(pts::BoundaryPoints{T},corner_kress::Bool=false) where {T<:Real} → cache::BoundaryGeomCache{T}

Constructs the geometric cache associated with a parametrized boundary
discretization.

## Description
The function precomputes pairwise distances, inverse distances, tangential
interaction factors, the periodic logarithmic Kress kernel, parametrization
speeds and the scaled curvature entering diagonal kernel limits.

If `corner_kress=true`, the computational parameter nodes are copied into
`original_ts` for use by the corner-graded Kress construction.

## Arguments
* `pts`: Parametrized boundary discretization.
* `corner_kress`: Whether the computational parameter grid should be retained for corner grading.

## Returns
* `cache`: A [`BoundaryGeomCache`](@ref) instance containing the precomputed geometry.
"""
function boundary_geom_cache(pts::BoundaryPoints{T},corner_kress::Bool=false) where {T<:Real}
    N=length(pts)
    length(pts.tangent)==N||throw(ArgumentError("BoundaryPoints does not contain tangent data"))
    length(pts.tangent_2)==N||throw(ArgumentError("BoundaryPoints does not contain tangent_2 data"))
    length(pts.ts)==N||throw(ArgumentError("BoundaryPoints does not contain ts data"))
    ts=pts.ts
    X=getindex.(pts.xy,1)
    Y=getindex.(pts.xy,2)
    dX=getindex.(pts.tangent,1)
    dY=getindex.(pts.tangent,2)
    ddX=getindex.(pts.tangent_2,1)
    ddY=getindex.(pts.tangent_2,2)
    ΔX=@. X-X'
    ΔY=@. Y-Y'
    R=hypot.(ΔX,ΔY)
    R[diagind(R)].=one(T)
    invR=inv.(R)
    invR[diagind(invR)].=zero(T)
    dX_row=reshape(dX,1,N)
    dY_row=reshape(dY,1,N)
    inner=@. dY_row*ΔX-dX_row*ΔY
    original_ts=corner_kress ? copy(ts) : T[]
    ΔT=ts.-ts'
    logterm=log.(4 .*sin.(ΔT./2).^2)
    logterm[diagind(logterm)].=zero(T)
    speed=@. hypot(dX,dY)
    κnum=-(dX.*ddY.-dY.*ddX)
    κden=dX.^2 .+dY.^2
    kappa=inv_two_pi.*(κnum./κden)
    return BoundaryGeomCache(R,invR,inner,logterm,speed,kappa,original_ts)
end