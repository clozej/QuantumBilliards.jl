"""
    MultiCornerAdaptedFourierBessel{T,B,Sy} <: AbsBasis

Direct-sum basis formed from several [`CornerAdaptedFourierBessel`](@ref) structs
centered at different corners.

If the individual blocks are

    B₁, B₂, ..., Bₘ,

with dimensions `d₁,d₂,...,dₘ`, the global basis is their ordered direct sum

    B = B₁ ⊕ B₂ ⊕ ... ⊕ Bₘ,

with total dimension

    dim = d₁ + d₂ + ... + dₘ.

Global basis indices run consecutively through the blocks. Thus block `j`
occupies the global index range

    offsets[j]+1 : offsets[j+1].

The basis also stores a common scaling origin used by the
Vergini-Saraceno method. This origin is independent of the local center of each
CAFB block. For a basis function `φ` belonging to any block, the MCAFB
[`dk_fun`](@ref) evaluates the common-dilation derivative

    ∂ₖᴠˢ φ(x) = ((x-c₀)⋅∇φ(x))/k,

where `c₀ = scaling_origin`.

For a block centered at `cⱼ`, this differs in general from its ordinary
wavenumber derivative `∂φ/∂k`. The common origin is required so that every block
participates in the same Vergini-Saraceno dilation.

## Arguments
- `T`: Real scalar type used by the basis.
- `B`: Tuple type containing the individual CAFB blocks.
- `Sy`: Type of the optional symmetry metadata.

## Returns
A `MultiCornerAdaptedFourierBessel` basis representing the direct sum of all
stored CAFB blocks.

## Fields
- `blocks`: Tuple containing the individual corner-adapted Fourier-Bessel
  blocks.
- `dim`: Total number of basis functions.
- `offsets`: Cumulative block dimensions. Block `j` occupies
  `offsets[j]+1:offsets[j+1]`.
- `weights`: Relative block dimensions used when automatically redistributing a
  requested total dimension during [`resize_basis`](@ref).
- `symmetries`: Optional symmetry metadata associated with the combined basis.
- `scaling_origin`: Common Cartesian origin `c₀` used for Vergini-Saraceno
  scaling derivatives.

## Notes
The boundary weights used by a Vergini-Saraceno solver must be constructed with
the same `scaling_origin`. In particular, the scaling factor

    rₙ = (x-c₀)⋅n

must use exactly the same `c₀` as the MCAFB `dk_fun`.
"""
struct MultiCornerAdaptedFourierBessel{T<:Real,B<:Tuple,Sy} <: AbsBasis
    blocks::B
    dim::Int
    offsets::Vector{Int}
    weights::Vector{T}
    symmetries::Sy
    scaling_origin::SVector{2,T}
end

# check for DM and VS since the origins need to mathc with evalute_points and the one used in dk_fun
@inline basis_origin_check(::AbsBasis,::Type{T}) where {T<:Real}=SVector{2,T}(zero(T),zero(T))
@inline basis_origin_check(basis::CornerAdaptedFourierBessel,::Type{T}) where {T<:Real}=SVector{2,T}(basis.cs.origin)
@inline basis_origin_check(basis::MultiCornerAdaptedFourierBessel,::Type{T}) where {T<:Real}=SVector{2,T}(basis.scaling_origin)

"""
    MultiCornerAdaptedFourierBessel(blocks::Vararg{CornerAdaptedFourierBessel,N};symmetries=nothing,scaling_origin=nothing) where {N}

Construct a [`MultiCornerAdaptedFourierBessel`](@ref) from existing
corner-adapted Fourier-Bessel blocks.

The blocks are concatenated in the supplied order. Their dimensions determine
the global basis indexing and the initial resizing weights.

## Arguments
- `blocks`: One or more [`CornerAdaptedFourierBessel`](@ref) blocks. All blocks
  must use the same underlying real type and must have positive dimension.
- `symmetries`: Optional symmetry information attached to the combined basis.
- `scaling_origin`: Common Cartesian dilation origin used by the
  Vergini-Saraceno derivative. If omitted, the origin `(0,0)` is used.

## Returns
A [`MultiCornerAdaptedFourierBessel`](@ref) whose total dimension is the sum of
the dimensions of the supplied blocks.

## Notes
`scaling_origin` is not a CAFB center. Each block retains its own local polar
coordinate system, while `scaling_origin` defines the single global dilation
generator used by Vergini-Saraceno:

    (x-c₀)⋅∇.

Consequently, changing `scaling_origin` changes [`dk_fun`](@ref) for the
multi-corner basis but does not change the basis functions themselves.
"""
function MultiCornerAdaptedFourierBessel(blocks::Vararg{CornerAdaptedFourierBessel,N};symmetries=nothing,scaling_origin=nothing) where {N}
    N>0||throw(ArgumentError("at least one CAFB block is required"))
    T=typeof(blocks[1].corner_angle)
    all(b->typeof(b.corner_angle)===T,blocks)||throw(ArgumentError("all CAFB blocks must use the same real type"))
    dims=Int[b.dim for b in blocks]
    all(>(0),dims)||throw(ArgumentError("all CAFB blocks must have positive dimension"))
    offsets=Vector{Int}(undef,N+1)
    offsets[1]=0
    @inbounds for j=1:N
        offsets[j+1]=offsets[j]+dims[j]
    end
    dim=offsets[end]
    weights=T.(dims)./T(dim)
    origin=scaling_origin===nothing ? SVector{2,T}(zero(T),zero(T)) : SVector{2,T}(scaling_origin)
    return MultiCornerAdaptedFourierBessel{T,typeof(blocks),typeof(symmetries)}(blocks,dim,offsets,weights,symmetries,origin)
end

"""
    MultiCornerAdaptedFourierBessel(blocks::Tuple;symmetries=nothing,scaling_origin=nothing)

Tuple-based convenience constructor for [`MultiCornerAdaptedFourierBessel`](@ref).

Equivalent to splatting `blocks` into the vararg constructor.

## Arguments
- `blocks`: Tuple of [`CornerAdaptedFourierBessel`](@ref) blocks.
- `symmetries`: Optional symmetry information attached to the combined basis.
- `scaling_origin`: Common Cartesian dilation origin used by the
  Vergini-Saraceno derivative.

## Returns
A [`MultiCornerAdaptedFourierBessel`](@ref) constructed from the supplied tuple
of blocks.
"""
MultiCornerAdaptedFourierBessel(blocks::Tuple;symmetries=nothing,scaling_origin=nothing)=MultiCornerAdaptedFourierBessel(blocks...;symmetries=symmetries,scaling_origin=scaling_origin)

"""
    length(basis::MultiCornerAdaptedFourierBessel)

Return the total number of basis functions in `basis`.

## Arguments
- `basis`: Multi-corner adapted Fourier-Bessel basis.

## Returns
- `dim::Int`: Total number of basis functions.
"""
Base.length(basis::MultiCornerAdaptedFourierBessel)=basis.dim

"""
    _global_to_local_basis_index(basis::MultiCornerAdaptedFourierBessel,i::Int)

Map global MCAFB basis index `i` to its containing CAFB block and the
corresponding block-local index.

If block `j` occupies

    offsets[j]+1 : offsets[j+1],

then the returned local index is

    local_index = i-offsets[j].

Throws `BoundsError` when `i` lies outside `1:basis.dim`.

## Arguments
- `basis`: Multi-corner adapted Fourier-Bessel basis.
- `i::Int`: Global basis index.

## Returns
- `block::Int`: Index of the CAFB block containing global basis function `i`.
- `local_index::Int`: Corresponding basis index inside that block.
"""
@inline function _global_to_local_basis_index(basis::MultiCornerAdaptedFourierBessel,i::Int)
    1<=i<=basis.dim||throw(BoundsError(basis,i))
    block=searchsortedlast(basis.offsets,i-1)
    local_index=i-basis.offsets[block]
    return block,local_index
end

"""
    _distribute_basis_dimensions(weights::AbstractVector{T},dim::Int) where {T<:Real}

Distribute a requested total basis dimension `dim` among several blocks
according to `weights`.

Every block receives at least one basis function. The remaining
`dim-length(weights)` functions are distributed proportionally to `weights`.
Integer rounding is performed using the largest-remainder rule so that

    sum(dims) == dim

exactly.

## Arguments
- `weights::AbstractVector{T}`: Relative block weights, normally inherited from
  the original block dimensions of a [`MultiCornerAdaptedFourierBessel`](@ref).
- `dim::Int`: Requested total basis dimension.

## Returns
- `dims::Vector{Int}`: Dimension assigned to each block.

## Throws
`ArgumentError` if `dim` is smaller than the number of blocks.
"""
function _distribute_basis_dimensions(weights::AbstractVector{T},dim::Int) where {T<:Real}
    nblocks=length(weights)
    dim>=nblocks||throw(ArgumentError("total dimension $dim must be at least the number of blocks $nblocks"))
    remaining=dim-nblocks
    target=remaining.*weights
    extra=floor.(Int,target)
    dims=extra.+1
    leftover=dim-sum(dims)
    if leftover>0
        fractions=target.-extra
        order=sortperm(fractions;rev=true)
        @inbounds for j=1:leftover
            dims[order[j]]+=1
        end
    end
    return dims
end

"""
    resize_basis(basis::MultiCornerAdaptedFourierBessel,billiard::Bi,dim::Int,k) where {Bi<:AbsBilliard}

Resize a multi-corner basis to total dimension `dim`.

The requested dimension is distributed among the CAFB blocks according to the
stored relative `basis.weights`, with at least one function assigned to every
block. Each block is then resized independently using its own [`resize_basis`](@ref)
method.

The original relative weights are preserved in the returned MCAFB basis so that
successive automatic resizings continue to use the same allocation proportions.

The common `scaling_origin` and `symmetries` metadata are preserved.

## Arguments
- `basis`: Multi-corner adapted Fourier-Bessel basis to resize.
- `billiard::Bi`: Billiard geometry passed to the individual block resizing
  methods.
- `dim::Int`: Requested total basis dimension.
- `k`: Wavenumber used by the underlying basis resizing procedure.

## Returns
- `newbasis::MultiCornerAdaptedFourierBessel`: Resized basis with total
  dimension `dim` and the original relative block weights.
"""
function resize_basis(basis::MultiCornerAdaptedFourierBessel,billiard::Bi,dim::Int,k) where {Bi<:AbsBilliard}
    basis.dim==dim&&return basis
    dims=_distribute_basis_dimensions(basis.weights,dim)
    blocks=ntuple(j->resize_basis(basis.blocks[j],billiard,dims[j],k),length(basis.blocks))
    newbasis=MultiCornerAdaptedFourierBessel(blocks;symmetries=basis.symmetries,scaling_origin=basis.scaling_origin)
    return MultiCornerAdaptedFourierBessel{typeof(newbasis.weights[1]),typeof(newbasis.blocks),typeof(newbasis.symmetries)}(newbasis.blocks,newbasis.dim,newbasis.offsets,copy(basis.weights),newbasis.symmetries,newbasis.scaling_origin)
end

"""
    resize_basis(basis::MultiCornerAdaptedFourierBessel,billiard::Bi,dims::AbstractVector{<:Integer},k) where {Bi<:AbsBilliard}

Resize each CAFB block to an explicitly specified dimension.

`dims[j]` gives the requested dimension of block `j`. The number of supplied
dimensions must equal the number of blocks and every value must be positive.

Unlike the total-dimension overload, this method defines a new block-size
distribution, so the returned MCAFB resizing weights are recomputed from the
requested dimensions.

The common `scaling_origin` and `symmetries` metadata are preserved.

## Arguments
- `basis`: Multi-corner adapted Fourier-Bessel basis to resize.
- `billiard::Bi`: Billiard geometry passed to the individual block resizing
  methods.
- `dims::AbstractVector{<:Integer}`: Requested dimension of each CAFB block.
- `k`: Wavenumber used by the underlying basis resizing procedure.

## Returns
- `newbasis::MultiCornerAdaptedFourierBessel`: Basis whose block dimensions are
  given by `dims`.
"""
function resize_basis(basis::MultiCornerAdaptedFourierBessel,billiard::Bi,dims::AbstractVector{<:Integer},k) where {Bi<:AbsBilliard}
    length(dims)==length(basis.blocks)||throw(DimensionMismatch("received $(length(dims)) dimensions for $(length(basis.blocks)) blocks"))
    all(>(0),dims)||throw(ArgumentError("all block dimensions must be positive"))
    blocks=ntuple(j->resize_basis(basis.blocks[j],billiard,Int(dims[j]),k),length(basis.blocks))
    return MultiCornerAdaptedFourierBessel(blocks;symmetries=basis.symmetries,scaling_origin=basis.scaling_origin)
end

"""
    basis_fun(basis::MultiCornerAdaptedFourierBessel,i::Int,k::T,pts::AbstractArray) where {T<:Real}

Evaluate global MCAFB basis function `i` at the points `pts`.

The global index is mapped to the corresponding CAFB block and block-local
basis index before evaluation.

## Arguments
- `basis`: Multi-corner adapted Fourier-Bessel basis.
- `i::Int`: Global basis index.
- `k::T`: Wavenumber.
- `pts::AbstractArray`: Cartesian evaluation points.

## Returns
- `B`: Values of global basis function `i` at `pts`.
"""
@inline function basis_fun(basis::MultiCornerAdaptedFourierBessel,i::Int,k::T,pts::AbstractArray) where {T<:Real}
    block,local_index=_global_to_local_basis_index(basis,i)
    return basis_fun(basis.blocks[block],local_index,k,pts)
end

"""
    dk_fun(basis::MultiCornerAdaptedFourierBessel,i::Int,k::T,pts::AbstractArray) where {T<:Real}

Evaluate the Vergini-Saraceno common-scaling derivative of global MCAFB basis
function `i`.

For the common scaling origin `c₀ = basis.scaling_origin`, this method returns

    ∂ₖᴠˢ φ(x) = ((x-c₀)⋅∇φ(x))/k.

This is the derivative associated with a global dilation about `c₀`, not in
general the ordinary partial derivative of the underlying CAFB block with
respect to `k`.

For a CAFB block centered at `cⱼ`,

    k ∂φ/∂k = (x-cⱼ)⋅∇φ,

whereas the common Vergini-Saraceno derivative satisfies

    k ∂ₖᴠˢφ
    = (x-c₀)⋅∇φ
    = k ∂φ/∂k + (cⱼ-c₀)⋅∇φ.

This correction is essential when several CAFB blocks have different centers.

## Arguments
- `basis`: Multi-corner adapted Fourier-Bessel basis.
- `i::Int`: Global basis index.
- `k::T`: Wavenumber.
- `pts::AbstractArray`: Cartesian evaluation points.

## Returns
- `dB_dk::Vector{T}`: Vergini-Saraceno common-scaling derivative evaluated at
  `pts`.

## Notes
A Vergini-Saraceno boundary discretization used with this basis must employ the
same origin `c₀` in its scaling weight

    rₙ = (x-c₀)⋅n.
"""
@inline function dk_fun(basis::MultiCornerAdaptedFourierBessel,i::Int,k::T,pts::AbstractArray) where {T<:Real}
    block,local_index=_global_to_local_basis_index(basis,i)
    dB_dx,dB_dy=gradient(basis.blocks[block],local_index,k,pts)
    M=length(pts)
    dB_dk=Vector{T}(undef,M)
    x0=basis.scaling_origin[1]
    y0=basis.scaling_origin[2]
    @inbounds @simd for j=1:M
        rx=pts[j][1]-x0
        ry=pts[j][2]-y0
        dB_dk[j]=(rx*dB_dx[j]+ry*dB_dy[j])/k
    end
    return dB_dk
end

"""
    gradient(basis::MultiCornerAdaptedFourierBessel,i::Int,k::T,pts::AbstractArray) where {T<:Real}

Evaluate the Cartesian gradient of global MCAFB basis function `i`.

## Arguments
- `basis`: Multi-corner adapted Fourier-Bessel basis.
- `i::Int`: Global basis index.
- `k::T`: Wavenumber.
- `pts::AbstractArray`: Cartesian evaluation points.

## Returns
- `dB_dx`: Derivative with respect to the global Cartesian coordinate `x`.
- `dB_dy`: Derivative with respect to the global Cartesian coordinate `y`.
"""
@inline function gradient(basis::MultiCornerAdaptedFourierBessel,i::Int,k::T,pts::AbstractArray) where {T<:Real}
    block,local_index=_global_to_local_basis_index(basis,i)
    return gradient(basis.blocks[block],local_index,k,pts)
end

"""
    basis_and_gradient(basis::MultiCornerAdaptedFourierBessel,i::Int,k::T,pts::AbstractArray) where {T<:Real}

Evaluate global MCAFB basis function `i` together with its Cartesian gradient.

## Arguments
- `basis`: Multi-corner adapted Fourier-Bessel basis.
- `i::Int`: Global basis index.
- `k::T`: Wavenumber.
- `pts::AbstractArray`: Cartesian evaluation points.

## Returns
- `B`: Basis-function values at `pts`.
- `dB_dx`: Derivative with respect to the global Cartesian coordinate `x`.
- `dB_dy`: Derivative with respect to the global Cartesian coordinate `y`.
"""
@inline function basis_and_gradient(basis::MultiCornerAdaptedFourierBessel,i::Int,k::T,pts::AbstractArray) where {T<:Real}
    block,local_index=_global_to_local_basis_index(basis,i)
    return basis_and_gradient(basis.blocks[block],local_index,k,pts)
end

"""
    _group_basis_indices_by_block(basis::MultiCornerAdaptedFourierBessel,indices)

Group requested global MCAFB indices according to their underlying CAFB block.

For every block `j`:

- `columns[j]` contains the output-matrix columns associated with that block,
- `local_indices[j]` contains the corresponding block-local basis indices.

This permits batched block evaluation while preserving the original ordering of
`indices` in the returned matrices.

## Arguments
- `basis`: Multi-corner adapted Fourier-Bessel basis.
- `indices`: Requested global basis indices.

## Returns
- `columns`: For each CAFB block, output-column positions associated with that
  block.
- `local_indices`: Corresponding block-local basis indices.
"""
function _group_basis_indices_by_block(basis::MultiCornerAdaptedFourierBessel,indices)
    nblocks=length(basis.blocks)
    columns=[Int[] for _=1:nblocks]
    local_indices=[Int[] for _=1:nblocks]
    @inbounds for (column,i) in enumerate(indices)
        block,local_index=_global_to_local_basis_index(basis,i)
        push!(columns[block],column)
        push!(local_indices[block],local_index)
    end
    return columns,local_indices
end

"""
    basis_fun(basis::MultiCornerAdaptedFourierBessel,indices::AbstractArray,k::T,pts::AbstractArray;multithreaded::Bool=true) where {T<:Real}

Evaluate several MCAFB basis functions at once.

If `M = length(pts)` and `N = length(indices)`, returns an `M × N` matrix `B`
with

    B[j,c] = φ_{indices[c]}(pts[j]).

Requested functions are grouped by CAFB block and evaluated in block batches.
The final column ordering agrees with the supplied `indices`.

## Arguments
- `basis`: Multi-corner adapted Fourier-Bessel basis.
- `indices::AbstractArray`: Requested global basis indices.
- `k::T`: Wavenumber.
- `pts::AbstractArray`: Cartesian evaluation points.
- `multithreaded::Bool`: Whether the underlying CAFB block evaluations may use
  multithreading.

## Returns
- `B::Matrix{T}`: Matrix of size `length(pts) × length(indices)` containing the
  requested basis-function values.
"""
@inline function basis_fun(basis::MultiCornerAdaptedFourierBessel,indices::AbstractArray,k::T,pts::AbstractArray;multithreaded::Bool=true) where {T<:Real}
    M=length(pts)
    N=length(indices)
    B=Matrix{T}(undef,M,N)
    columns,local_indices=_group_basis_indices_by_block(basis,indices)
    @inbounds for block in eachindex(basis.blocks)
        isempty(columns[block])&&continue
        block_values=basis_fun(basis.blocks[block],local_indices[block],k,pts;multithreaded=multithreaded)
        @views B[:,columns[block]].=block_values
    end
    return B
end

"""
    dk_fun(basis::MultiCornerAdaptedFourierBessel,indices::AbstractArray,k::T,pts::AbstractArray;multithreaded::Bool=true) where {T<:Real}

Evaluate the Vergini-Saraceno common-scaling derivatives of several MCAFB
basis functions.

For every requested basis function `φ` and common scaling origin
`c₀ = basis.scaling_origin`, the returned derivative is

    ∂ₖᴠˢ φ(x) = ((x-c₀)⋅∇φ(x))/k.

If `M = length(pts)` and `N = length(indices)`, the result is an `M × N`
matrix whose columns follow the ordering of `indices`.

The gradients are evaluated blockwise and then contracted with the same global
dilation vector `x-c₀` for every CAFB block.

## Arguments
- `basis`: Multi-corner adapted Fourier-Bessel basis.
- `indices::AbstractArray`: Requested global basis indices.
- `k::T`: Wavenumber.
- `pts::AbstractArray`: Cartesian evaluation points.
- `multithreaded::Bool`: Whether the underlying CAFB gradient evaluations may
  use multithreading.

## Returns
- `dB_dk::Matrix{T}`: Matrix of size
  `length(pts) × length(indices)` containing the common-scaling derivatives.

## Notes
This is intentionally different from applying the ordinary CAFB `dk_fun` to
each block independently. Independent block derivatives correspond to
dilations about different CAFB centers and therefore do not define a common
Vergini-Saraceno scaling transformation.
"""
function dk_fun(basis::MultiCornerAdaptedFourierBessel,indices::AbstractArray,k::T,pts::AbstractArray;multithreaded::Bool=true) where {T<:Real}
    M=length(pts)
    N=length(indices)
    dB_dk=Matrix{T}(undef,M,N)
    columns,local_indices=_group_basis_indices_by_block(basis,indices)
    x0=basis.scaling_origin[1]
    y0=basis.scaling_origin[2]
    @inbounds for block in eachindex(basis.blocks)
        isempty(columns[block])&&continue
        block_dx,block_dy=gradient(basis.blocks[block],local_indices[block],k,pts;multithreaded=multithreaded)
        cols=columns[block]
        @views block_dk=dB_dk[:,cols]
        @inbounds for c in eachindex(cols)
            @simd for j=1:M
                rx=pts[j][1]-x0
                ry=pts[j][2]-y0
                block_dk[j,c]=(rx*block_dx[j,c]+ry*block_dy[j,c])/k
            end
        end
    end
    return dB_dk
end

"""
    gradient(basis::MultiCornerAdaptedFourierBessel,indices::AbstractArray,k::T,pts::AbstractArray;multithreaded::Bool=true) where {T<:Real}

Evaluate the Cartesian gradients of several MCAFB basis functions.

If `M = length(pts)` and `N = length(indices)`, returns two `M × N` matrices

    dB_dx, dB_dy,

with

    dB_dx[j,c] = ∂ₓφ_{indices[c]}(pts[j]),
    dB_dy[j,c] = ∂ᵧφ_{indices[c]}(pts[j]).

The requested functions are grouped and evaluated blockwise while preserving
the original column ordering.

## Arguments
- `basis`: Multi-corner adapted Fourier-Bessel basis.
- `indices::AbstractArray`: Requested global basis indices.
- `k::T`: Wavenumber.
- `pts::AbstractArray`: Cartesian evaluation points.
- `multithreaded::Bool`: Whether the underlying CAFB gradient evaluations may
  use multithreading.

## Returns
- `dB_dx::Matrix{T}`: Derivatives with respect to the global Cartesian
  coordinate `x`.
- `dB_dy::Matrix{T}`: Derivatives with respect to the global Cartesian
  coordinate `y`.
"""
function gradient(basis::MultiCornerAdaptedFourierBessel,indices::AbstractArray,k::T,pts::AbstractArray;multithreaded::Bool=true) where {T<:Real}
    M=length(pts)
    N=length(indices)
    dB_dx=Matrix{T}(undef,M,N)
    dB_dy=Matrix{T}(undef,M,N)
    columns,local_indices=_group_basis_indices_by_block(basis,indices)
    @inbounds for block in eachindex(basis.blocks)
        isempty(columns[block])&&continue
        block_dx,block_dy=gradient(basis.blocks[block],local_indices[block],k,pts;multithreaded=multithreaded)
        @views dB_dx[:,columns[block]].=block_dx
        @views dB_dy[:,columns[block]].=block_dy
    end
    return dB_dx,dB_dy
end

"""
    basis_and_gradient(basis::MultiCornerAdaptedFourierBessel,indices::AbstractArray,k::T,pts::AbstractArray;multithreaded::Bool=true) where {T<:Real}

Evaluate several MCAFB basis functions together with their Cartesian gradients.

If `M = length(pts)` and `N = length(indices)`, returns

    B, dB_dx, dB_dy,

where each array has size `M × N` and

    B[j,c]     = φ_{indices[c]}(pts[j]),
    dB_dx[j,c] = ∂ₓφ_{indices[c]}(pts[j]),
    dB_dy[j,c] = ∂ᵧφ_{indices[c]}(pts[j]).

Evaluation is performed blockwise using the underlying CAFB implementations,
while the returned columns preserve the ordering of `indices`.

## Arguments
- `basis`: Multi-corner adapted Fourier-Bessel basis.
- `indices::AbstractArray`: Requested global basis indices.
- `k::T`: Wavenumber.
- `pts::AbstractArray`: Cartesian evaluation points.
- `multithreaded::Bool`: Whether the underlying CAFB evaluations may use
  multithreading.

## Returns
- `B::Matrix{T}`: Basis-function values.
- `dB_dx::Matrix{T}`: Derivatives with respect to the global Cartesian
  coordinate `x`.
- `dB_dy::Matrix{T}`: Derivatives with respect to the global Cartesian
  coordinate `y`.
"""
function basis_and_gradient(basis::MultiCornerAdaptedFourierBessel,indices::AbstractArray,k::T,pts::AbstractArray;multithreaded::Bool=true) where {T<:Real}
    M=length(pts)
    N=length(indices)
    B=Matrix{T}(undef,M,N)
    dB_dx=Matrix{T}(undef,M,N)
    dB_dy=Matrix{T}(undef,M,N)
    columns,local_indices=_group_basis_indices_by_block(basis,indices)
    @inbounds for block in eachindex(basis.blocks)
        isempty(columns[block])&&continue
        block_values,block_dx,block_dy=basis_and_gradient(basis.blocks[block],local_indices[block],k,pts;multithreaded=multithreaded)
        @views B[:,columns[block]].=block_values
        @views dB_dx[:,columns[block]].=block_dx
        @views dB_dy[:,columns[block]].=block_dy
    end
    return B,dB_dx,dB_dy
end