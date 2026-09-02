################################################################################
########################## SYMMETRY-REDUCED BOUNDARY MAPS ######################
################################################################################
# Boundary-integral solvers with discrete symmetries are assembled on a reduced
# fundamental domain rather than on the full boundary.
#
# If G is the symmetry group and b denotes one fundamental boundary node, the
# corresponding full-boundary orbit is
#
#     O_b={g·b : g∈G}.
#
# In a chosen irreducible representation, the boundary unknown satisfies
#
#     u(g·b)=χ(g)u(b),
#
# where χ(g) is the parity factor for reflections or the character factor for
# rotations. A full-boundary source sum can therefore be folded into one reduced
# source column:
#
#     K_red(i,b)=Σ_{g∈G} χ(g)K(i,g·b).
#
# The purpose of the orbit maps below is to encode this reduction once in exact
# integer index space. They provide both
#
#     fundamental index -> all full indices in its symmetry orbit
#
# and
#
#     full index -> corresponding fundamental index and irrep factor.
#
# The boundary discretizations are constructed to respect the symmetry exactly:
# the number of nodes is a multiple of |G|, the fundamental domain is the first
# contiguous block, and symmetry-related components have identical node counts
# and known ordering. Consequently no geometric point matching or floating-point
# tolerance is required; all symmetry actions are exact permutations of boundary
# indices.
#
# The same SymmetryOrbitMap representation is used by DLP, DLP-Kress and CFIE-Kress,
# so matrix sizing, reduced kernel assembly, and reconstruction can share one
# common symmetry interface.
################################################################################

@inline symmetry_order(::BilliardGeometry.XAxisReflection)=2
@inline symmetry_order(::BilliardGeometry.YAxisReflection)=2
@inline symmetry_order(::BilliardGeometry.XYAxisReflection)=4
@inline symmetry_order(::DiagonalReflection)=2
@inline symmetry_order(::AntiDiagonalReflection)=2
@inline symmetry_order(sym::BilliardGeometry.NFoldRotation)=sym.order
@inline symmetry_order(syms::AbstractVector{<:BilliardGeometry.NFoldRotation})=isempty(syms) ? 1 : syms[1].order

@inline symmetry_node_multiple(::BilliardGeometry.XAxisReflection)=4
@inline symmetry_node_multiple(::BilliardGeometry.YAxisReflection)=4
@inline symmetry_node_multiple(::BilliardGeometry.XYAxisReflection)=4
@inline symmetry_node_multiple(::DiagonalReflection)=8
@inline symmetry_node_multiple(::AntiDiagonalReflection)=8
@inline symmetry_node_multiple(sym::BilliardGeometry.CompositeReflection)=foldl(lcm,(symmetry_node_multiple(ref) for ref in sym.reflections);init=1)
@inline symmetry_node_multiple(sym::BilliardGeometry.NFoldRotation)=sym.order
@inline symmetry_node_multiple(syms::AbstractVector{<:BilliardGeometry.NFoldRotation})=isempty(syms) ? 1 : syms[1].order

"""
    SymmetryOrbitMap{T}
Store the exact correspondence between a symmetry-reduced boundary and the full
boundary discretization.
For a symmetry group `G` of order `ng`, each reduced boundary degree of freedom
represents one full-boundary orbit
    O_b={g·b : g∈G}.
The matrices `fund_to_full` and `fund_to_scale` describe this expansion directly:
    q = fund_to_full[g,b],
    χ = fund_to_scale[g,b],
meaning that group image `g` of reduced node `b` is full node `q` and
    u[q]=χ*u[Ifund[b]].
The inverse arrays `full_to_fund` and `full_to_scale` provide the corresponding
full-to-reduced lookup.
This representation is used when constructing symmetry-reduced boundary-integral
matrices, where a complete source orbit is folded into one reduced source column.
## Attributes
* `Ifund`: Full-boundary indices chosen as representatives of the fundamental domain.
* `full_to_fund`: Reduced-orbit index associated with every full-boundary node.
* `full_to_scale`: Irrep factor relating every full-boundary node to its representative.
* `fund_to_full`: Full-boundary indices in each orbit; rows are symmetry images and columns are reduced nodes.
* `fund_to_scale`: Irrep factors corresponding to `fund_to_full`.
"""
struct SymmetryOrbitMap{T<:Real}
    Ifund::Vector{Int}
    full_to_fund::Vector{Int}
    full_to_scale::Vector{Complex{T}}
    fund_to_full::Matrix{Int}
    fund_to_scale::Matrix{Complex{T}}
end
Base.length(orbits::SymmetryOrbitMap)=length(orbits.Ifund)
@inline fundamental_size(orbits::SymmetryOrbitMap)=length(orbits.Ifund)
@inline full_size(orbits::SymmetryOrbitMap)=length(orbits.full_to_fund)
@inline orbit_size(orbits::SymmetryOrbitMap)=size(orbits.fund_to_full,1)

"""
    symmetry_orbit(orbits,b)
Return the full-boundary indices and irrep factors represented by reduced node
`b`. If qs,χs=symmetry_orbit(orbits,b),
then the reduced boundary value `u_b` generates the full orbit according to u[qs[l]]=χs[l]*u_b.
The returned arrays are views into `orbits` and therefore allocate no copies.
"""
@inline function symmetry_orbit(orbits::SymmetryOrbitMap,b::Int)
    return @view(orbits.fund_to_full[:,b]),@view(orbits.fund_to_scale[:,b])
end

"""
    _build_periodic_symmetry_orbit_map(::Type{T},perms,scales)
Build a [`SymmetryOrbitMap`](@ref) from exact boundary-index permutations.
Each entry `perms[g]` represents one symmetry-group element: perms[g][q]
is the full-boundary node obtained by applying group element `g` to node `q`.
The corresponding representation factor is `scales[g]`.

The boundary discretization is assumed to have been constructed consistently
with the symmetry:

* `N` is divisible by the group order `ng=length(perms)`;
* every symmetry image is represented by an exact index permutation;
* every node belongs to a complete orbit of size `ng`.

Fundamental representatives are selected directly from the integer permutations
and therefore need not form a contiguous index block.

## Arguments
* `T`: Real scalar type used for the complex irrep factors.
* `perms`: Full-boundary index permutation for each group element.
* `scales`: Parity or character factor associated with each group element.
## Returns
A [`SymmetryOrbitMap{T}`](@ref).
"""
function _build_periodic_symmetry_orbit_map(::Type{T},perms::Vector{Vector{Int}},scales::Vector{Complex{T}}) where {T<:Real}
    ng=length(perms)
    ng>0||throw(ArgumentError("At least one symmetry permutation is required"))
    length(scales)==ng||throw(DimensionMismatch("Received $ng permutations but $(length(scales)) irrep factors"))
    N=length(perms[1])
    all(length(p)==N for p in perms)||throw(DimensionMismatch("All symmetry permutations must have length $N"))
    N%ng==0||throw(ArgumentError("Node count N=$N must be divisible by symmetry-group order $ng"))
    nf=N÷ng
    Ifund=Vector{Int}(undef,nf)
    fund_to_full=Matrix{Int}(undef,ng,nf)
    fund_to_scale=Matrix{Complex{T}}(undef,ng,nf)
    full_to_fund=Vector{Int}(undef,N)
    full_to_scale=Vector{Complex{T}}(undef,N)
    seen=falses(N)
    b=0
    @inbounds for q in 1:N
        seen[q]&&continue
        b+=1
        b<=nf||throw(ArgumentError("Symmetry permutations generate more than $nf boundary orbits"))
        Ifund[b]=q
        for g in 1:ng
            qi=perms[g][q]
            1<=qi<=N||throw(BoundsError(perms[g],qi))
            seen[qi]&&g>1&&throw(ArgumentError("Symmetry orbit of node $q has fewer than $ng distinct nodes"))
            χ=scales[g]
            fund_to_full[g,b]=qi
            fund_to_scale[g,b]=χ
            full_to_fund[qi]=b
            full_to_scale[qi]=χ
            seen[qi]=true
        end
    end
    b==nf||throw(ArgumentError("Expected $nf symmetry orbits, constructed $b"))
    return SymmetryOrbitMap{T}(Ifund,full_to_fund,full_to_scale,fund_to_full,fund_to_scale)
end

# Build the complete exact `Cₙ` boundary-orbit map for one periodically ordered
# closed boundary component.
# For rotation order `n=sym.order`, the `l`-th group image acts on the boundary
# indices as q -> q + orientation*l*N/n    (mod N).
# orientation=1 corresponds to the canonical outer-boundary orientation, while
# orientation=-1 is used for hole boundaries whose ordering has been reversed.
# The irrep character of the `l`-th rotational image is χ_l=exp(2πim*sym.sector*l/n).
function _rotation_orbits(::Type{T},N::Int,sym::BilliardGeometry.NFoldRotation,orientation::Int) where {T<:Real}
    n=sym.order
    n>=2||throw(ArgumentError("Rotation order must be at least two; received n=$n"))
    N%n==0||throw(ArgumentError("C$n symmetry requires N divisible by $n; received N=$N"))
    perms=Vector{Vector{Int}}(undef,n)
    scales=Vector{Complex{T}}(undef,n)
    @inbounds for l in 0:n-1
        p=Vector{Int}(undef,N)
        for q in 1:N
            p[q]=_idx_rotate(q,N,n,orientation*l)
        end
        perms[l+1]=p
        scales[l+1]=cis(T(2*pi)*T(sym.sector*l)/T(n))
    end
    return _build_periodic_symmetry_orbit_map(T,perms,scales)
end

function _reflection_orbits(::Type{T},N::Int,χ::Complex{T},idxmap) where {T<:Real}
    N%4==0||throw(ArgumentError("Reflection symmetry requires N divisible by four; received N=$N"))
    id=collect(1:N)
    refl=Vector{Int}(undef,N)
    @inbounds for q in 1:N
        refl[q]=idxmap(q,N)
    end
    return _build_periodic_symmetry_orbit_map(T,[id,refl],[one(Complex{T}),χ])
end

"""
    _composite_reflection_orbits(::Type{T},N::Int,sym::BilliardGeometry.CompositeReflection,reversed::Bool=false) where {T<:Real} → SymmetryOrbitMap{T}

Build the complete boundary-orbit map generated by a
[`CompositeReflection`](@ref).

The constituent reflection constraints are first expanded into primitive index
generators. `XAxisReflection` and `YAxisReflection` each contribute one
generator, while `XYAxisReflection` contributes both coordinate-axis
reflections. Diagonal and anti-diagonal reflections contribute their
corresponding periodic index maps.

Starting from the identity, the complete finite symmetry group is generated by
closure under composition. The irrep factor of a generated action is the product
of the factors of its generators. If the same geometric action is generated
with two different factors, the requested reflection parities are inconsistent
and an `ArgumentError` is thrown.

For a canonically ordered outer boundary, `reversed=false`. Interior-hole
components are reversed during CFIE boundary construction, so `reversed=true`
interchanges the periodic index formulas associated with the diagonal and
anti-diagonal reflections. The physical reflections themselves are unchanged.

`XYAxisReflection` contributes the x- and y-reflection generators rather than
an explicit `π` rotation. Their product, together with all other required group
elements, is generated automatically by the closure.

## Arguments
* `T`: Real scalar type used for complex irrep factors.
* `N`: Number of nodes on the closed periodic boundary component.
* `sym`: Composite collection of reflection constraints.
* `reversed`: Whether the component node ordering has been reversed.

## Returns
* `orbits::SymmetryOrbitMap{T}`: Exact orbit map for the complete generated symmetry group.
"""
function _composite_reflection_orbits(::Type{T},N::Int,sym::BilliardGeometry.CompositeReflection,reversed::Bool=false) where {T<:Real}
    isempty(sym.reflections)&&throw(ArgumentError("CompositeReflection requires at least one reflection"))
    genperms=Vector{Vector{Int}}()
    genscales=Complex{T}[]
    for ref in sym.reflections
        if ref isa BilliardGeometry.XAxisReflection
            push!(genperms,[_idx_reflect_x(q,N) for q in 1:N])
            push!(genscales,Complex{T}(ref.parity_y))
        elseif ref isa BilliardGeometry.YAxisReflection
            push!(genperms,[_idx_reflect_y(q,N) for q in 1:N])
            push!(genscales,Complex{T}(ref.parity_x))
        elseif ref isa BilliardGeometry.XYAxisReflection
            push!(genperms,[_idx_reflect_x(q,N) for q in 1:N])
            push!(genscales,Complex{T}(ref.parity_y))
            push!(genperms,[_idx_reflect_y(q,N) for q in 1:N])
            push!(genscales,Complex{T}(ref.parity_x))
        elseif ref isa BilliardGeometry.DiagonalReflection
            f=reversed ? _idx_reflect_diag_minus : _idx_reflect_diag_plus
            push!(genperms,[f(q,N) for q in 1:N])
            push!(genscales,BilliardGeometry.symmetry_irrep_character(T,ref))
        elseif ref isa BilliardGeometry.AntiDiagonalReflection
            f=reversed ? _idx_reflect_diag_plus : _idx_reflect_diag_minus
            push!(genperms,[f(q,N) for q in 1:N])
            push!(genscales,BilliardGeometry.symmetry_irrep_character(T,ref))
        else
            throw(ArgumentError("Unsupported reflection type $(typeof(ref)) in CompositeReflection"))
        end
    end
    perms=Vector{Int}[collect(1:N)]
    scales=Complex{T}[one(Complex{T})]
    head=1
    while head<=length(perms)
        p=perms[head]
        χ=scales[head]
        for g in eachindex(genperms)
            gp=genperms[g]
            pref=[gp[p[q]] for q in 1:N]
            χnew=genscales[g]*χ
            j=findfirst(==(pref),perms)
            if isnothing(j)
                push!(perms,pref)
                push!(scales,χnew)
            elseif scales[j]!=χnew
                throw(ArgumentError("Composite reflection parities are inconsistent"))
            end
        end
        head+=1
    end
    return _build_periodic_symmetry_orbit_map(T,perms,scales)
end

################################################################################
######################## PERIODIC SINGLE-BOUNDARY MAPS #########################
################################################################################

# The formulas below assume the canonical periodic boundary convention: the
# component begins on the positive x-axis and follows the canonical orientation.
# All actions are exact integer permutations of midpoint-node indices.
# Reflection maps reverse the periodic ordering and differ only by their phase
# shift. Rotations preserve the ordering and act by exact cyclic shifts.
@inline _idx_reflect_x(q::Int,N::Int)=mod1(N-q+1,N)
@inline _idx_reflect_y(q::Int,N::Int)=mod1(N÷2-q+1,N)
@inline _idx_reflect_diag_plus(q::Int,N::Int)=mod1(N÷4-q+1,N)
@inline _idx_reflect_diag_minus(q::Int,N::Int)=mod1(3*N÷4-q+1,N)
@inline _idx_rotate_pi(q::Int,N::Int)=mod1(q+N÷2,N)
@inline _idx_rotate(q::Int,N::Int,n::Int,l::Int)=mod1(q+l*(N÷n),N)

"""
    symmetry_index_orbits(::Type{T},N::Int,sym::BilliardGeometry.XAxisReflection) where {T<:Real} → SymmetryOrbitMap{T}

Build the exact two-element boundary orbits generated by reflection across the
x-axis.

For the canonical periodic boundary ordering, reflection acts exactly as

    q -> N-q+1    (mod N).

The reflected member of each orbit carries the parity `sym.parity_y`.
Fundamental representatives are extracted directly from the integer permutation
and need not form a contiguous index block.

## Returns
* `orbits::SymmetryOrbitMap{T}`: Orbit map with orbit size `2` and reduced size `N/2`.
"""
function symmetry_index_orbits(::Type{T},N::Int,sym::BilliardGeometry.XAxisReflection) where {T<:Real}
    iseven(N)||throw(ArgumentError("XAxisReflection requires an even node count; received N=$N"))
    id=collect(1:N)
    refl=Vector{Int}(undef,N)
    @inbounds for q in 1:N
        refl[q]=_idx_reflect_x(q,N)
    end
    χ=BilliardGeometry.symmetry_irrep_character(T,sym)
    return _build_periodic_symmetry_orbit_map(T,[id,refl],[one(Complex{T}),χ])
end

"""
    symmetry_index_orbits(::Type{T},N::Int,sym::BilliardGeometry.YAxisReflection) where {T<:Real} → SymmetryOrbitMap{T}

Build the exact two-element boundary orbits generated by reflection across the
y-axis.

For the canonical periodic boundary ordering, reflection acts exactly as

    q -> N/2-q+1    (mod N).

The reflected member of each orbit carries the parity `sym.parity_x`.
Fundamental representatives are extracted directly from the exact integer
permutation.

## Returns
* `orbits::SymmetryOrbitMap{T}`: Orbit map with orbit size `2` and reduced size `N/2`.
"""
function symmetry_index_orbits(::Type{T},N::Int,sym::BilliardGeometry.YAxisReflection) where {T<:Real}
    iseven(N)||throw(ArgumentError("YAxisReflection requires an even node count; received N=$N"))
    id=collect(1:N)
    refl=Vector{Int}(undef,N)
    @inbounds for q in 1:N
        refl[q]=_idx_reflect_y(q,N)
    end
    χ=BilliardGeometry.symmetry_irrep_character(T,sym)
    return _build_periodic_symmetry_orbit_map(T,[id,refl],[one(Complex{T}),χ])
end

"""
    symmetry_index_orbits(::Type{T},N::Int,sym::BilliardGeometry.XYAxisReflection) where {T<:Real} → SymmetryOrbitMap{T}

Build the complete four-element `D₂` orbit map generated by reflections across
the coordinate axes.

For each boundary node the four group images are

    I,
    R_x,
    R_y,
    R_xR_y,

with irrep factors

    1,
    parity_x,
    parity_y,
    parity_x*parity_y.

The node count must be divisible by four. Fundamental representatives are
selected directly from the exact integer permutations.

## Returns
* `orbits::SymmetryOrbitMap{T}`: Orbit map with orbit size `4` and reduced size `N/4`.
"""
function symmetry_index_orbits(::Type{T},N::Int,sym::BilliardGeometry.XYAxisReflection) where {T<:Real}
    N%4==0||throw(ArgumentError("XYAxisReflection requires N divisible by four; received N=$N"))
    id=collect(1:N)
    rx=Vector{Int}(undef,N)
    ry=Vector{Int}(undef,N)
    rxy=Vector{Int}(undef,N)
    @inbounds for q in 1:N
        rx[q]=_idx_reflect_x(q,N)
        ry[q]=_idx_reflect_y(q,N)
        rxy[q]=_idx_rotate_pi(q,N)
    end
    χx=Complex{T}(sym.parity_x)
    χy=Complex{T}(sym.parity_y)
    return _build_periodic_symmetry_orbit_map(T,[id,rx,ry,rxy],Complex{T}[1,χx,χy,χx*χy])
end

"""
    symmetry_index_orbits(::Type{T},N::Int,sym::DiagonalReflection) where {T<:Real} → SymmetryOrbitMap{T}

Build the exact two-element orbit map generated by reflection across the
diagonal `y=x`.

For the canonical periodic boundary ordering the index action is

    q -> N/4-q+1    (mod N).

The node count must be divisible by four. The reflected member of each orbit
carries the irrep parity stored by `sym`.

## Returns
* `orbits::SymmetryOrbitMap{T}`: Orbit map with orbit size `2` and reduced size `N/2`.
"""
function symmetry_index_orbits(::Type{T},N::Int,sym::DiagonalReflection) where {T<:Real}
    χ=symmetry_irrep_character(T,sym)
    return _reflection_orbits(T,N,χ,_idx_reflect_diag_plus)
end

"""
    symmetry_index_orbits(::Type{T},N::Int,sym::AntiDiagonalReflection) where {T<:Real} → SymmetryOrbitMap{T}

Build the exact two-element orbit map generated by reflection across the
anti-diagonal `y=-x`.

For the canonical periodic boundary ordering the index action is

    q -> 3N/4-q+1    (mod N).

The node count must be divisible by four. The reflected member of each orbit
carries the irrep parity stored by `sym`.

## Returns
* `orbits::SymmetryOrbitMap{T}`: Orbit map with orbit size `2` and reduced size `N/2`.
"""
function symmetry_index_orbits(::Type{T},N::Int,sym::AntiDiagonalReflection) where {T<:Real}
    χ=symmetry_irrep_character(T,sym)
    return _reflection_orbits(T,N,χ,_idx_reflect_diag_minus)
end

"""
    symmetry_index_orbits(::Type{T},N::Int,sym::BilliardGeometry.CompositeReflection) where {T<:Real} → SymmetryOrbitMap{T}

Build the complete symmetry-orbit map generated by a composite collection of
reflection constraints on one canonically ordered periodic boundary component.

For example, combining `XYAxisReflection` with `DiagonalReflection` on a
square-symmetric boundary generates the full eight-element `D₄` action.

The supplied parity factors must define a consistent one-dimensional
representation of the generated group.

## Returns
* `orbits::SymmetryOrbitMap{T}`: Exact orbit map for the generated reflection group.
"""
function symmetry_index_orbits(::Type{T},N::Int,sym::BilliardGeometry.CompositeReflection) where {T<:Real}
    return _composite_reflection_orbits(T,N,sym)
end

"""
    symmetry_index_orbits(::Type{T},pts::BoundaryPoints,sym::BilliardGeometry.CompositeReflection) where {T<:Real} → SymmetryOrbitMap{T}

Build the composite-reflection orbit map for one complete periodic [`BoundaryPoints`](@ref) component.
"""
@inline symmetry_index_orbits(::Type{T},pts::BoundaryPoints,sym::BilliardGeometry.CompositeReflection) where {T<:Real}=symmetry_index_orbits(T,length(pts),sym)

"""
    symmetry_index_orbits(::Type{T},N::Int,sym::BilliardGeometry.NFoldRotation) where {T<:Real} → SymmetryOrbitMap{T}

Build the complete `Cₙ` boundary-orbit decomposition associated with an
[`NFoldRotation`](@ref).

Although `sym` may describe one nontrivial rotational image, reduction uses the
complete cyclic group

    I,R,R²,...,Rⁿ⁻¹.

For irrep sector `s=sym.sector`, image `R^l` carries

    χ_l=exp(2πim*s*l/n).

The node count must be divisible by `n=sym.order`. Fundamental representatives
are extracted directly from the exact cyclic index permutations.

## Returns
* `orbits::SymmetryOrbitMap{T}`: Orbit map with orbit size `n` and reduced size `N/n`.
"""
function symmetry_index_orbits(::Type{T},N::Int,sym::BilliardGeometry.NFoldRotation) where {T<:Real}
    return _rotation_orbits(T,N,sym,1)
end

"""
    symmetry_index_orbits(::Type{T},N,syms::AbstractVector{<:BilliardGeometry.NFoldRotation})
Build a complete cyclic orbit map from the nontrivial rotation images returned by
[`Cn_symmetry`](@ref).
All supplied images must belong to the same cyclic group and irrep sector. Since
the complete orbit is determined entirely by the common order and sector, one
representative image is sufficient to construct the map.
"""
function symmetry_index_orbits(::Type{T},N::Int,syms::AbstractVector{<:BilliardGeometry.NFoldRotation}) where {T<:Real}
    isempty(syms)&&throw(ArgumentError("Rotation image collection cannot be empty"))
    order=syms[1].order
    sector=syms[1].sector
    all(s->s.order==order&&s.sector==sector,syms)||throw(ArgumentError("All NFoldRotation images must have identical order and irrep sector"))
    return symmetry_index_orbits(T,N,syms[1])
end

################################################################################
######################## MULTICOMPONENT EXACT ORBITS ############################
################################################################################
# A multicomponent boundary consists of one outer physical boundary followed by
# zero or more interior holes:
#
#     pts = [outer,hole₁,hole₂,...].
#
# Each BoundaryPoints object represents one complete closed physical component,
# even when that component was geometrically constructed from several curve
# pieces. Under the supported symmetries every physical component is individually
# invariant; symmetry does not permute distinct boundary components.
#
# The full multicomponent orbit map is therefore obtained by constructing the
# ordinary periodic symmetry map independently on every component and then
# concatenating the resulting maps in flattened boundary order.
#
# Different components may have different node counts. Only each individual
# component must satisfy the divisibility requirements of the active symmetry.

"""
    _combine_component_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},maps::Vector{SymmetryOrbitMap{T}}) where {T<:Real} → SymmetryOrbitMap{T}

Combine independently symmetry-reduced closed boundary components into one
global [`SymmetryOrbitMap`](@ref).

Each entry of `pts` is one complete connected physical boundary component. The
corresponding entry of `maps` is its local periodic symmetry-orbit map. Distinct
components may have different node counts and reduced dimensions, but must have
the same symmetry-group order.

Local full-boundary indices are shifted by the component offsets and local
fundamental indices are concatenated in component order.

## Arguments
* `T`: Real scalar type used by the complex irrep factors.
* `pts::Vector{BoundaryPoints{T}}`: Complete connected boundary components.
* `maps::Vector{SymmetryOrbitMap{T}}`: Local symmetry-orbit maps, one per component.

## Returns
* `orbits::SymmetryOrbitMap{T}`: Global orbit map for the flattened multicomponent boundary.
"""
function _combine_component_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},maps::Vector{SymmetryOrbitMap{T}}) where {T<:Real}
    nc=length(pts)
    length(maps)==nc||throw(DimensionMismatch("Received $(length(maps)) symmetry maps for $nc boundary components"))
    nc>0||throw(ArgumentError("Boundary cannot be empty"))
    ng=orbit_size(maps[1])
    all(m->orbit_size(m)==ng,maps)||throw(DimensionMismatch("All boundary components must use the same symmetry-group order"))
    offs=component_offsets(pts)
    Ntot=offs[end]-1
    Nred=sum(fundamental_size,maps)
    Ifund=Vector{Int}(undef,Nred)
    full_to_fund=Vector{Int}(undef,Ntot)
    full_to_scale=Vector{Complex{T}}(undef,Ntot)
    fund_to_full=Matrix{Int}(undef,ng,Nred)
    fund_to_scale=Matrix{Complex{T}}(undef,ng,Nred)
    bred=0
    @inbounds for a in 1:nc
        map=maps[a]
        off=offs[a]-1
        ma=fundamental_size(map)
        for b in 1:ma
            bg=bred+b
            Ifund[bg]=off+map.Ifund[b]
            for g in 1:ng
                q=off+map.fund_to_full[g,b]
                χ=map.fund_to_scale[g,b]
                fund_to_full[g,bg]=q
                fund_to_scale[g,bg]=χ
                full_to_fund[q]=bg
                full_to_scale[q]=χ
            end
        end
        bred+=ma
    end
    return SymmetryOrbitMap{T}(Ifund,full_to_fund,full_to_scale,fund_to_full,fund_to_scale)
end

"""
    symmetry_index_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},sym::BilliardGeometry.XAxisReflection) where {T<:Real} → SymmetryOrbitMap{T}

Build the global x-axis-reflection orbit map for one outer boundary and any
number of interior holes.

Each connected boundary component is reflected onto itself and reduced using
its own periodic node ordering. Components may have different node counts.
"""
function symmetry_index_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},sym::BilliardGeometry.XAxisReflection) where {T<:Real}
    maps=[symmetry_index_orbits(T,p,sym) for p in pts]
    return _combine_component_orbits(T,pts,maps)
end

"""
    symmetry_index_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},sym::BilliardGeometry.YAxisReflection) where {T<:Real} → SymmetryOrbitMap{T}

Build the global y-axis-reflection orbit map for one outer boundary and any
number of interior holes.

Each connected boundary component is reflected onto itself and reduced using
its own periodic node ordering. Components may have different node counts.
"""
function symmetry_index_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},sym::BilliardGeometry.YAxisReflection) where {T<:Real}
    maps=[symmetry_index_orbits(T,p,sym) for p in pts]
    return _combine_component_orbits(T,pts,maps)
end

"""
    symmetry_index_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},sym::BilliardGeometry.XYAxisReflection) where {T<:Real} → SymmetryOrbitMap{T}

Build the global `D₂` orbit map for one outer boundary and any number of
interior holes.

Every connected component is individually invariant under both coordinate
reflections. Components may have different node counts, provided each local
count is compatible with the four-element symmetry group.
"""
function symmetry_index_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},sym::BilliardGeometry.XYAxisReflection) where {T<:Real}
    maps=[symmetry_index_orbits(T,p,sym) for p in pts]
    return _combine_component_orbits(T,pts,maps)
end

"""
    symmetry_index_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},sym::DiagonalReflection) where {T<:Real} → SymmetryOrbitMap{T}

Build the global `y=x` diagonal-reflection orbit map for one outer boundary and
any number of interior holes.

Every connected physical boundary component is assumed to be individually
invariant under the diagonal reflection and is reduced independently using its
exact periodic index action. Distinct components may have different node counts.

The resulting local orbit maps are concatenated in flattened boundary order.

## Returns
* `orbits::SymmetryOrbitMap{T}`: Global two-element reflection orbit map.
"""
function symmetry_index_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},sym::DiagonalReflection) where {T<:Real}
    isempty(pts)&&throw(ArgumentError("Boundary cannot be empty"))
    χ=symmetry_irrep_character(T,sym)
    maps=Vector{SymmetryOrbitMap{T}}(undef,length(pts))
    maps[1]=_reflection_orbits(T,length(pts[1]),χ,_idx_reflect_diag_plus)
    @inbounds for a in 2:length(pts)
        maps[a]=_reflection_orbits(T,length(pts[a]),χ,_idx_reflect_diag_minus)
    end
    return _combine_component_orbits(T,pts,maps)
end

"""
    symmetry_index_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},sym::AntiDiagonalReflection) where {T<:Real} → SymmetryOrbitMap{T}

Build the global `y=-x` anti-diagonal-reflection orbit map for one outer boundary
and any number of interior holes.

Every connected physical boundary component is assumed to be individually
invariant under the anti-diagonal reflection and is reduced independently using
its exact periodic index action. Distinct components may have different node
counts.

The resulting local orbit maps are concatenated in flattened boundary order.

## Returns
* `orbits::SymmetryOrbitMap{T}`: Global two-element reflection orbit map.
"""
function symmetry_index_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},sym::AntiDiagonalReflection) where {T<:Real}
    isempty(pts)&&throw(ArgumentError("Boundary cannot be empty"))
    χ=symmetry_irrep_character(T,sym)
    maps=Vector{SymmetryOrbitMap{T}}(undef,length(pts))
    maps[1]=_reflection_orbits(T,length(pts[1]),χ,_idx_reflect_diag_minus)
    @inbounds for a in 2:length(pts)
        maps[a]=_reflection_orbits(T,length(pts[a]),χ,_idx_reflect_diag_plus)
    end
    return _combine_component_orbits(T,pts,maps)
end

"""
    symmetry_index_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},sym::BilliardGeometry.CompositeReflection) where {T<:Real} → SymmetryOrbitMap{T}

Build the global composite-reflection orbit map for one outer boundary and any
number of interior holes.

Every connected physical boundary component is assumed to be individually
invariant under the symmetry group generated by `sym`. A complete local group
closure is constructed independently for every component and the resulting
orbit maps are combined in flattened boundary order.

## Returns
* `orbits::SymmetryOrbitMap{T}`: Global orbit map for the complete generated reflection group.
"""
function symmetry_index_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},sym::BilliardGeometry.CompositeReflection) where {T<:Real}
    isempty(pts)&&throw(ArgumentError("Boundary cannot be empty"))
    maps=Vector{SymmetryOrbitMap{T}}(undef,length(pts))
    maps[1]=_composite_reflection_orbits(T,length(pts[1]),sym)
    @inbounds for a in 2:length(pts)
        maps[a]=_composite_reflection_orbits(T,length(pts[a]),sym,true)
    end
    return _combine_component_orbits(T,pts,maps)
end

"""
    symmetry_index_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},sym::BilliardGeometry.NFoldRotation) where {T<:Real} → SymmetryOrbitMap{T}

Build the global `Cₙ` orbit map for one outer boundary and any number of
interior holes.

Every connected component is individually invariant under the rotation group
and is reduced independently using its periodic node ordering. Components may
have different node counts, provided each count is divisible by the rotation
order.

## Returns
* `orbits::SymmetryOrbitMap{T}`: Global cyclic orbit map.
"""
function symmetry_index_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},sym::BilliardGeometry.NFoldRotation) where {T<:Real}
    isempty(pts)&&throw(ArgumentError("Boundary cannot be empty"))
    maps=Vector{SymmetryOrbitMap{T}}(undef,length(pts))
    maps[1]=_rotation_orbits(T,length(pts[1]),sym,1)
    @inbounds for a in 2:length(pts)
        maps[a]=_rotation_orbits(T,length(pts[a]),sym,-1)
    end
    return _combine_component_orbits(T,pts,maps)
end

"""
    symmetry_index_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},syms::AbstractVector{<:BilliardGeometry.NFoldRotation}) where {T<:Real} → SymmetryOrbitMap{T}

Build the global cyclic orbit map for [`Cn_symmetry`](@ref).

All supplied rotations must belong to the same cyclic group and irrep sector.
Each connected physical boundary component is reduced independently.
"""
function symmetry_index_orbits(::Type{T},pts::Vector{BoundaryPoints{T}},syms::AbstractVector{<:BilliardGeometry.NFoldRotation}) where {T<:Real}
    isempty(syms)&&throw(ArgumentError("Rotation image collection cannot be empty"))
    order=syms[1].order
    sector=syms[1].sector
    all(s->s.order==order&&s.sector==sector,syms)||throw(ArgumentError("All NFoldRotation images must have identical order and irrep sector"))
    return symmetry_index_orbits(T,pts,syms[1])
end

@inline symmetry_index_orbits(::Type{T},pts::BoundaryPoints,sym::BilliardGeometry.XAxisReflection) where {T<:Real}=symmetry_index_orbits(T,length(pts),sym)
@inline symmetry_index_orbits(::Type{T},pts::BoundaryPoints,sym::BilliardGeometry.YAxisReflection) where {T<:Real}=symmetry_index_orbits(T,length(pts),sym)
@inline symmetry_index_orbits(::Type{T},pts::BoundaryPoints,sym::BilliardGeometry.XYAxisReflection) where {T<:Real}=symmetry_index_orbits(T,length(pts),sym)
@inline symmetry_index_orbits(::Type{T},pts::BoundaryPoints,sym::DiagonalReflection) where {T<:Real}=symmetry_index_orbits(T,length(pts),sym)
@inline symmetry_index_orbits(::Type{T},pts::BoundaryPoints,sym::AntiDiagonalReflection) where {T<:Real}=symmetry_index_orbits(T,length(pts),sym)
@inline symmetry_index_orbits(::Type{T},pts::BoundaryPoints,sym::BilliardGeometry.NFoldRotation) where {T<:Real}=symmetry_index_orbits(T,length(pts),sym)
@inline symmetry_index_orbits(::Type{T},pts::BoundaryPoints,syms::AbstractVector{<:BilliardGeometry.NFoldRotation}) where {T<:Real}=symmetry_index_orbits(T,length(pts),syms)