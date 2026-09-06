"""
    BoundaryPoints{T} <: Any

`BoundaryPoints` is a concrete type that collects the boundary discretization data
(sampled points, normals, arc-length coordinates and quadrature weights) used by
the solvers to construct the boundary matrices.

## Description
Instances are produced either by [`boundary_coords`](@ref), which samples the full
boundary (including curves marked with `QuantumSolverIgnore`) to compute `xy`,
`normal`, `s` and `ds`, or by a solver's `evaluate_points` method, which samples
only the physical boundary and computes solver-specific quadrature weights such as
`w_vs` (see [`VerginiSaracenoSolver`](@ref)) or `w_dm` (see
[`DecompositionMethodSolver`](@ref)). Fields that are not populated by a given
constructor call default to empty vectors. The inner constructor validates that
all non-empty vector fields have the same length as `xy`.

## Attributes
* `xy`: Boundary points in Cartesian coordinates.
* `normal`: Outward unit normal vectors at each boundary point.
* `kappa`: Curvature of the boundary at each point (reserved, currently unused).
* `s`: Arc-length coordinate of each boundary point, measured continuously along the full composite boundary.
* `ds`: Arc-length quadrature element at each boundary point.
* `rdotn`: Dot product of the position vector with the normal, `r ⋅ n` (reserved, currently unused).
* `w_vs`: Quadrature weights for the Vergini–Saraceno method, see [`VerginiSaracenoSolver`](@ref).
* `w_dm`: Quadrature weights for the decomposition method, see [`DecompositionMethodSolver`](@ref).
* `xy_int`: Interior points (reserved, currently unused).

## API
The following functions can be evaluated for this type:
- [`boundary_coords`](@ref)
- `Base.length`
- `Base.isempty`
"""
struct BoundaryPoints{T<:Real}
    xy::Vector{SVector{2,T}}
    normal::Vector{SVector{2,T}}
    kappa::Vector{T}
    s::Vector{T}
    ds::Vector{T}
    rdotn::Vector{T}
    w_vs::Vector{T}
    w_dm::Vector{T}
    xy_int::Vector{SVector{2,T}}
    
    # Inner constructor with validation
    function BoundaryPoints{T}(xy, normal, kappa, s, ds, rdotn, w_vs, w_dm, xy_int) where T<:Real
        n = length(xy)
        # Validate that non-empty vectors have consistent lengths
        for (name, vec) in [(:normal, normal), (:ds, ds), (:rdotn, rdotn), 
                             (:w_vs, w_vs), (:w_dm, w_dm)]
            if !isempty(vec) && length(vec) != n
                error("Length of $name ($(length(vec))) must match xy ($n)")
            end
        end
        new{T}(xy, normal, kappa, s, ds, rdotn, w_vs, w_dm, xy_int)
    end
end

# 2. Convenience constructor to infer T from xy
"""
    BoundaryPoints(xy::Vector{SVector{2,T}}; normal = SVector{2,T}[], kappa = T[], s = T[], ds = T[], rdotn = T[], w_vs = T[], w_dm = T[], xy_int = SVector{2,T}[]) where T<:Real → bp::BoundaryPoints{T}

Constructs a [`BoundaryPoints`](@ref) instance from the boundary points `xy`,
inferring the element type `T` from `xy` and defaulting all other fields to empty
vectors when not supplied as keyword arguments.

## Arguments
* `xy`: Vector of boundary points in Cartesian coordinates.

## Keyword arguments
* `normal::Vector{SVector{2,T}} = SVector{2,T}[]`: Outward unit normal vectors at each boundary point.
* `kappa::Vector{T} = T[]`: Curvature of the boundary at each point.
* `s::Vector{T} = T[]`: Arc-length coordinate of each boundary point.
* `ds::Vector{T} = T[]`: Arc-length quadrature element at each boundary point.
* `rdotn::Vector{T} = T[]`: Dot product of the position vector with the normal, `r ⋅ n`.
* `w_vs::Vector{T} = T[]`: Quadrature weights for the Vergini–Saraceno method.
* `w_dm::Vector{T} = T[]`: Quadrature weights for the decomposition method.
* `xy_int::Vector{SVector{2,T}} = SVector{2,T}[]`: Interior points.

## Returns
* `bp`: A [`BoundaryPoints{T}`](@ref) instance holding the supplied boundary data.
"""
function BoundaryPoints(xy::Vector{SVector{2,T}}; 
                        normal=SVector{2,T}[], 
                        kappa=T[],
                        s=T[], 
                        ds=T[], 
                        rdotn=T[], 
                        w_vs=T[], 
                        w_dm=T[], 
                        xy_int=SVector{2,T}[]) where T<:Real
    return BoundaryPoints{T}(xy, normal, kappa, s, ds, rdotn, w_vs, w_dm, xy_int)
end

# 3. Add useful methods
"""
    length(bp::BoundaryPoints) → n::Int

Returns the number of sampled boundary points, `n = length(bp.xy)`.
"""
Base.length(bp::BoundaryPoints) = length(bp.xy)

"""
    isempty(bp::BoundaryPoints) → flag::Bool

Returns `true` if `bp` contains no boundary points, i.e. `isempty(bp.xy)`.
"""
Base.isempty(bp::BoundaryPoints) = isempty(bp.xy)

function _determine_bp_sizes(curves, bs, k)
    Ns = Vector{Int64}(undef,length(curves)) # store the data to indexwise access. This needs to be this way b/c we dont know beforehand which curves are real and which are abstract. Use sizehint! to give an idea as to not need to resize b/c it could the that real and abstract curves and intermingled
    @inbounds for i in eachindex(curves) # make an initial size calculation of the resulting vectors
        crv=curves[i]
        Ns[i] =max(20,round(Int,k*crv.length*bs[i]/2*pi))
    end
    return Ns
end


"""
    boundary_coords(billiard::Bi, samplers::Vector{AbsSampler}, Ns::Vector{Int64}) where {Bi<:AbsBilliard} → bp::BoundaryPoints

Samples the full boundary of `billiard` (including curves marked with
`QuantumSolverIgnore`) and returns the sampled points, outward normals and
arc-length coordinates as a [`BoundaryPoints`](@ref) instance.

## Description
Each boundary curve, as returned by [`get_boundary_curves_with_ignored`](@ref), is
sampled independently with its own sampler and number of points using
`sample_points`. The per-curve arc-length coordinates are offset by the cumulative
length of the preceding curves, so that `s` runs continuously over the whole
composite boundary.

## Arguments
* `billiard`: The billiard whose boundary is sampled.
* `samplers`: Vector of samplers, one per boundary curve (including ignored curves).
* `Ns`: Vector with the number of sample points to generate for each boundary curve.

## Returns
* `bp`: A [`BoundaryPoints`](@ref) instance with the `xy`, `normal`, `s` and `ds` fields populated.
"""
function boundary_coords(billiard::Bi, samplers::Vector{AbsSampler}, Ns::Vector{Int64}) where {Bi<:AbsBilliard}
    curves = get_boundary_curves_with_ignored(billiard)
    T = typeof(curves[1].length)
    M = length(Ns)
    xy_all = Vector{Vector{SVector{2,T}}}(undef, M)
    normal_all = Vector{Vector{SVector{2,T}}}(undef, M)
    s_all = Vector{Vector{T}}(undef, M)
    ds_all = Vector{Vector{T}}(undef, M)
    L0 = zero(T)
    for i in eachindex(curves)
        crv = curves[i]
        L = crv.length
        sampler = samplers[i]
        t, dt = sample_points(sampler, Ns[i])
        ds = L*dt #this needs modification!!!
        xy = curve(crv,t)
        normal = domain_gradient_vector(crv, xy)
        normal .= normal./norm(normal)
        #rn = dot.(xy, normal)
        xy_all[i] = xy
        normal_all[i] = normal
        s_all[i] = arc_length(crv,t) .+ L0 #arc_lengt(crv, xy)
        ds_all[i] = ds  
        #w_n_all[i] = (ds.*rn)./(2.0*k.^2)
        L0 += L
    end

    return BoundaryPoints(vcat(xy_all...); normal = vcat(normal_all...), s=vcat(s_all...), ds = vcat(ds_all...) )
end

"""
    boundary_coords(billiard::Bi, fourier_sampler::FourierNodes, M) where {Bi<:AbsBilliard} → bp::BoundaryPoints

Samples the full boundary of `billiard` (including curves marked with
`QuantumSolverIgnore`) using a single `FourierNodes` sampler that
distributes a total of `M` points across all boundary curves. Otherwise behaves
like [`boundary_coords(::AbsBilliard, ::Vector{AbsSampler}, ::Vector{Int64})`](@ref).

## Arguments
* `billiard`: The billiard whose boundary is sampled.
* `fourier_sampler`: The `FourierNodes` sampler used to distribute the `M` points among all boundary curves.
* `M`: Total number of sample points to generate over the whole boundary.

## Returns
* `bp`: A [`BoundaryPoints`](@ref) instance with the `xy`, `normal`, `s` and `ds` fields populated.
"""
function boundary_coords(billiard::Bi, fourier_sampler::FourierNodes, M) where {Bi<:AbsBilliard}
    curves = get_boundary_curves_with_ignored(billiard)
    T = typeof(curves[1].length)
    n_curves = length(curves)

    ts,dts = sample_points(fourier_sampler, M)
    xy_all = Vector{Vector{SVector{2,T}}}(undef, n_curves)
    normal_all = Vector{Vector{SVector{2,T}}}(undef, n_curves)
    s_all = Vector{Vector{T}}(undef, n_curves)
    ds_all = Vector{Vector{T}}(undef, n_curves)
    #w_n_all = Vector{Vector{T}}(undef, M)
    L0 = zero(T)
    for i in eachindex(curves)
        crv = curves[i]
        L = crv.length
        t = ts[i]
        dt = dts[i]
        ds = L*dt #this needs modification!!!
        xy = curve(crv,t)
        normal = domain_gradient_vector(crv, xy)
        normal .= normal./norm(normal)
        #rn = dot.(xy, normal)
        xy_all[i] = xy
        normal_all[i] = normal
        s_all[i] = arc_length(crv,t) .+ L0 #arc_lengt(crv, xy)
        ds_all[i] = ds  
        #w_n_all[i] = (ds.*rn)./(2.0*k.^2)
        L0 += L
    end

    return BoundaryPoints(vcat(xy_all...); normal = vcat(normal_all...), s=vcat(s_all...), ds = vcat(ds_all...))
end

"""
    get_boundary_curves_with_ignored(domain::D) where D<:AbsSimpleDomain → boundary::Vector{AbsCurve}

Returns the connected boundary curves of `domain` used for full-boundary sampling,
including both `SpecularReflection` curves and curves marked with
`QuantumSolverIgnore`.

## Description
This differs from `get_boundary_curves`, which only retains `SpecularReflection`
curves used for constructing the solver matrices: `QuantumSolverIgnore` curves are
excluded there but are needed here so that quantities such as arc length and
boundary points, via [`boundary_coords`](@ref), can be computed over the entire
physical boundary of the domain.

## Arguments
* `domain`: A simple domain whose boundary curves are collected.

## Returns
* `boundary`: The connected vector of boundary curves.
"""
function get_boundary_curves_with_ignored(domain::D) where D<:AbsSimpleDomain
    is_outer(crv) = (typeof(crv.bc) <: SpecularReflection || typeof(crv.bc) <: QuantumSolverIgnore)
    boundary = filter(is_outer, domain.boundary)
    return connect_curves(boundary)
end


"""
    get_boundary_curves_with_ignored(composite_domain::D) where D<:AbsCompositeDomain → boundary::Vector{AbsCurve}

Returns the connected boundary curves gathered over all subdomains of
`composite_domain`, including both `SpecularReflection` and `QuantumSolverIgnore`
curves. See [`get_boundary_curves_with_ignored`](@ref) for details.

## Arguments
* `composite_domain`: A composite domain whose subdomains' boundary curves are collected.

## Returns
* `boundary`: The connected vector of boundary curves.
"""
function get_boundary_curves_with_ignored(composite_domain::D) where D<:AbsCompositeDomain
    boundary = Vector{AbsCurve}()
    for domain in composite_domain.subdomains
        subboundary = get_boundary_curves(domain)
        append!(boundary,subboundary)
    end
    return connect_curves(boundary)
end

"""
    get_boundary_curves_with_ignored(billiard::B) where B<:AbsBilliard → boundary::Vector{AbsCurve}

Returns the connected boundary curves of the fundamental domain of `billiard`,
including both `SpecularReflection` and `QuantumSolverIgnore` curves. See
[`get_boundary_curves_with_ignored`](@ref) for details.

## Arguments
* `billiard`: The billiard whose fundamental domain's boundary curves are collected.

## Returns
* `boundary`: The connected vector of boundary curves.
"""
function get_boundary_curves_with_ignored(billiard::B) where B<:AbsBilliard
    return get_boundary_curves_with_ignored(billiard.fundamental_domain)
end

"""
    random_interior_points(billiard::Bi, N::Int; grd::Int = 1000) where {Bi<:AbsBilliard} → pts::Vector{SVector{2,T}}

Generates `N` points sampled uniformly (by rejection) from the interior of
`billiard`, for use e.g. by [`ParticularSolutionsMethod`](@ref)'s
[`evaluate_points`](@ref).

## Description
Candidate points are drawn uniformly from the padded bounding box of
`billiard`'s boundary curves (see [`boundary_limits`](@ref)) and accepted only
if [`is_inside`](@ref) confirms they lie within the billiard, until `N`
interior points have been collected.

## Arguments
* `billiard`: The billiard whose interior is sampled.
* `N`: The number of interior points to generate.

## Keyword arguments
* `grd::Int = 1000`: Target sampling density used to determine the bounding box, see [`boundary_limits`](@ref).

## Returns
* `pts`: A `Vector{SVector{2,T}}` of `N` interior points.
"""
function random_interior_points(billiard::Bi, N::Int; grd::Int=1000) where {Bi<:AbsBilliard}
    xlim, ylim = boundary_limits(get_boundary_curves(billiard); grd=grd)
    dx = xlim[2] - xlim[1]
    dy = ylim[2] - ylim[1]
    T = typeof(dx)
    pts = Vector{SVector{2,T}}(undef, N)
    n = 0
    while n < N
        x = dx*rand() + xlim[1]
        y = dy*rand() + ylim[1]
        pt = SVector(x, y)
        if is_inside(billiard, pt)
            n += 1
            pts[n] = pt
        end
    end
    return pts
end