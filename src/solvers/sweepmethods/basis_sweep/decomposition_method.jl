"""
    DecompositionMethodSolver{T} <: SweepSolver

`DecompositionMethodSolver` is a concrete [`SweepSolver`](@ref) implementing the
boundary decomposition method for computing quantum billiard spectra by sweeping
over individual wavenumbers.

## Description
For a fixed wavenumber `k`, the method constructs the matrices `F` and `G` (see
[`construct_matrices`](@ref)) from a boundary quadrature with weights `ds` and
`w_n` (see [`evaluate_points`](@ref)) and extracts a tension `t = 1/λ0` from the
largest generalized eigenvalue λ0 of `F * x = λ * G * x` (see [`solve`](@ref)). A
sequence of tensions over a range of wavenumbers is minimized/scanned by
[`solve_wavenumber`](@ref) or [`k_sweep`](@ref) to locate the eigenvalues of the
billiard.

## Attributes
* `dim_scaling_factor`: Scaling factor used to determine the basis dimension from the boundary length and wavenumber.
* `pts_scaling_factor`: Vector of scaling factors, one per fundamental boundary curve, used to determine the number of boundary sampling points.
* `sampler`: Vector of samplers, one per fundamental boundary curve, used to generate boundary points.
* `eps`: Relative tolerance used to filter small eigenvalues in the generalized eigenvalue decomposition.
* `min_dim`: Minimum basis dimension.
* `min_pts`: Minimum number of boundary sampling points.
* `rellich_origin`: Origin `c₀` used in the Rellich boundary weight.

## API
The following functions can be evaluated for this type:
- [`evaluate_points`](@ref)
- [`construct_matrices`](@ref)
- [`solve`](@ref)
- [`solve_vect`](@ref)
- [`solve_wavenumber`](@ref)
- [`k_sweep`](@ref)
"""
struct DecompositionMethodSolver{T} <: SweepSolver where {T<:Real}
    dim_scaling_factor::T
    pts_scaling_factor::Vector{T}
    sampler::Vector
    eps::T
    min_dim::Int64
    min_pts::Int64
    rellich_origin::SVector{2,T}
end

"""
    DecompositionMethodSolver(dim_scaling_factor::T,pts_scaling_factor::Union{T,Vector{T}};min_dim::Int=100,min_pts::Int=500) where {T<:Real} → solver::DecompositionMethodSolver{T}

Constructs a [`DecompositionMethodSolver`](@ref) with a single `GaussLegendreNodes`
sampler shared by every fundamental boundary curve.

## Arguments
* `dim_scaling_factor`: Scaling factor used to determine the basis dimension.
* `pts_scaling_factor`: Scaling factor, or vector thereof (one per fundamental boundary curve), used to determine the number of boundary sampling points.

## Keyword arguments
* `min_dim::Int = 100`: Minimum basis dimension.
* `min_pts::Int = 500`: Minimum number of boundary sampling points.
* `rellich_origin`: Origin `c₀` used in the Rellich boundary weight.

## Returns
* `solver`: A [`DecompositionMethodSolver{T}`](@ref) instance.
"""
function DecompositionMethodSolver(dim_scaling_factor::T,pts_scaling_factor::Union{T,Vector{T}};min_dim::Int=100,min_pts::Int=500,rellich_origin::SVector{2,T}=SVector{2,T}(zero(T),zero(T))) where {T<:Real}
    d=dim_scaling_factor
    bs=pts_scaling_factor isa T ? [pts_scaling_factor] : pts_scaling_factor
    sampler=[BilliardGeometry.GaussLegendreNodes()]
    return DecompositionMethodSolver(d,bs,sampler,eps(T),min_dim,min_pts,rellich_origin)
end

"""
    DecompositionMethodSolver(dim_scaling_factor::T,pts_scaling_factor::Union{T,Vector{T}},samplers::Vector{<:AbsSampler};min_dim::Int=100,min_pts::Int=500) where {T<:Real} → solver::DecompositionMethodSolver{T}

Constructs a [`DecompositionMethodSolver`](@ref) with a user-supplied sampler for
each fundamental boundary curve.

## Arguments
* `dim_scaling_factor`: Scaling factor used to determine the basis dimension.
* `pts_scaling_factor`: Scaling factor, or vector thereof (one per fundamental boundary curve), used to determine the number of boundary sampling points.
* `samplers`: Vector of samplers, one per fundamental boundary curve.

## Keyword arguments
* `min_dim::Int = 100`: Minimum basis dimension.
* `min_pts::Int = 500`: Minimum number of boundary sampling points.
* `rellich_origin`: Origin `c₀` used in the Rellich boundary weight.

## Returns
* `solver`: A [`DecompositionMethodSolver{T}`](@ref) instance.
"""
function DecompositionMethodSolver(dim_scaling_factor::T,pts_scaling_factor::Union{T,Vector{T}},samplers::Vector{<:BilliardGeometry.AbsSampler};min_dim::Int=100,min_pts::Int=500,rellich_origin::SVector{2,T}=SVector{2,T}(zero(T),zero(T))) where {T<:Real}
    d=dim_scaling_factor
    bs=pts_scaling_factor isa T ? [pts_scaling_factor] : pts_scaling_factor
    return DecompositionMethodSolver(d,bs,samplers,eps(T),min_dim,min_pts,rellich_origin)
end

"""
    evaluate_points(solver::DecompositionMethodSolver,billiard::Bi,k) where {Bi<:BilliardGeometry.AbsBilliard} → pts::BoundaryPoints

Samples the boundary of `billiard` and computes the boundary decomposition method
quadrature weights needed to construct the matrices in [`construct_matrices`](@ref).

## Description
The scaling factors and samplers are first adjusted to match the number of
fundamental boundary curves with [`adjust_scaling_and_samplers`](@ref). Each curve
is then sampled with its own sampler using [`boundary_coords`](@ref), which
provides the boundary coordinates, outward normals, arc-length coordinates and
arc-length quadrature elements. The normal-derivative decomposition weight is

    w_n = ds * ((x-c₀) ⋅ n) / (2 k²),

where `c₀=solver.rellich_origin`, `x` is the boundary point and `n` its outward unit normal.

## Arguments
* `solver`: The [`DecompositionMethodSolver`](@ref) used to determine the sampling parameters.
* `billiard`: The billiard whose boundary is sampled.
* `k`: The wavenumber used to determine the number of boundary sampling points and the `w_n` weights.

## Returns
* `pts`: A [`BoundaryPoints`](@ref) instance with the `xy`, `normal`, `s`, `ds` and `w_n` fields populated.
"""
function evaluate_points(solver::DecompositionMethodSolver,billiard::Bi,k) where {Bi<:BilliardGeometry.AbsBilliard}
    bs,samplers=adjust_scaling_and_samplers(solver,billiard)
    curves=BilliardGeometry.get_boundary_curves(billiard)
    T=eltype(solver.pts_scaling_factor)
    Ns=_determine_bp_sizes(curves,bs,k)
    M=length(Ns)
    xy_all=Vector{Vector{SVector{2,T}}}(undef,M)
    normal_all=Vector{Vector{SVector{2,T}}}(undef,M)
    s_all=Vector{Vector{T}}(undef,M)
    ds_all=Vector{Vector{T}}(undef,M)
    w_n_all=Vector{Vector{T}}(undef,M)
    c0=solver.rellich_origin
    L0=zero(T)
    @inbounds for i in eachindex(curves)
        xy,normal,s,ds=boundary_coords(curves[i],samplers[i],Ns[i])
        rn=Vector{T}(undef,length(xy))
        @inbounds @simd for j in eachindex(xy)
            rn[j]=dot(xy[j]-c0,normal[j])
        end
        rnmin=minimum(rn)
        rnmin>zero(T)||throw(ArgumentError("DM Rellich origin $c0 is outside the star kernel: minimum (x-c0)⋅n=$rnmin"))
        xy_all[i]=xy
        normal_all[i]=normal
        s_all[i]=s.+L0
        ds_all[i]=ds
        w_n_all[i]=(ds.*rn)./(T(2)*T(k)^2)
        L0+=curves[i].length
    end
    return BoundaryPoints(vcat(xy_all...);normal=vcat(normal_all...),s=vcat(s_all...),ds=vcat(ds_all...),w_n=vcat(w_n_all...))
end

"""
    construct_matrices(solver::DecompositionMethodSolver,basis::Ba,pts::BoundaryPoints,k;multithreaded::Bool=true) where {Ba<:AbsBasis} → (F::Matrix,G::Matrix)

Constructs the boundary decomposition method matrices `F` and `G` used to compute
the generalized eigenvalue problem `F * x = λ * G * x` at wavenumber `k`.

## Description
`F = B' * W * B`, where `B` is the [`basis_and_gradient_matrices`](@ref) basis
matrix at `k` and `W` is the diagonal quadrature weight matrix built from
`pts.ds`. `G = Bn' * Wn * Bn`, where `Bn = nx * dB/dx + ny * dB/dy` is the normal
derivative of the basis (built from `pts.normal`) and `Wn` is the diagonal
quadrature weight matrix built from `pts.w_n`. Both weight matrices are further
normalized by the number of basis symmetries `nsym`. Both `F` and `G` are
assembled with BLAS `syrk!` rank-k updates on the upper triangle, which is then
mirrored to the lower triangle, to minimize memory allocations.

## Arguments
* `solver`: The [`DecompositionMethodSolver`](@ref) whose matrices are constructed.
* `basis`: The basis used to evaluate `B` and its gradient.
* `pts`: The [`BoundaryPoints`](@ref) with sampled boundary points, normals and quadrature weights `ds`, `w_n`.
* `k`: The wavenumber at which the basis and its gradient are evaluated.

## Keyword arguments
* `multithreaded::Bool = true`: Whether the matrix construction is multithreaded.

## Returns
* `F`: The `F = B' * W * B` matrix.
* `G`: The `G = Bn' * Wn * Bn` matrix.
"""
function construct_matrices(solver::DecompositionMethodSolver,basis::Ba,pts::BoundaryPoints,k;multithreaded::Bool=true) where {Ba<:AbsBasis}
    @timeit_debug "construct_matrices" begin
        xy=pts.xy
        w=pts.ds
        wn=pts.w_n
        N=basis.dim
        M=length(xy)
        nsym=one(eltype(w))
        @debug "Matrix construction started" N M k
        @timeit_debug "basis_and_gradient_matrices" begin
            # the alogrithm consctructs B and the normal derivative Bn with syrk to minimize the allocation cost. It does this with the trick of putting sqrt(w_n) into both the rows of B and the rows of B' so that we can use syrk on sqrt(W)*B to get B'*(W*B) without forming W*B as a temporary matrix (posible b/c W is diagonal)
            @blas_1 B,dX,dY=basis_and_gradient_matrices(basis,k,xy;multithreaded)
            # Form F = B'*(W*B) by inplace scaling the rows of B by sqrt(w) (inplace to B) and use syrk to perform the Rank-k update of a symmetric matrix
            _scale_rows_sqrtw!(B,w,nsym) # trick of putting sqrt(w_n) into the rows of the transposed and original B to get (sqrt(W)*B)' * (sqrt(W)*B) so we can use syrk on sqrt(W)*B
        end
        @debug "Basis and gradient matrix computed" size=size(B)
        @timeit_debug "compute_F" begin
            F=Matrix{eltype(B)}(undef,N,N) # preallocate F
            @blas_multi_then_1 MAX_BLAS_THREADS BLAS.syrk!('U','T',one(eltype(B)),B,zero(eltype(B)),F) # F[u ∈ upper]+=1.0*B'*B, no need to fill(F,0) since the additive constant in C is 0
            _symmetrize_from_upper!(F) # since we chose "U" in syrk, we need to mirror upper -> lower
            # Build Bn into dX: dX <- nx*dX + ny*dY
            _build_Bn_inplace!(dX,dY,pts.normal)
            # Form G = Bn'*(Wn*Bn) by first scaling the rows of Bn (dX) by sqrt(w_n) (inplace to dX) and use syrk to perform the Rank-k update of a symmetric matrix
            _scale_rows_sqrtw!(dX,wn,nsym) # like for F form sqrt(Wn*Bn) with row scaling with the same trick of putting sqrt(w_n) into the rows of the transposed and original dX to get (sqrt(Wn)*Bn)' * (sqrt(Wn)*Bn) so we can use syrk on dX
        end
        @debug "F computed" size=size(F)
        @timeit_debug "compute_G" begin
            G=Matrix{eltype(B)}(undef,N,N) # preallocate G, no need to fill with zeros since we use zero(eltype(B)) for the additive constant in syrk
            @blas_multi_then_1 MAX_BLAS_THREADS BLAS.syrk!('U','T',one(eltype(B)),dX,zero(eltype(B)),G) # G[u ∈ upper]+=1.0*dX'*dX where dX is now sqrt(Wn)*Bn due to inplace scalings above
            _symmetrize_from_upper!(G) # since we chose "U" in syrk, we need to mirror upper -> lower
        end
        @debug "G computed" size=size(G)
        return F,G
    end
end

"""
    solve(solver::DecompositionMethodSolver,basis::Ba,pts::BoundaryPoints,k;multithreaded::Bool=true) where {Ba<:AbsBasis} → t::Real

Computes the boundary decomposition method tension `t` at wavenumber `k` for
`basis` on the boundary points `pts`.

## Description
The matrices `F` and `G` are built with [`construct_matrices`](@ref), and the
generalized eigenvalues `mu` of `F * x = λ * G * x` are computed with
[`generalized_eigvals`](@ref) (truncated using `solver.eps`). The tension is
`t = 1 / mu[end]`, where `mu[end]` is the largest generalized eigenvalue.

## Arguments
* `solver`: The [`DecompositionMethodSolver`](@ref) used to solve the eigenvalue problem.
* `basis`: The basis used to approximate the eigenstate.
* `pts`: The [`BoundaryPoints`](@ref) with sampled boundary points and quadrature weights `ds`, `w_n`.
* `k`: The wavenumber at which the tension is evaluated.

## Keyword arguments
* `multithreaded::Bool = true`: Whether the matrix construction is multithreaded.

## Returns
* `t`: The tension at wavenumber `k`, `t = 1 / mu[end]`.
"""
function solve(solver::DecompositionMethodSolver,basis::Ba,pts::BoundaryPoints,k;multithreaded::Bool=true) where {Ba<:AbsBasis}
    F,G=construct_matrices(solver,basis,pts,k;multithreaded=multithreaded)
    @blas_multi_then_1 MAX_BLAS_THREADS mu=generalized_eigvals(Symmetric(F),Symmetric(G);eps=solver.eps)
    lam0=mu[end]
    t=1.0/lam0
    return t
end

"""
    solve(solver::DecompositionMethodSolver,F,G) → t::Real

Computes the boundary decomposition method tension `t` directly from precomputed
matrices `F` and `G` (see [`construct_matrices`](@ref)), instead of constructing
them from a basis and boundary points. See [`solve`](@ref) for details.

## Arguments
* `solver`: The [`DecompositionMethodSolver`](@ref) used to solve the eigenvalue problem.
* `F`: Precomputed `F` matrix, see [`construct_matrices`](@ref).
* `G`: Precomputed `G` matrix, see [`construct_matrices`](@ref).

## Returns
* `t`: The tension, `t = 1 / mu[end]`, where `mu` are the generalized eigenvalues of `F * x = λ * G * x`.
"""
function solve(solver::DecompositionMethodSolver,F,G)
    @blas_multi_then_1 MAX_BLAS_THREADS mu=generalized_eigvals(Symmetric(F),Symmetric(G);eps=solver.eps)
    lam0=mu[end]
    t=1.0/lam0
    return t
end

"""
    solve_vect(solver::DecompositionMethodSolver,basis::AbsBasis,pts::BoundaryPoints,k;multithreaded::Bool=true) → (t::Real,x::Vector)

Computes the boundary decomposition method tension `t` and the corresponding
eigenvector `x` (expressed in the original basis) at wavenumber `k`.

## Description
The matrices `F` and `G` are built with [`construct_matrices`](@ref), and the
generalized eigenproblem `F * x = λ * G * x` is solved with
[`generalized_eigen`](@ref) (truncated using `solver.eps`) to obtain the largest
eigenvalue `mu[end]`, its eigenvector `Z[:,end]` in the reduced space, and the
change-of-basis matrix `C`. The eigenvector is transformed back into the original
basis as `x = C * Z[:,end]`, normalized by `sqrt(mu[end])`, and the tension is
`t = 1 / mu[end]`.

## Arguments
* `solver`: The [`DecompositionMethodSolver`](@ref) used to solve the eigenvalue problem.
* `basis`: The basis used to approximate the eigenstate.
* `pts`: The [`BoundaryPoints`](@ref) with sampled boundary points and quadrature weights `ds`, `w_n`.
* `k`: The wavenumber at which the eigenstate is evaluated.

## Keyword arguments
* `multithreaded::Bool = true`: Whether the matrix construction is multithreaded.

## Returns
* `t`: The tension at wavenumber `k`, `t = 1 / mu[end]`.
* `x`: The eigenvector expressed in the original basis, normalized by `sqrt(mu[end])`.
"""
function solve_vect(solver::DecompositionMethodSolver,basis::AbsBasis,pts::BoundaryPoints,k;multithreaded::Bool=true)
    F,G=construct_matrices(solver,basis,pts,k;multithreaded=multithreaded)
    @blas_multi MAX_BLAS_THREADS mu,Z,C=generalized_eigen(Symmetric(F),Symmetric(G);eps=solver.eps)
    x=Z[:,end]
    @blas_multi_then_1 MAX_BLAS_THREADS x=C*x
    lam0=mu[end]
    t=1.0/lam0
    return t,x./sqrt(lam0)
end