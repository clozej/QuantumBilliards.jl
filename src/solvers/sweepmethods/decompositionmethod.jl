#include("../../abstracttypes.jl")
#include("../decompositions.jl")
#include("../matrixconstructors.jl")
#include("../samplers.jl")
#include("../src/billiards/triangle.jl")
#include("../../utils/benchmarkutils.jl")
using LinearAlgebra, StaticArrays, TimerOutputs

"""
A struct for the decomposition method solver.

# Description
`DecompositionMethod` stores parameters for a solver that uses decomposition methods. It includes scaling factors, samplers, precision, and minimum dimensions.

# Fields
- `dim_scaling_factor::T`: Scaling factor for the problem's dimensionality.
- `pts_scaling_factor::Vector{T}`: Scaling factors for the points.
- `sampler::Vector`: A vector of samplers used for point generation.
- `eps::T`: Precision used for numerical operations.
- `min_dim::Int64`: Minimum allowed dimensionality.
- `min_pts::Int64`: Minimum number of points.
"""
struct DecompositionMethod{T} <: SweepSolver where {T<:Real}
    dim_scaling_factor::T
    pts_scaling_factor::Vector{T}
    sampler::Vector
    eps::T
    min_dim::Int64
    min_pts::Int64
end

"""
Create a `DecompositionMethod` solver with default sampler.

# Description
Constructs a `DecompositionMethod` solver with specified dimensional and point scaling factors, along with default parameters for minimum dimension and minimum points.
- The default sampler for points is the Gauss - Legendre for all -> that is why it is implemented as a vector object as per the struct's definition.

# Arguments
- `dim_scaling_factor::T`: Scaling factor for the problem's dimensionality.
- `pts_scaling_factor::Union{T,Vector{T}}`: Scaling factor(s) for the points.
- `min_dim::Int64`: Minimum allowed dimensionality (default `100`).
- `min_pts::Int64`: Minimum number of points (default `500`).

# Returns
- A `DecompositionMethod` object.
"""
function DecompositionMethod(dim_scaling_factor::T, pts_scaling_factor::Union{T,Vector{T}}; min_dim = 100, min_pts = 500) where T<:Real 
    d = dim_scaling_factor
    bs = typeof(pts_scaling_factor) == T ? [pts_scaling_factor] : pts_scaling_factor
    sampler = [GaussLegendreNodes()]
return DecompositionMethod(d, bs, sampler, eps(T), min_dim, min_pts)
end

"""
Create a `DecompositionMethod` solver with a vector of samplers.

# Description
Constructs a `DecompositionMethod` solver with specified dimensional and point scaling factors, along with default parameters for minimum dimension and minimum points.
- The samplers are not the default Gauss-Legendre anymore.

# Arguments
- `dim_scaling_factor::T`: Scaling factor for the problem's dimensionality.
- `pts_scaling_factor::Union{T,Vector{T}}`: Scaling factor(s) for the points.
- `min_dim::Int64`: Minimum allowed dimensionality (default `100`).
- `min_pts::Int64`: Minimum number of points (default `500`).

# Returns
- A `DecompositionMethod` object.
"""
function DecompositionMethod(dim_scaling_factor::T, pts_scaling_factor::Union{T,Vector{T}}, samplers::Vector{AbsSampler}; min_dim = 100, min_pts = 500) where {T<:Real} 
    d = dim_scaling_factor
    bs = typeof(pts_scaling_factor) == T ? [pts_scaling_factor] : pts_scaling_factor
    return DecompositionMethod(d, bs, samplers, eps(T), min_dim, min_pts)
end

"""
Stores boundary points and weights for the decomposition method (DM).

# Description
`BoundaryPointsDM` holds coordinates, normal vectors, and weightings for points on the boundary of a billiard, used specifically for decomposition methods.

# Fields
- `xy::Vector{SVector{2,T}}`: Coordinates of boundary points.
- `normal::Vector{SVector{2,T}}`: Normal vectors at the boundary points.
- `w::Vector{T}`: Tension weights.
- `w_n::Vector{T}`: Normalization weights.
"""
struct BoundaryPointsDM{T} <: AbsPoints where {T<:Real}
    xy::Vector{SVector{2,T}}
    normal::Vector{SVector{2,T}} #normal vectors in points
    w::Vector{T} # tension weights
    w_n::Vector{T} #normalization weights
end

"""
Evaluate boundary points for a given billiard using the decomposition method.

# Description
Generates boundary points, normal vectors, and corresponding weights for use with the `DecompositionMethod` solver.

# Arguments
- `solver::DecompositionMethod`: The decomposition method solver.
- `billiard::Bi`: The billiard on which the points are evaluated.
- `k`: Wavenumber used for evaluation.

# Returns
- A `BoundaryPointsDM` object containing evaluated points, normals and weights.
"""
function evaluate_points(solver::DecompositionMethod, billiard::Bi, k) where {Bi<:AbsBilliard}
    bs, samplers = adjust_scaling_and_samplers(solver, billiard)
    curves = billiard.fundamental_boundary
    type = eltype(solver.pts_scaling_factor)
    xy_all = Vector{SVector{2,type}}()
    normal_all = Vector{SVector{2,type}}()
    w_all = Vector{type}()
    w_n_all = Vector{type}()

    for i in eachindex(curves)
        crv = curves[i]
        if typeof(crv) <: AbsRealCurve
            L = crv.length
            N = max(solver.min_pts,round(Int, k*L*bs[i]/(2*pi)))
            sampler = samplers[i]
            t, dt = sample_points(sampler,N)
            ds = L*dt #modify this
            xy = curve(crv,t)
            normal = normal_vec(crv,t)
            rn = dot.(xy, normal)
            w = ds
            w_n = (w.*rn)./(2.0*k.^2) 
            append!(xy_all, xy)
            append!(normal_all, normal)
            append!(w_all, w)
            append!(w_n_all, w_n)
        end
    end
    return BoundaryPointsDM{type}(xy_all,normal_all, w_all, w_n_all)
end

"""
Construct matrices with benchmarking for the decomposition method.

# Description
Constructs and benchmarks the construction of matrices `F` and `G` used in the decomposition method. Times for each step are recorded and printed.

# Arguments
- `solver::DecompositionMethod`: The decomposition method solver.
- `basis::Ba`: The basis functions used.
- `pts::BoundaryPointsDM`: The boundary points and weights.
- `k`: Wavenumber used for matrix construction.

# Returns
- The matrices `F` and `G`.
"""
function construct_matrices_benchmark(solver::DecompositionMethod, basis::Ba, pts::BoundaryPointsDM, k) where {Ba<:AbsBasis}
    to = TimerOutput()
    w = pts.w
    w_n = pts.w_n
    symmetries=basis.symmetries
    if ~isnothing(symmetries)
        norm = (length(symmetries)+1.0)
        w = w.*norm
        w_n = w_n.*norm
    end
    #basis and gradient matrices
    @timeit to "basis_and_gradient_matrices" B, dX, dY = basis_and_gradient_matrices(basis, k, pts.xy)
    N = basis.dim
    type = eltype(B)
    F = zeros(type,(N,N))
    G = similar(F)
    
    @timeit to "F construction" begin 
        @timeit to "weights" T = (w .* B) #reused later
        @timeit to "product" mul!(F,B',T) #boundary norm matrix
    end

    @timeit to "G construction" begin 
        @timeit to "normal derivative" nx = getindex.(pts.normal,1)
        @timeit to "normal derivative" ny = getindex.(pts.normal,2)
        #inplace modifications
        @timeit to "normal derivative" dX = nx .* dX 
        @timeit to "normal derivative" dY = ny .* dY
        #reuse B
        @timeit to "normal derivative" B = dX .+ dY
    
    #B is now normal derivative matrix (u function)
        @timeit to "weights" T = (w_n .* B) #apply integration weights
        @timeit to "product" mul!(G,B',T)#norm matrix
    end
    print_timer(to)
    return F, G    
end

"""
Construct matrices `F` and `G` for the decomposition method.

# Description
Constructs the matrices `F` and `G` using the boundary points, basis functions, and wavenumber for the `DecompositionMethod` solver.

# Arguments
- `solver::DecompositionMethod`: The decomposition method solver.
- `basis::Ba`: The basis functions used.
- `pts::BoundaryPointsDM`: The boundary points and weights.
- `k`: Wavenumber used for matrix construction.

# Returns
- The matrices `F` and `G`.
"""
function construct_matrices(solver::DecompositionMethod, basis::Ba, pts::BoundaryPointsDM, k) where {Ba<:AbsBasis}
    #basis and gradient matrices
    w = pts.w
    w_n = pts.w_n
    symmetries=basis.symmetries
    if ~isnothing(symmetries)
        norm = (length(symmetries)+1.0)
        w = w.*norm
        w_n = w_n.*norm
    end
    B, dX, dY = basis_and_gradient_matrices(basis, k, pts.xy)
    type = eltype(B)
    #allocate matrices
    N = basis.dim
    F = zeros(type,(N,N))
    G = similar(F)
    #apply weights
    T = (w .* B) #reused later
    mul!(F,B',T) #boundary norm matrix
    nx = getindex.(pts.normal,1)
    ny = getindex.(pts.normal,2)
    #inplace modifications
    dX = nx .* dX 
    dY = ny .* dY
    #reuse B
    B = dX .+ dY
    #B is now normal derivative matrix (u function)
    T = (w_n .* B) #apply integration weights
    mul!(G,B',T)#norm matrix
    return F, G    
end

"""
Solve for the leading eigenvalue using the decomposition method.

# Description
Solves the generalized eigenvalue problem for matrices `F` and `G` to obtain the largest eigenvalue (up to the `solver`'s tolerance).

# Arguments
- `solver::DecompositionMethod`: The decomposition method solver.
- `basis::Ba`: The basis functions used.
- `pts::BoundaryPointsDM`: The boundary points and weights.
- `k`: Wavenumber used for solving.

# Returns
- The computed value `t` corresponding to the inverse of the largest eigenvalue.
"""
function solve(solver::DecompositionMethod,basis::Ba, pts::BoundaryPointsDM, k) where {Ba<:AbsBasis}
    F, G = construct_matrices(solver, basis, pts, k)
    mu = generalized_eigvals(Symmetric(F),Symmetric(G);eps=solver.eps)
    lam0 = mu[end]
    t = 1.0/lam0
    return  t
end

"""
Solve for the largest eigenvalue given matrices `F` and `G` (up to the `solver`'s tolerance).

# Description
Solves the generalized eigenvalue problem for precomputed matrices `F` and `G` to obtain the leading eigenvalue.

# Arguments
- `solver::DecompositionMethod`: The decomposition method solver.
- `F`: The matrix `F`.
- `G`: The matrix `G`.

# Returns
- The computed value `t` corresponding to the inverse of the largest eigenvalue.
"""
function solve(solver::DecompositionMethod,F,G)
    #F, G = construct_matrices(solver, basis, pts, k)
    mu = generalized_eigvals(Symmetric(F),Symmetric(G);eps=solver.eps)
    lam0 = mu[end]
    t = 1.0/lam0
    return  t
end

"""
Solve for the leading eigenvalue and corresponding eigenvector (up to the `solver`'s tolerance).

# Description
Solves the generalized eigenvalue problem for matrices `F` and `G` to obtain the leading eigenvalue and its corresponding eigenvector. The result is transformed back to the original (normalized) basis.

# Arguments
- `solver::DecompositionMethod`: The decomposition method solver.
- `basis::AbsBasis`: The basis functions used.
- `pts::BoundaryPointsDM`: The boundary points and weights.
- `k`: Wavenumber used for solving.

# Returns
- The computed value `t` and the eigenvector `x`.
"""
function solve_vect(solver::DecompositionMethod,basis::AbsBasis, pts::BoundaryPointsDM, k)
    F, G = construct_matrices(solver, basis, pts, k)
    mu, Z, C = generalized_eigen(Symmetric(F),Symmetric(G);eps=solver.eps)
    x = Z[:,end]
    x = C*x #transform into original basis 
    lam0 = mu[end]
    t = 1.0/lam0
    return  t, x./sqrt(lam0)
end
