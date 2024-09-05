#include("../../abstracttypes.jl")
#include("../decompositions.jl")
#include("../matrixconstructors.jl")
#include("../src/billiards/triangle.jl")
#include("../../utils/benchmarkutils.jl")

"""
This module implements various scaling methods and associated functions for solving boundary value problems using accelerated solvers. The key components include:

- `ScalingMethodA` and `ScalingMethodB`: Two different scaling strategies, each defined by a set of parameters and a sampling strategy.
- `BoundaryPointsSM`: A structure to store points on the boundary of a domain and associated weights.
- Methods to evaluate points on a boundary, construct matrices for solving, and solve for eigenvalues and eigenvectors using the defined scaling methods.

The module leverages advanced linear algebra, numerical methods, and benchmarking utilities for performance analysis.
"""

using LinearAlgebra, StaticArrays
using TimerOutputs
abstract type AbsScalingMethod <: AcceleratedSolver 
end

#TODO List of samplers needed
"""
Struct for implementing a scaling method in boundary value problems.

# Fields
- `dim_scaling_factor::T`: Scaling factor for the problem's dimensions.
- `pts_scaling_factor::Vector{T}`: Vector of scaling factors for sampling points.
- `sampler::Vector`: Vector of sampling strategies used for point generation.
- `eps::T`: Numerical tolerance or regularization parameter.
- `min_dim::Int64`: Minimum dimension for the problem.
- `min_pts::Int64`: Minimum number of sampling points.

`T` is a subtype of `Real`, ensuring consistency across all numerical fields.
"""
struct ScalingMethodA{T} <: AbsScalingMethod where {T<:Real}
    dim_scaling_factor::T
    pts_scaling_factor::Vector{T}
    sampler::Vector
    eps::T
    min_dim::Int64
    min_pts::Int64
end

"""
Construct an instance of `ScalingMethodA`. It defines the discretization of the boundary using the input sampler.

# Arguments
- `dim_scaling_factor::T`: A scaling factor for the dimensions, where `T` is a subtype of `Real`.
- `pts_scaling_factor::Union{T,Vector{T}}`: Scaling factor(s) for the points. Can be a single value or a vector of values.
- `samplers::Vector{AbsSampler}`: A vector of sampling strategies (optional). If not given the Gauss Legendre nodes are used. This can be seen from the return statement
- `min_dim::Int64`: Minimum dimension for the problem. Defaults to 100.
- `min_pts::Int64`: Minimum number of points for sampling. Defaults to 500.

# Returns
- An instance of `ScalingMethodA` initialized with the given parameters.
"""
function ScalingMethodA(dim_scaling_factor::T, pts_scaling_factor::Union{T,Vector{T}}; min_dim = 100, min_pts = 500) where T<:Real 
    d = dim_scaling_factor
    bs = typeof(pts_scaling_factor) == T ? [pts_scaling_factor] : pts_scaling_factor
    sampler = [GaussLegendreNodes()]
return ScalingMethodA(d, bs, sampler, eps(T), min_dim, min_pts)
end

function ScalingMethodA(dim_scaling_factor::T, pts_scaling_factor::Union{T,Vector{T}}, samplers::Vector{AbsSampler}; min_dim = 100, min_pts = 500) where {T<:Real} 
    d = dim_scaling_factor
    bs = typeof(pts_scaling_factor) == T ? [pts_scaling_factor] : pts_scaling_factor
    return ScalingMethodA(d, bs, samplers, eps(T), min_dim, min_pts)
end

struct ScalingMethodB{T} <: AbsScalingMethod where {T<:Real}
    dim_scaling_factor::T
    pts_scaling_factor::Vector{T}
    sampler::Vector
    eps::T
    min_dim::Int64
    min_pts::Int64
end

function ScalingMethodB(dim_scaling_factor::T, pts_scaling_factor::Union{T,Vector{T}}; min_dim = 100, min_pts = 500) where T<:Real 
    d = dim_scaling_factor
    bs = typeof(pts_scaling_factor) == T ? [pts_scaling_factor] : pts_scaling_factor
    sampler = [GaussLegendreNodes()]
return ScalingMethodB(d, bs, sampler, eps(T), min_dim, min_pts)
end

function ScalingMethodB(dim_scaling_factor::T, pts_scaling_factor::Union{T,Vector{T}}, samplers::Vector{AbsSampler}; min_dim = 100, min_pts = 500) where {T<:Real} 
    d = dim_scaling_factor
    bs = typeof(pts_scaling_factor) == T ? [pts_scaling_factor] : pts_scaling_factor
    return ScalingMethodB(d, bs, samplers, eps(T), min_dim, min_pts)
end
struct BoundaryPointsSM{T} <: AbsPoints where {T<:Real}
    xy::Vector{SVector{2,T}}
    w::Vector{T}
end

"""
Evaluate points on the boundary of a domain using the given scaling method. This constructs the (x,y) points on the boundary and the diagonal wights for the construction of the the tensions matrix and it's derivative.

# Logic
- We adjust scaling factors (b - point density per wavelength) such that all curve have an adequate number of points
- For each curve in the fundamental boundary of the billiard we check if it is a subtype of AbsRealCurve (like LineSegment, CircleSegments, etc.) because they have curve and normal_vec methods that gets (x,y) boundary coords based on the bs from the previous points and the normal vector for the points (based on the sampler of the AbsScalingMethod)
- For each point (x,y) one the boundary also calcualte the weights for the F and dFdk matrix construction (1/r_n)

# Arguments
- `solver::AbsScalingMethod`: The scaling method to use for evaluation.
- `billiard::Bi`: The billiard domain, where `Bi` is a subtype of `AbsBilliard`.
- `k`: A parameter related to the wave number or frequency.

# Returns
- `BoundaryPointsSM{type}`: A structure containing the evaluated boundary points and associated weights. The struct is defined by the 
"""
function evaluate_points(solver::AbsScalingMethod, billiard::Bi, k) where {Bi<:AbsBilliard}
    bs, samplers = adjust_scaling_and_samplers(solver, billiard)
    curves = billiard.fundamental_boundary
    type = eltype(solver.pts_scaling_factor)
    xy_all = Vector{SVector{2,type}}()
    w_all = Vector{type}()
    
    for i in eachindex(curves)
        crv = curves[i]
        if typeof(crv) <: AbsRealCurve
            L = crv.length
            N = max(solver.min_pts,round(Int, k*L*bs[i]/(2*pi)))
            sampler = samplers[i]
            t, dt = sample_points(sampler, N)
            
            ds = L*dt #this needs modification!!!
            xy = curve(crv,t)
            normal = normal_vec(crv,t)
            rn = dot.(xy, normal)
            w = ds ./ rn
            append!(xy_all, xy)
            append!(w_all, w)
        end
    end
    return BoundaryPointsSM{type}(xy_all, w_all)
end

#generalize for other types
"""
Construct matrices required for solving a boundary value problem using `ScalingMethodA`, with benchmarking.

# Arguments
- `solver::ScalingMethodA`: The scaling method instance.
- `basis::Ba`: The basis functions, where `Ba` is a subtype of `AbsBasis`.
- `pts::BoundaryPointsSM`: The evaluated boundary points and their associated weights.
- `k`: A parameter related to the wave number or frequency.

# Returns
- `F`: The boundary norm matrix.
- `Fk`: The derivative matrix, symmetrized for the problem.

# Logic
This function constructs the matrices `F` and `Fk` required for solving boundary value problems:
1. **Symmetry Adjustment**: If symmetries are present in the basis, the weights `w` are adjusted by multiplying with the factor `n`, where `n` is the number of symmetries plus one. This is to ensure that the number of weights is correct when using the desymmetrized/virtual boundaries
2. **Matrix Dimensions**: `N`, the dimension of the basis, is determined.
3. **Basis Matrix (`B`) Construction**: The basis matrix `B` is computed using the `basis_matrix` function.
4. **Boundary Norm Matrix (`F`) Construction**: The matrix `F` is computed by multiplying the transposed basis matrix with the weights, and then multiplying with the basis matrix again. This calculates the F matrix using Barnett's signature: https://users.flatironinstitute.org/~ahb/thesis_html/node80.html : Solving for the scaling eigenfunctions
5. **Derivative Matrix (`Fk`) Construction**: The matrix `Fk` is computed similarly, but using the derivative of the basis functions, obtained through `dk_matrix`.
6. **Symmetrization**: The matrix `Fk` is symmetrized by adding it to its transpose. This calculates the dFdk matrix using Barnett's signature: https://users.flatironinstitute.org/~ahb/thesis_html/node80.html : Solving for the scaling eigenfunctions (take note of the fact that we use the dFdk + transpose like Barnett's work)

The function also uses `TimerOutput` to benchmark the performance of different stages, providing detailed timing information.
"""
function construct_matrices_benchmark(solver::ScalingMethodA, basis::Ba, pts::BoundaryPointsSM, k) where {Ba<:AbsBasis}
    to = TimerOutput()
    symmetries = basis.symmetries 
    #type = eltype(pts.w)
    xy, w = pts.xy, pts.w
    symmetries = basis.symmetries 
    if ~isnothing(symmetries)
        n = (length(symmetries)+1.0)
        w = w.*n
    end
    #M =  length(xy)
    N = basis.dim
    #basis matrix
    @timeit to "basis_matrix" B = basis_matrix(basis, k, xy)
    type = eltype(B)
    F = zeros(type,(N,N))
    Fk = similar(F)
    @timeit to "F construction" begin 
        @timeit to "weights" T = (w.*B) #reused later
        #@timeit to "copy" Bt = copy(B')
        @timeit to "product" mul!(F,B',T) #boundary norm matrix
    end
    #reuse B
    @timeit to "dk_matrix" B = dk_matrix(basis,k, xy)
    @timeit to "Fk construction" begin 
        @timeit to "product" mul!(Fk,B',T) #B is now derivative matrix
        #symmetrize matrix
        @timeit to "addition" Fk = Fk+Fk'
    end
    print_timer(to)    
    return F, Fk        
end

"""
Construct matrices required for solving a boundary value problem using `ScalingMethodA`.

# Arguments
- `solver::ScalingMethodA`: The scaling method instance.
- `basis::Ba`: The basis functions, where `Ba` is a subtype of `AbsBasis`.
- `pts::BoundaryPointsSM`: The evaluated boundary points and their associated weights.
- `k`: A parameter related to the wave number or frequency.

# Returns
- `F`: The boundary norm matrix.
- `Fk`: The derivative matrix, symmetrized for the problem.

# Logic
This function constructs the matrices `F` and `Fk` required for solving boundary value problems:

1. **Symmetry Adjustment**: If symmetries are present in the basis, the weights `w` are adjusted by multiplying with the factor `n`, where `n` is the number of symmetries plus one. This is to ensure that the number of weights is correct when using the desymmetrized/virtual boundaries
2. **Matrix Dimensions**: `N`, the dimension of the basis, is determined.
3. **Basis Matrix (`B`) Construction**: The basis matrix `B` is computed using the `basis_matrix` function.
4. **Boundary Norm Matrix (`F`) Construction**: The matrix `F` is computed by multiplying the transposed basis matrix with the weights, and then multiplying with the basis matrix again. This calculates the F matrix using Barnett's signature: https://users.flatironinstitute.org/~ahb/thesis_html/node80.html : Solving for the scaling eigenfunctions
5. **Derivative Matrix (`Fk`) Construction**: The matrix `Fk` is computed similarly, but using the derivative of the basis functions, obtained through `dk_matrix`.
6. **Symmetrization**: The matrix `Fk` is symmetrized by adding it to its transpose. This calculates the dFdk matrix using Barnett's signature: https://users.flatironinstitute.org/~ahb/thesis_html/node80.html : Solving for the scaling eigenfunctions (take note of the fact that we use the dFdk + transpose like Barnett's work)

This function is a streamlined version focused solely on matrix construction without the overhead of performance benchmarking.
"""
function construct_matrices(solver::ScalingMethodA, basis::Ba, pts::BoundaryPointsSM, k) where {Ba<:AbsBasis}
    xy = pts.xy
    w = pts.w
    symmetries=basis.symmetries
    if ~isnothing(symmetries)
        n = (length(symmetries)+1.0)
        w = w.*n
    end
    N = basis.dim
    #basis matrix
    B = basis_matrix(basis, k, xy)
    type = eltype(B)
    F = zeros(type,(N,N))
    Fk = similar(F)
    T = (w .* B) #reused later
    mul!(F,B',T) #boundary norm matrix
    #reuse B
    B = dk_matrix(basis,k, xy)
    mul!(Fk,B',T) #B is now derivative matrix
    #symmetrize matrix
    Fk = Fk + Fk' 
    return F, Fk    
end

function construct_matrices_benchmark(solver::ScalingMethodB, basis::Ba, pts::BoundaryPointsSM, k) where {Ba<:AbsBasis}
    to = TimerOutput()
    #type = eltype(pts.w)
    xy, w = pts.xy, pts.w
    symmetries=basis.symmetries
    if ~isnothing(symmetries)
        n = (length(symmetries)+1.0)
        w = w.*n
    end
    #M =  length(xy)
    #basis and gradient matrices
    @timeit to "basis_and_gradient_matrices" B, dX, dY = basis_and_gradient_matrices(basis, k, xy)
    N = basis.dim
    type = eltype(B)
    F = zeros(type,(N,N))
    Fk = similar(F)
    @timeit to "F construction" begin 
        @timeit to "weights" T = (w.*B) #reused later
        #@timeit to "copy" Bt = copy(B')
        @timeit to "product" mul!(F,B',T) #boundary norm matrix
    end
    #reuse B
    @timeit to "Fk construction" begin 
        @timeit to "dilation derivative" x = getindex.(xy,1)
        @timeit to "dilation derivative" y = getindex.(xy,2)
        #inplace modifications
        @timeit to "dilation derivative" dX = x .* dX 
        @timeit to "dilation derivative" dY = y .* dY
        #reuse B
        @timeit to "dilation derivative" B = dX .+ dY
        @timeit to "product" mul!(Fk,B',T) #B is now derivative matrix
        #symmetrize matrix
        @timeit to "addition" Fk = (Fk+Fk') ./ k
    end
    print_timer(to)    
    return F, Fk        
end

function construct_matrices(solver::ScalingMethodB, basis::Ba, pts::BoundaryPointsSM, k) where {Ba<:AbsBasis}
    xy = pts.xy
    w = pts.w
    symmetries=basis.symmetries
    if ~isnothing(symmetries)
        n = (length(symmetries)+1.0)
        w = w.*n
    end
    N = basis.dim
    #basis matrix
    B, dX, dY = basis_and_gradient_matrices(basis, k, pts.xy)
    type = eltype(B)
    F = zeros(type,(N,N))
    Fk = similar(F)
    T = (w .* B) #reused later
    mul!(F,B',T) #boundary norm matrix
    x = getindex.(xy,1)
    y = getindex.(xy,2)
    #inplace modifications
    dX = x .* dX 
    dY = y .* dY
    #reuse B
    B = dX .+ dY
    mul!(Fk,B',T) #B is now derivative matrix
    #symmetrize matrix
    Fk = (Fk+Fk') ./ k
    return F, Fk    
end

"""
Compute the first-order correction to the reference `k` values and determine the tension using Alex Barnett's method.

first order correction: https://users.flatironinstitute.org/~ahb/thesis_html/node81.html : Generalized eigenproblem
tension determination logic: https://users.flatironinstitute.org/~ahb/thesis_html/node60.html : Truncating the singular generalized eigenproblem

# Arguments
- `mu`: The generalized eigenvalues from the previous computation.
- `k`: The reference wave number or frequency.

# Returns
- `ks`: The corrected `k` values incorporating the first-order correction.
- `ten`: The corresponding tension values computed using Barnett's formula.
"""
function sm_results(mu,k)
    ks = k .- 2 ./mu .+ 2/k ./(mu.^2) 
    ten = 2.0 .*(2.0 ./ mu).^2
    return ks, ten
end

#=
function sm_vects_results(mu,k)
    ks = k .- 2 ./mu .+ 2/k ./(mu.^2) 
    ten = 2.0 .*(2.0 ./ mu).^2
    #does not sort the results
    return ks, ten
end
=#

"""
Solve the generalized eigenproblem for potential `k` values in the Vergini-Saraceno method.

    logic: https://users.flatironinstitute.org/~ahb/thesis_html/node81.html : Generalized eigenproblem

# Arguments
- `solver::AbsScalingMethod`: The scaling method used to set up the problem.
- `basis::Ba`: The basis functions, where `Ba` is a subtype of `AbsBasis`.
- `pts::BoundaryPointsSM`: The evaluated boundary points and their associated weights.
- `k`: The reference wave number or frequency.
- `dk`: The tolerance for selecting eigenvalues close to the reference `k`.

# Returns
- `ks`: Sorted vector of corrected `k` values derived from the generalized eigenproblem, filtered by proximity to `k`.
- `ten`: Sorted vector of corresponding tension values.

# Logic
1. **Matrix Construction**: Use the `construct_matrices` function to build the boundary norm matrix `F` and the derivative matrix `Fk`.
2. **Generalized Eigenproblem**: Solve the generalized eigenproblem, `F * v = mu * Fk * v`, to obtain eigenvalues `mu`. These eigenvalues represent potential `k` values.
3. **First-Order Corrections**: Apply the `sm_results` function to calculate the corrected `k` values (`ks`) and associated tensions (`ten`) using the eigenvalues `mu` and the reference `k`.
4. **Filter Results**: Retain the `ks` values that lie within a specified tolerance `dk` of the reference `k`, and filter the corresponding tensions accordingly.
5. **Sort Results**: Sort the corrected `ks` values and their associated tensions in ascending order.
"""
function solve(solver::AbsScalingMethod, basis::Ba, pts::BoundaryPointsSM, k, dk) where {Ba<:AbsBasis}
    F, Fk = construct_matrices(solver, basis, pts, k)
    mu = generalized_eigvals(Symmetric(F),Symmetric(Fk);eps=solver.eps)
    ks, ten = sm_results(mu,k)
    idx = abs.(ks.-k) .< dk
    ks = ks[idx]
    ten = ten[idx]
    p = sortperm(ks)
    return ks[p], ten[p]
end

"""
Solve the generalized eigenproblem for potential `k` values in the Vergini-Saraceno method. Here we explicitely input the matrices F and Fk and do not use the basis for construction

    logic: https://users.flatironinstitute.org/~ahb/thesis_html/node81.html : Generalized eigenproblem

# Arguments
- `solver::AbsScalingMethod`: The scaling method used to set up the problem.
- `F`: The tension matrix as described in the reference
- `Fk`: The derivative matrix as described in the reference
- `pts::BoundaryPointsSM`: The evaluated boundary points and their associated weights.
- `k`: The reference wave number or frequency.
- `dk`: The tolerance for selecting eigenvalues close to the reference `k`.

# Returns
- `ks`: Sorted vector of corrected `k` values derived from the generalized eigenproblem, filtered by proximity to `k`.
- `ten`: Sorted vector of corresponding tension values.

# Logic
1. **Generalized Eigenproblem**: Solve the generalized eigenproblem, `F * v = mu * Fk * v`, to obtain eigenvalues `mu`. These eigenvalues represent potential `k` values.
2. **First-Order Corrections**: Apply the `sm_results` function to calculate the corrected `k` values (`ks`) and associated tensions (`ten`) using the eigenvalues `mu` and the reference `k`.
3. **Filter Results**: Retain the `ks` values that lie within a specified tolerance `dk` of the reference `k`, and filter the corresponding tensions accordingly.
4. **Sort Results**: Sort the corrected `ks` values and their associated tensions in ascending order.
"""
function solve(solver::AbsScalingMethod,F,Fk, k, dk)
    #F, Fk = construct_matrices(solver, basis, pts, k)
    mu = generalized_eigvals(Symmetric(F),Symmetric(Fk);eps=solver.eps)
    ks, ten = sm_results(mu,k)
    idx = abs.(ks.-k) .< dk
    ks = ks[idx]
    ten = ten[idx]
    p = sortperm(ks)
    return ks[p], ten[p]
end

"""
Solve the generalized eigenproblem for potential `k` values in the Vergini-Saraceno method.

    logic: https://users.flatironinstitute.org/~ahb/thesis_html/node81.html : Generalized eigenproblem
    eigenvector construction: https://users.flatironinstitute.org/~ahb/thesis_html/node60.html : Truncating the singular generalized eigenproblem

# Arguments
- `solver::AbsScalingMethod`: The scaling method used to set up the problem.
- `basis::Ba`: The basis functions, where `Ba` is a subtype of `AbsBasis`.
- `pts::BoundaryPointsSM`: The evaluated boundary points and their associated weights.
- `k`: The reference wave number or frequency.
- `dk`: The tolerance for selecting eigenvalues close to the reference `k`.

# Returns
- `ks`: Sorted vector of corrected `k` values derived from the generalized eigenproblem, filtered by proximity to `k`.
- `ten`: Sorted vector of corresponding tension values.

# Logic
1. **Matrix Construction**: Use the `construct_matrices` function to build the boundary norm matrix `F` and the derivative matrix `Fk`.
2. **Generalized Eigenproblem**: Solve the generalized eigenproblem, `F * v = mu * Fk * v`, to obtain eigenvalues `mu`. These eigenvalues represent potential `k` values.
3. **First-Order Corrections**: Apply the `sm_results` function to calculate the corrected `k` values (`ks`) and associated tensions (`ten`) using the eigenvalues `mu` and the reference `k`.
4. **Filter Results**: Retain the `ks` values that lie within a specified tolerance `dk` of the reference `k`, and filter the corresponding tensions accordingly.
5. **Eigenvector construction**: Construct the eigenvectors for the ks (tension values) that survive the specified tolerance. For further reference check the link
5. **Sort Results**: Sort the corrected `ks` values and their associated tensions in ascending order.
"""
function solve_vectors(solver::AbsScalingMethod, basis::Ba, pts::BoundaryPointsSM, k, dk) where {Ba<:AbsBasis}
    F, Fk = construct_matrices(solver, basis, pts, k)
    mu, Z, C = generalized_eigen(Symmetric(F),Symmetric(Fk);eps=solver.eps)
    ks, ten = sm_results(mu,k)
    idx = abs.(ks.-k) .< dk
    ks = ks[idx]
    ten = ten[idx]
    Z = Z[:,idx]
    X = C*Z #transform into original basis 
    X = (sqrt.(ten))' .* X
    p = sortperm(ks)
    return  ks[p], ten[p], X[:,p]
end

