using LinearAlgebra, StaticArrays, TimerOutputs

#TODO More info: Reviving the Method of Particular Solutions, Timo Betcke & Lloyd N. Trefethen

"""
A struct for the Particular Solutions Method solver.

# Description
`ParticularSolutionsMethod` stores parameters for a solver that uses the Particular Solutions Method. It includes scaling factors for boundary and interior points, samplers, precision, and minimum dimensions.

# Fields
- `dim_scaling_factor::T`: Scaling factor for the problem's dimensionality.
- `pts_scaling_factor::Vector{T}`: Scaling factors for boundary points.
- `int_pts_scaling_factor::T`: Scaling factor for interior points.
- `sampler::Vector`: A vector of samplers used for point generation.
- `eps::T`: Precision used for numerical operations.
- `min_dim::Int64`: Minimum allowed dimensionality.
- `min_pts::Int64`: Minimum number of boundary points.
- `min_int_pts::Int64`: Minimum number of interior points.
"""
struct ParticularSolutionsMethod{T} <: SweepSolver where {T<:Real}
    dim_scaling_factor::T
    pts_scaling_factor::Vector{T}
    int_pts_scaling_factor::T
    sampler::Vector
    eps::T
    min_dim::Int64
    min_pts::Int64
    min_int_pts::Int64
end

"""
Create a `ParticularSolutionsMethod` solver with default samplers. This is just a vector of Gauss-Legendre default samplers.

# Description
Constructs a `ParticularSolutionsMethod` solver with specified dimensional scaling factors for boundary and interior points, along with default parameters for minimum dimensions and points.

# Arguments
- `dim_scaling_factor::T`: Scaling factor for the problem's dimensionality.
- `pts_scaling_factor::Union{T,Vector{T}}`: Scaling factor(s) for boundary points.
- `int_pts_scaling_factor::T`: Scaling factor for interior points.
- `min_dim::Int64`: Minimum allowed dimensionality (default `100`).
- `min_pts::Int64`: Minimum number of boundary points (default `500`).
- `min_int_pts::Int64`: Minimum number of interior points (default `500`).

# Returns
- A `ParticularSolutionsMethod` object.
"""
function ParticularSolutionsMethod(dim_scaling_factor::T, pts_scaling_factor::Union{T,Vector{T}}, int_pts_scaling_factor::T; min_dim = 100, min_pts = 500, min_int_pts=500) where T<:Real 
    d = dim_scaling_factor
    bs = typeof(pts_scaling_factor) == T ? [pts_scaling_factor] : pts_scaling_factor
    sampler = [GaussLegendreNodes()]
    return ParticularSolutionsMethod(d, bs,int_pts_scaling_factor, sampler, eps(T), min_dim, min_pts, min_int_pts)
end

"""
Create a `ParticularSolutionsMethod` solver with inputted vector of samplers.

# Description
Constructs a `ParticularSolutionsMethod` solver with specified dimensional scaling factors for boundary and interior points, along with default parameters for minimum dimensions and points.

# Arguments
- `dim_scaling_factor::T`: Scaling factor for the problem's dimensionality.
- `pts_scaling_factor::Union{T,Vector{T}}`: Scaling factor(s) for boundary points.
- `int_pts_scaling_factor::T`: Scaling factor for interior points.
- `samplers::Vector{AbsSampler}`: A vector of samplers used for point generation
- `min_dim::Int64`: Minimum allowed dimensionality (default `100`).
- `min_pts::Int64`: Minimum number of boundary points (default `500`).
- `min_int_pts::Int64`: Minimum number of interior points (default `500`).

# Returns
- A `ParticularSolutionsMethod` object.
"""
function ParticularSolutionsMethod(dim_scaling_factor::T, pts_scaling_factor::Union{T,Vector{T}}, int_pts_scaling_factor::T, samplers::Vector{AbsSampler}; min_dim = 100, min_pts = 500, min_int_pts=500) where T<:Real 
    d = dim_scaling_factor
    bs = typeof(pts_scaling_factor) == T ? [pts_scaling_factor] : pts_scaling_factor
    return ParticularSolutionsMethod(d, bs,int_pts_scaling_factor, samplers, eps(T), min_dim, min_pts, min_int_pts)
end

"""
Stores boundary and interior points for the Particular Solutions Method.

# Description
`PointsPSM` holds coordinates for boundary and interior points used specifically for the Particular Solutions Method.

# Fields
- `xy_boundary::Vector{SVector{2,T}}`: Coordinates of boundary points.
- `xy_interior::Vector{SVector{2,T}}`: Coordinates of interior points.
"""
struct PointsPSM{T} <: AbsPoints where {T<:Real}
    xy_boundary::Vector{SVector{2,T}}
    xy_interior::Vector{SVector{2,T}} #normal vectors in points
end

"""
Evaluate boundary and interior points for a given billiard using the Particular Solutions Method.

# Description
Generates boundary and interior points for use with the `ParticularSolutionsMethod` solver, based on the provided wavenumber `k`.

# Arguments
- `solver::ParticularSolutionsMethod`: The Particular Solutions Method solver.
- `billiard::Bi`: The billiard instance on which the points are evaluated.
- `k`: Wavenumber used for evaluation.

# Returns
- A `PointsPSM` object containing the evaluated boundary and interior points.
"""
function evaluate_points(solver::ParticularSolutionsMethod, billiard::Bi, k) where {Bi<:AbsBilliard}
    bs, samplers = adjust_scaling_and_samplers(solver, billiard)
    b_int = solver.int_pts_scaling_factor
    curves = billiard.fundamental_boundary
    type = eltype(solver.pts_scaling_factor)
    xy_all = Vector{SVector{2,type}}()
    xy_int_all = Vector{SVector{2,type}}()
    
    for i in eachindex(curves)
        crv = curves[i]
        if typeof(crv) <: AbsRealCurve
            L = crv.length
            N = max(solver.min_pts,round(Int, k*L*bs[i]/(2*pi)))
            sampler = samplers[i]
            t, dt = sample_points(sampler, N)
            xy = curve(crv,t)
            append!(xy_all, xy)
        end
    end
    L = billiard.length
    M = max(solver.min_int_pts,round(Int, k*L*b_int/(2*pi)))
    xy_int_all = random_interior_points(billiard,M)
    return PointsPSM{type}(xy_all, xy_int_all)
end

"""
Construct matrices with benchmarking for the Particular Solutions Method.

# Description
Constructs and benchmarks the construction of matrices `B` and `B_int` used in the Particular Solutions Method. The time taken for each step is recorded and printed.

# Arguments
- `solver::ParticularSolutionsMethod`: The Particular Solutions Method solver.
- `basis::Ba`: The basis functions used.
- `pts::PointsPSM`: The boundary and interior points.
- `k`: Wavenumber used for matrix construction.

# Returns
- The matrices `B` and `B_int`.
"""
function construct_matrices_benchmark(solver::ParticularSolutionsMethod, basis::Ba, pts::PointsPSM, k) where {Ba<:AbsBasis}
    to = TimerOutput()
    pts_bd = pts.xy_boundary
    pts_int = pts.xy_interior
    #basis and gradient matrices
    @timeit to "basis_matrices" begin
        @timeit to "boundary" B = basis_matrix(basis,k,pts_bd)
        @timeit to "interior" B_int = basis_matrix(basis,k,pts_int)
    end
    print_timer(to)
    return B, B_int  
end

"""
Construct matrices `B` and `B_int` for the Particular Solutions Method.

# Description
Constructs the matrices `B` and `B_int` using the boundary and interior points, basis functions, and wavenumber for the `ParticularSolutionsMethod` solver.

# Arguments
- `solver::ParticularSolutionsMethod`: The Particular Solutions Method solver.
- `basis::Ba`: The basis functions used.
- `pts::PointsPSM`: The boundary and interior points.
- `k`: Wavenumber used for matrix construction.

# Returns
- The matrices `B` and `B_int`.
"""
function construct_matrices(solver::ParticularSolutionsMethod, basis::Ba, pts::PointsPSM, k) where {Ba<:AbsBasis}
    pts_bd = pts.xy_boundary
    pts_int = pts.xy_interior
    B = basis_matrix(basis,k,pts_bd)
    B_int = basis_matrix(basis,k,pts_int)
    return B, B_int  
end

"""
Solve for the smallest singular value using the Particular Solutions Method.

# Description
Solves for the smallest singular value using the SVD of the matrices `B` and `B_int` constructed for the `ParticularSolutionsMethod`.
- The matrix B and B_int is done internally using the `construct_matrices` method
# Arguments
- `solver::ParticularSolutionsMethod`: The Particular Solutions Method solver.
- `basis::Ba`: The basis functions used.
- `pts::PointsPSM`: The boundary and interior points.
- `k`: Wavenumber used for solving.

# Returns
- The minimum singular value.
"""
function solve(solver::ParticularSolutionsMethod, basis::Ba, pts::PointsPSM, k) where {Ba<:AbsBasis}
    B, B_int = construct_matrices(solver, basis, pts, k)
    solution = svdvals(B, B_int)
    return minimum(solution)
end

"""
Solve for the smallest singular value given matrices `B` and `B_int`.

# Description
Solves for the smallest singular value using the SVD of the precomputed matrices `B` and `B_int` for the `ParticularSolutionsMethod`.

# Arguments
- `solver::ParticularSolutionsMethod`: The Particular Solutions Method solver.
- `B`: The boundary points matrix `B`.
- `B_int`: The interior points matrix `B_int`.

# Returns
- The minimum singular value.
"""
function solve(solver::ParticularSolutionsMethod, B, B_int)
    solution = svdvals(B, B_int)
    return minimum(solution)
end

"""
Solve for the smallest singular value and corresponding vector using the Particular Solutions Method.
For futher reading and explanation of the logic please follow Betcke's paper: Reviving the Method of Particular Solutions, Timo Betcke & Lloyd N. Trefethen

# Description
Solves for the smallest singular value and its corresponding vector using the SVD of the matrices `B` and `B_int` constructed for the `ParticularSolutionsMethod`. The result is returned as both the singular value and the corresponding vector.

# Arguments
- `solver::ParticularSolutionsMethod`: The Particular Solutions Method solver.
- `basis::Ba`: The basis functions used.
- `pts::PointsPSM`: The boundary and interior points.
- `k`: Wavenumber used for solving.

# Returns
- The minimum singular value and the corresponding vector.
"""
function solve_vect(solver::ParticularSolutionsMethod, basis::Ba, pts::PointsPSM, k) where {Ba<:AbsBasis}
    B, B_int = construct_matrices(solver, basis, pts, k)
    F = svd(B, B_int)
    H = F.R*F.Q'
    idx = 1:F.k + F.l #inidices containing the singular values we need
    sv = F.a[idx] ./ F.b[idx] #generalized singular values
    X = H[idx,:]
    i_min = argmin(sv)
    return  sv[i_min], X[i_min,:] 
end