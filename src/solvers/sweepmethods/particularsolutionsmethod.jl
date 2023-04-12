using LinearAlgebra, StaticArrays, TimerOutputs

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

function ParticularSolutionsMethod(dim_scaling_factor::T, pts_scaling_factor::Union{T,Vector{T}}, int_pts_scaling_factor::T; min_dim = 100, min_pts = 500, min_int_pts=500) where T<:Real 
    d = dim_scaling_factor
    bs = typeof(pts_scaling_factor) == T ? [pts_scaling_factor] : pts_scaling_factor
    sampler = [GaussLegendreNodes()]
    return ParticularSolutionsMethod(d, bs,int_pts_scaling_factor, sampler, eps(T), min_dim, min_pts, min_int_pts)
end

function ParticularSolutionsMethod(dim_scaling_factor::T, pts_scaling_factor::Union{T,Vector{T}}, int_pts_scaling_factor::T, samplers::Vector{AbsSampler}; min_dim = 100, min_pts = 500, min_int_pts=500) where T<:Real 
    d = dim_scaling_factor
    bs = typeof(pts_scaling_factor) == T ? [pts_scaling_factor] : pts_scaling_factor
    return ParticularSolutionsMethod(d, bs,int_pts_scaling_factor, samplers, eps(T), min_dim, min_pts, min_int_pts)
end

struct PointsPSM{T} <: AbsPoints where {T<:Real}
    xy_boundary::Vector{SVector{2,T}}
    xy_interior::Vector{SVector{2,T}} #normal vectors in points
end

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

function construct_matrices(solver::ParticularSolutionsMethod, basis::Ba, pts::PointsPSM, k) where {Ba<:AbsBasis}
    pts_bd = pts.xy_boundary
    pts_int = pts.xy_interior
    B = basis_matrix(basis,k,pts_bd)
    B_int = basis_matrix(basis,k,pts_int)
    return B, B_int  
end

function solve(solver::ParticularSolutionsMethod, basis::Ba, pts::PointsPSM, k) where {Ba<:AbsBasis}
    B, B_int = construct_matrices(solver, basis, pts, k)
    solution = svdvals(B, B_int)
    return minimum(solution)
end

function solve(solver::ParticularSolutionsMethod, B, B_int)
    solution = svdvals(B, B_int)
    return minimum(solution)
end

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