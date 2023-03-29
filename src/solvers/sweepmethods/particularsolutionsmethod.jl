using LinearAlgebra, StaticArrays, TimerOutputs

struct ParticularSolutionsMethod{T} <: SweepSolver where T<:Real
    dim_scaling_factor::T
    pts_scaling_factor::T
    int_pts_scaling_factor::T
    eps::T
end

ParticularSolutionsMethod(d,b,b_int) = ParticularSolutionsMethod(d,b,b_int,eps(typeof(d)))
ParticularSolutionsMethod(d,b) = ParticularSolutionsMethod(d,b, 0.2*b)

struct PointsPSM{T} <: AbsPoints where {T<:Real}
    xy_boundary::Vector{SVector{2,T}}
    xy_interior::Vector{SVector{2,T}} #normal vectors in points
end

function evaluate_points(solver::ParticularSolutionsMethod, billiard::Bi, sampler::Function, k) where {Bi<:AbsBilliard}
    b = solver.pts_scaling_factor
    b_int = solver.int_pts_scaling_factor
    type = typeof(solver.pts_scaling_factor)
    xy_all = Vector{SVector{2,type}}()
    xy_int_all = Vector{SVector{2,type}}()

    for crv in billiard.boundary
        if typeof(crv) <: AbsRealCurve
            L = crv.length
            N = round(Int, k*L*b/(2*pi))
            t, dt = sampler(N)
            xy = curve(crv,t)
            append!(xy_all, xy)
        end
    end
    L = real_length(billiard)
    M = round(Int, k*L*b_int/(2*pi))
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
    B, B_int = construct_matrices(solver,basis, pts, k)
    solution = svdvals(B, B_int)
    return minimum(solution)
end

function solve(solver::ParticularSolutionsMethod, B, B_int)
    solution = svdvals(B, B_int)
    return minimum(solution)
end

function solve_vect(solver::ParticularSolutionsMethod,basis::Ba, pts::PointsPSM, k) where {Ba<:AbsBasis}
    B, B_int = construct_matrices(solver,basis, pts, k)
    F = svd(B, B_int)
    H = F.R*F.Q'
    idx = 1:F.k + F.l #inidices containing the singular values we need
    sv = F.a[idx] ./ F.b[idx] #generalized singular values
    X = H[idx,:]
    i_min = argmin(sv)
    return  sv[i_min], X[i_min,:] 
end