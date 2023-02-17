include("../../abstracttypes.jl")
include("../../billiards/billiard.jl")
include("../samplers.jl")
#include("../src/billiards/triangle.jl")
using LinearAlgebra, Optim
abstract type AbsPSMSolver <: SweepSolver end

struct ParticularSolutionsMethod <: AbsPSMSolver
    dim_scaling_factor::Float64
    pts_scaling_factor::Float64
    int_pts_scaling_factor::Float64
    eps:: Float64
end

struct DistributedParticularSolutionsMethod <: AbsPSMSolver
    dim_scaling_factor::Vector{Float64}
    pts_scaling_factor::Vector{Float64}
    int_pts_scaling_factor::Float64
    eps:: Float64
end

struct AcceleratedParticularSolutionsMethod <: AbsPSMSolver
    dim_scaling_factor::Float64
    pts_scaling_factor::Float64
    int_pts_scaling_factor::Float64
    eps:: Float64
end

ParticularSolutionsMethod(d,b,b_int) = ParticularSolutionsMethod(d,b,b_int,1e-15)
ParticularSolutionsMethod(d::Vector{Float64},b::Vector{Float64},b_int::Float64) = DistributedParticularSolutionsMethod(d,b,b_int,1e-15)
ParticularSolutionsMethod(d) = ParticularSolutionsMethod(d,d,d,1e-15)
ParticularSolutionsMethod(d::Vector{Float64},b_int::Float64) = DistributedParticularSolutionsMethod(d,d,b_int,1e-15)

struct PSMPoints <: AbsPoints
    x::Vector{Float64}
    y::Vector{Float64}
    x_int::Vector{Float64}
    y_int::Vector{Float64}
end


function evaluate_points(solver::ParticularSolutionsMethod, curve::AbsCurve, sampler::Function, k)
    b = solver.pts_scaling_factor
    L = curve.length
    N = round(Int, k*L*b/(2*pi))
    t, _ = sampler(N)
    x, y = curve.r(t)
    return x, y
end


function evaluate_points(solver::ParticularSolutionsMethod, billiard::AbsBilliard, sampler::Function, k)
    b = solver.pts_scaling_factor
    b_int = solver.int_pts_scaling_factor
    x_all = Float64[]
    y_all = Float64[]
    N_int = round(Int, k*real_length(billiard)*b_int/(2*pi))
    x_int, y_int = random_interior_points(billiard, N_int; grd = 1000)

    for curve in billiard.boundary
        if typeof(curve) <: AbsRealCurve
            L = curve.length
            N = round(Int, k*L*b/(2*pi))
            t, _ = sampler(N)
            x, y = curve.r(t)
            append!(x_all, x)
            append!(y_all, y)
        end
    end
    return PSMPoints(x_all, y_all, x_int, y_int)
end

function evaluate_points(solver::DistributedParticularSolutionsMethod, billiard::AbsBilliard, sampler::Function, k)
    b = solver.pts_scaling_factor
    b_int = solver.int_pts_scaling_factor
    pts_vec = []#Vector{PSMPoints}
    N_int = round(Int, k*real_length(billiard)*b_int/(2*pi))
    x_int, y_int = random_interior_points(billiard, N_int; grd = 1000)

    for (i,curve) in enumerate(billiard.boundary)
        if typeof(curve) <: AbsRealCurve
            b = solver.pts_scaling_factor[i]
            L = curve.length
            N = round(Int, k*L*b/(2*pi))
            t, _ = sampler(N)
            x, y = curve.r(t)
            push!(pts_vec, PSMPoints(x, y, x_int, y_int))
        end
    end
    return pts_vec
end

function construct_matrices(solver::ParticularSolutionsMethod, basis::AbsBasis, pts::PSMPoints, k)#{T} where T<:Real
    x, y, x_int, y_int = pts.x, pts.y, pts.x_int, pts.y_int
    N = basis.dim
    M = length(x)
    M_int = length(x_int)
    B = Array{Float64}(undef,M,N)  #basis matrix
    B_int = Array{Float64}(undef,M_int,N)
    for i in 1:N
        B[:,i] = basis_fun(basis, i, k, x, y)
        B_int[:,i] = basis_fun(basis, i, k, x_int, y_int)
    end 
    
    return B, B_int    
end

function construct_matrices(solver::DistributedParticularSolutionsMethod, basis_vec, pts_vec, k)#{T} where T<:Real
    sl = ParticularSolutionsMethod(1.0) #only needed for dispatch otherwise irrelevant. Find better solution.
    B_col = []
    B_int_col = []
    for (basis,pts) in zip(basis_vec,pts_vec)
        B, B_int = construct_matrices(sl, basis, pts, k)
        push!(B_col,B)
        push!(B_int_col,B_int)
    end
    return reduce(directsum, B_col), reduce(hcat, B_int_col)
end

function solve(solver::ParticularSolutionsMethod, basis::AbsBasis, pts::PSMPoints, k)
    B, B_int = construct_matrices(solver,basis, pts, k)
    solution = svdvals(B, B_int)
    return minimum(solution)
end

function solve(solver::DistributedParticularSolutionsMethod, basis::AbsBasis, pts::PSMPoints, k)
    B, B_int = construct_matrices(solver,basis, pts, k)
    solution = svdvals(B, B_int)
    return minimum(solution)
end

#=
function solve_wavenumber(solver::ParticularSolutionsMethod, basis::AbsBasis, pts::PSMPoints, k, dk)
    function f(k)
        B, B_int = construct_matrices(solver,basis, pts, k)
        solution = svdvals(B, B_int)
        return minimum(solution)
    end
    res =  optimize(f, k-0.5*dk, k+0.5*dk)
    k0,t0 = res.minimizer, res.minimum
    return k0, t0
end

function solve_wavenumber(solver::DistributedParticularSolutionsMethod, basis_vec, pts_vec, k, dk)
    function f(k)
        B, B_int = construct_matrices(solver, basis_vec, pts_vec, k)
        solution = svdvals(B, B_int)
        return minimum(solution)
    end
    res =  optimize(f, k-0.5*dk, k+0.5*dk)
    k0,t0 = res.minimizer, res.minimum
    return k0, t0
end

function k_sweep(solver::ParticularSolutionsMethod, basis::AbsBasis, pts::PSMPoints, ks)
    res = similar(ks)
    for (i,k) in enumerate(ks)
        B, B_int = construct_matrices(solver,basis, pts, k)
        solution = svdvals(B, B_int)
        res[i] = minimum(solution)
    end
    return res
end

function k_sweep(solver::DistributedParticularSolutionsMethod,  basis_vec, pts_vec, ks)
    res = similar(ks)
    for (i,k) in enumerate(ks)
        B, B_int = construct_matrices(solver, basis_vec, pts_vec, k)
        solution = svdvals(B, B_int)
        res[i] = minimum(solution)
    end
    return res
end
=#