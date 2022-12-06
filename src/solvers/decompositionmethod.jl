include("../abstracttypes.jl")
include("decompositions.jl")
include("samplers.jl")
#include("../src/billiards/triangle.jl")
using LinearAlgebra

struct DecompositionMethod <: SweepSolver 
    dim_scaling_factor::Float64
    pts_scaling_factor::Float64
    eps::Float64
end

DecompositionMethod(dim_scaling_factor, pts_scaling_factor) = DecompositionMethod(dim_scaling_factor, pts_scaling_factor, 10.0*eps())

struct DMPoints{T<:Number} <: AbsPoints 
    x::Vector{T}
    y::Vector{T}
    n_x::Vector{T}
    n_y::Vector{T}
    w::Vector{T} # tension weights
    w_n::Vector{T} #normalization weights
end

function evaluate_points(solver::DecompositionMethod, billiard::AbsBilliard, sampler::Function, k)
    b = solver.pts_scaling_factor
    x_all = Float64[]
    y_all = Float64[]
    n_x_all = Float64[]
    n_y_all = Float64[]
    w_all = Float64[]
    w_n_all = Float64[]

    for curve in billiard.boundary
        if typeof(curve) <: AbsRealCurve
            L = curve.length
            N = round(Int, k*L*b/(2*pi))
            t, dt = sampler(N)
            ds = L*dt
            x, y = curve.r(t)
            nx, ny = curve.n(t)
            rn = (x .* nx .+ y .* ny)
            w = ds
            w_n = (w.*rn)./(2.0*k.^2) 
            append!(x_all, x)
            append!(y_all, y)
            append!(n_x_all, nx)
            append!(n_y_all, ny)
            append!(w_all, w)
            append!(w_n_all, w_n)
            
        end
    end
    return DMPoints(x_all, y_all,n_x_all,n_y_all, w_all, w_n_all)
end

function construct_matrices(solver::DecompositionMethod, basis::AbsBasis, pts::DMPoints, k)#{T} where T<:Real
    x, y, nx, ny, w, w_n = pts.x, pts.y, pts.n_x, pts.n_y, pts.w, pts.w_n
    M =  length(x)
    N = basis.dim
    B = Array{Float64}(undef,M,N)  #basis matrix
    for i in 1:N
        B[:,i] = basis_fun(basis, i, k, x, y)
    end 
    T = (w .* B) #reused later
    F = B' * T #boundary norm matrix
    #reuse B
    for i in 1:N #@inbounds Threads.@threads 
        dx, dy = grad_fun(basis, i, k, x, y)
        u = (dx.*nx) .+ (dy.*ny)
        B[:,i] = u
    end
    #B is now normal derivative matrix (u function)
    T = (w_n .* B) #apply integration weights
    G = B' * T #norm matrix
    return F, G    
end

function solve(solver::DecompositionMethod,basis::AbsBasis, pts::DMPoints, k)
    F, G = construct_matrices(solver, basis, pts, k)
    mu = generalized_eigvals(Symmetric(F),Symmetric(G);eps=solver.eps)
    lam0 = mu[end]
    t = 1.0/lam0
    return  t
end

function solve_vect(solver::DecompositionMethod,basis::AbsBasis, pts::DMPoints, k)
    F, G = construct_matrices(solver, basis, pts, k)
    mu, Z, C = generalized_eigen(Symmetric(F),Symmetric(G);eps=solver.eps)
    x = Z[:,end]
    x = C*x #transform into original basis 
    lam0 = mu[end]
    t = 1.0/lam0
    return  t, x./sqrt(lam0)
end


