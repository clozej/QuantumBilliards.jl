include("../abstracttypes.jl")
include("decompositions.jl")
#include("../src/billiards/triangle.jl")
using LinearAlgebra

struct ScalingMethod <: AcceleratedSolver 
    dim_scaling_factor::Float64
    pts_scaling_factor::Float64
    eps::Float64
end

ScalingMethod(dim_scaling_factor, pts_scaling_factor) = ScalingMethod(dim_scaling_factor, pts_scaling_factor, 10.0*eps())

struct SMPoints <: AbsPoints
    x::Vector{Float64}
    y::Vector{Float64}
    w::Vector{Float64}
end

function evaluate_points(solver::ScalingMethod, billiard::AbsBilliard, sampler::Function, k)
    b = solver.pts_scaling_factor
    x_all = Float64[]
    y_all = Float64[]
    w_all = Float64[]

    for curve in billiard.boundary
        if typeof(curve) <: AbsRealCurve
            L = curve.length
            N = round(Int, k*L*b/(2*pi))
            t, dt = sampler(N)
            ds = L*dt
            x, y = curve.r(t)
            nx, ny = curve.n(t)
            rn = (x .* nx .+ y .* ny)
            w = ds ./ rn
            append!(x_all, x)
            append!(y_all, y)
            append!(w_all, w)
        end
    end
    return SMPoints(x_all, y_all, w_all)
end

#generalize for other types
function construct_matrices(solver::ScalingMethod, basis::AbsBasis, pts::SMPoints, k)#{T} where T<:Real
    x, y, w = pts.x, pts.y, pts.w
    M =  length(x)
    N = basis.dim
    B = Array{Float64}(undef,M,N)  #basis matrix
    @inbounds Threads.@threads for i in 1:N
        B[:,i] = basis_fun(basis, i, k, x, y)
    end 
    T = (w .* B) #reused later
    F = B' * T #boundary norm matrix
    #reuse B
    @inbounds Threads.@threads for i in 1:N
        B[:,i] = dk_fun(basis, i, k, x, y)
    end
    Fk = B' * T #B is now derivative matrix
    #symmetrize matrix
    Fk = Fk + Fk'
    return F, Fk    
end

function sm_results(mu,k)
    ks = k .- 2 ./mu #.+ 2/k ./(mu.^2) 
    ten = 2.0 .*(2.0 ./ mu).^2
    p = sortperm(ks)
    return ks[p], ten[p]
end

function solve(solver::ScalingMethod,basis::AbsBasis, pts::SMPoints, k, dk)
    F, Fk = construct_matrices(solver, basis, pts, k)
    mu = generalized_eigvals(Symmetric(F),Symmetric(Fk);eps=solver.eps)
    ks, ten = sm_results(mu,k)
    idx = abs.(ks.-k) .< dk
    return ks[idx], ten[idx]
end

function solve(solver::ScalingMethod,F,Fk, k, dk)
    #F, Fk = construct_matrices(solver, basis, pts, k)
    mu = generalized_eigvals(Symmetric(F),Symmetric(Fk);eps=solver.eps)
    ks, ten = sm_results(mu,k)
    idx = abs.(ks.-k) .< dk
    return ks[idx], ten[idx]
end


#=
function construct_matrices(solver::ScalingMethod, basis::AbsBasis, k)#{T} where T<:Real
    x, y, w = pts.x, pts.y, pts.w
    M =  length(x)
    N = basis.dim
    B = Array{Float64}(undef,M,N)  #basis matrix
    for i in 1:N
        B[:,i] = basis_fun(basis, i, k, x, y)
    end 
    T = (w .* B) #reused later
    F = B' * T #boundary norm matrix
    #reuse B
    for i in 1:N
        B[:,i] = dk_fun(basis, i, k, x, y)
    end
    Fk = B' * T #B is now derivative matrix
    #symmetrize matrix
    Fk = Fk + Fk'
    return F, Fk    
end
=#

#=
#non allocating version
function construct_matrices!(out,basis::AbsBasis, pts::SMPoints, k) where T<:Real
    N = basis.dim
      #basis matrix
    for i in 1:N
        out[:,i] = basis_fun(basis, i, k, pts.x, pts.y)
    end 
end
=#