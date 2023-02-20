include("../../abstracttypes.jl")
include("../decompositions.jl")
include("../matrixconstructors.jl")
#include("../samplers.jl")
#include("../src/billiards/triangle.jl")
using LinearAlgebra, StaticArrays

struct DecompositionMethod{T} <: SweepSolver where T<:Real
    dim_scaling_factor::T
    pts_scaling_factor::T
    eps::T
end

DecompositionMethod(dim_scaling_factor, pts_scaling_factor) = DecompositionMethod(dim_scaling_factor, pts_scaling_factor, 10.0*eps())

struct BoundaryPointsDM{T} <: AbsPoints where {T<:Real}
    xy::Vector{SVector{2,T}}
    normal::Vector{SVector{2,T}} #normal vectors in points
    w::Vector{T} # tension weights
    w_n::Vector{T} #normalization weights
end

function evaluate_points(solver::DecompositionMethod, billiard::AbsBilliard, sampler::Function, k)
    b = solver.pts_scaling_factor
    type = typeof(solver.pts_scaling_factor)
    xy_all = Vector{SVector{2,type}}()
    normal_all = Vector{SVector{2,type}}()
    w_all = Vector{type}()
    w_n_all = Vector{type}()

    for crv in billiard.boundary
        if typeof(crv) <: AbsRealCurve
            L = crv.length
            N = round(Int, k*L*b/(2*pi))
            t, dt = sampler(N)
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



function construct_matrices(solver::DecompositionMethod, basis::AbsBasis, pts::BoundaryPointsDM, k)#{T} where T<:Real
    #type = eltype(pts.w)
    #=
    B = basis_matrix(basis, k, pts.xy)
    =#
    B, dX, dY = basis_and_gradient_matrices(basis, k, pts.xy)
    T = (pts.w .* B) #reused later
    F = B' * T #boundary norm matrix
    nx = getindex.(pts.normal,1)
    ny = getindex.(pts.normal,2)
    #println(nx,ny)
    #println(size(dX))
    #inplace modifications
    dX = nx .* dX 
    dY = ny .* dY
    #reuse B
    B = dX .+ dY
    #B is now normal derivative matrix (u function)
    T = (pts.w_n .* B) #apply integration weights
    G = B' * T #norm matrix
    return F, G    
end

function solve(solver::DecompositionMethod,basis::AbsBasis, pts::BoundaryPointsDM, k)
    F, G = construct_matrices(solver, basis, pts, k)
    mu = generalized_eigvals(Symmetric(F),Symmetric(G);eps=solver.eps)
    lam0 = mu[end]
    t = 1.0/lam0
    return  t
end

function solve_vect(solver::DecompositionMethod,basis::AbsBasis, pts::BoundaryPointsDM, k)
    F, G = construct_matrices(solver, basis, pts, k)
    mu, Z, C = generalized_eigen(Symmetric(F),Symmetric(G);eps=solver.eps)
    x = Z[:,end]
    x = C*x #transform into original basis 
    lam0 = mu[end]
    t = 1.0/lam0
    return  t, x./sqrt(lam0)
end


