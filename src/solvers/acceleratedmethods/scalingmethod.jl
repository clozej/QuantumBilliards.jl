include("../../abstracttypes.jl")
include("../decompositions.jl")
include("../matrixconstructors.jl")
#include("../src/billiards/triangle.jl")
using LinearAlgebra, StaticArrays

struct ScalingMethod{T} <: AcceleratedSolver where T<:Number
    dim_scaling_factor::T
    pts_scaling_factor::T
    eps::T
end

ScalingMethod(dim_scaling_factor, pts_scaling_factor) = ScalingMethod(dim_scaling_factor, pts_scaling_factor, 10.0*eps(typeof(dim_scaling_factor)))

struct BoundaryPointsSM{T} <: AbsPoints where {T<:Number}
    xy::Vector{SVector{2,T}}
    w::Vector{T}
end

function evaluate_points(solver::ScalingMethod, billiard::AbsBilliard, sampler::Function, k)
    b = solver.pts_scaling_factor
    type = typeof(solver.pts_scaling_factor)
    xy_all = Vector{SVector{2,type}}()
    w_all = Vector{type}()

    for crv in billiard.boundary
        if typeof(crv) <: AbsRealCurve
            L = crv.length
            N = round(Int, k*L*b/(2*pi))
            t, dt = sampler(N)
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
function construct_matrices(solver::ScalingMethod, basis::AbsBasis, pts::BoundaryPointsSM, k)#{T} where T<:Real
    type = eltype(pts.w)
    xy, w = pts.xy, pts.w
    M =  length(xy)
    N = basis.dim
    #basis matrix
    B = basis_matrix(basis, k, xy)
    T = (w .* B) #reused later
    F = B' * T #boundary norm matrix
    #reuse B
    B = dk_matrix(basis ,k, xy)
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

function solve(solver::ScalingMethod,basis::AbsBasis, pts::BoundaryPointsSM, k, dk)
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