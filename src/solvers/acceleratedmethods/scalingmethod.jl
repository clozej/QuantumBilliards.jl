#include("../../abstracttypes.jl")
#include("../decompositions.jl")
#include("../matrixconstructors.jl")
#include("../src/billiards/triangle.jl")
#include("../../utils/benchmarkutils.jl")
using LinearAlgebra, StaticArrays
using TimerOutputs
abstract type AbsScalingMethod <: AcceleratedSolver 
end
struct ScalingMethodA{T} <: AbsScalingMethod where T<:Real
    dim_scaling_factor::T
    pts_scaling_factor::T
    eps::T
end

ScalingMethodA(dim_scaling_factor, pts_scaling_factor) = ScalingMethodA(dim_scaling_factor, pts_scaling_factor, eps(typeof(dim_scaling_factor)))


struct ScalingMethodB{T} <: AbsScalingMethod where T<:Real
    dim_scaling_factor::T
    pts_scaling_factor::T
    eps::T
end

ScalingMethodB(dim_scaling_factor, pts_scaling_factor) = ScalingMethodB(dim_scaling_factor, pts_scaling_factor, eps(typeof(dim_scaling_factor)))


struct BoundaryPointsSM{T} <: AbsPoints where {T<:Real}
    xy::Vector{SVector{2,T}}
    w::Vector{T}
end

function evaluate_points(solver::AbsScalingMethod, billiard::Bi, sampler::Function, k) where {Bi<:AbsBilliard}
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
function construct_matrices_benchmark(solver::ScalingMethodA, basis::Ba, pts::BoundaryPointsSM, k) where {Ba<:AbsBasis}
    to = TimerOutput()
    #type = eltype(pts.w)
    xy, w = pts.xy, pts.w
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

function construct_matrices(solver::ScalingMethodA, basis::Ba, pts::BoundaryPointsSM, k) where {Ba<:AbsBasis}
    xy = pts.xy
    w = pts.w
    basis=basis
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
    basis=basis
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

#missing solve_vect and solve_vectors


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