#include("../../abstracttypes.jl")
#include("../decompositions.jl")
#include("../matrixconstructors.jl")
#include("../samplers.jl")
#include("../src/billiards/triangle.jl")
#include("../../utils/benchmarkutils.jl")
using LinearAlgebra, StaticArrays, TimerOutputs

struct DecompositionMethod{T,F} <: SweepSolver where {T<:Real,F<:Function}
    dim_scaling_factor::T
    pts_scaling_factor::T
    sampler::F
    eps::T
end

DecompositionMethod(dim_scaling_factor, pts_scaling_factor) = DecompositionMethod(dim_scaling_factor, pts_scaling_factor, gauss_legendre_nodes, eps(typeof(dim_scaling_factor)))

struct BoundaryPointsDM{T} <: AbsPoints where {T<:Real}
    xy::Vector{SVector{2,T}}
    normal::Vector{SVector{2,T}} #normal vectors in points
    w::Vector{T} # tension weights
    w_n::Vector{T} #normalization weights
end

function evaluate_points(solver::DecompositionMethod, billiard::Bi, k) where {Bi<:AbsBilliard}
    sampler = solver.sampler
    b = solver.pts_scaling_factor
    type = typeof(solver.pts_scaling_factor)
    xy_all = Vector{SVector{2,type}}()
    normal_all = Vector{SVector{2,type}}()
    w_all = Vector{type}()
    w_n_all = Vector{type}()

    for crv in billiard.fundamental_boundary
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

function construct_matrices_benchmark(solver::DecompositionMethod, basis::Ba, pts::BoundaryPointsDM, k) where {Ba<:AbsBasis}
    to = TimerOutput()
    w = pts.w
    w_n = pts.w_n
    symmetries=basis.symmetries
    if ~isnothing(symmetries)
        norm = (length(symmetries)+1.0)
        w = w.*norm
        w_n = w_n.*norm
    end
    #basis and gradient matrices
    @timeit to "basis_and_gradient_matrices" B, dX, dY = basis_and_gradient_matrices(basis, k, pts.xy)
    N = basis.dim
    type = eltype(B)
    F = zeros(type,(N,N))
    G = similar(F)
    
    @timeit to "F construction" begin 
        @timeit to "weights" T = (w .* B) #reused later
        @timeit to "product" mul!(F,B',T) #boundary norm matrix
    end

    @timeit to "G construction" begin 
        @timeit to "normal derivative" nx = getindex.(pts.normal,1)
        @timeit to "normal derivative" ny = getindex.(pts.normal,2)
        #inplace modifications
        @timeit to "normal derivative" dX = nx .* dX 
        @timeit to "normal derivative" dY = ny .* dY
        #reuse B
        @timeit to "normal derivative" B = dX .+ dY
    
    #B is now normal derivative matrix (u function)
        @timeit to "weights" T = (w_n .* B) #apply integration weights
        @timeit to "product" mul!(G,B',T)#norm matrix
    end
    print_timer(to)
    return F, G    
end

function construct_matrices(solver::DecompositionMethod, basis::Ba, pts::BoundaryPointsDM, k) where {Ba<:AbsBasis}
    #basis and gradient matrices
    w = pts.w
    w_n = pts.w_n
    symmetries=basis.symmetries
    if ~isnothing(symmetries)
        norm = (length(symmetries)+1.0)
        w = w.*norm
        w_n = w_n.*norm
    end
    B, dX, dY = basis_and_gradient_matrices(basis, k, pts.xy)
    type = eltype(B)
    #allocate matrices
    N = basis.dim
    F = zeros(type,(N,N))
    G = similar(F)
    #apply weights
    T = (w .* B) #reused later
    mul!(F,B',T) #boundary norm matrix
    nx = getindex.(pts.normal,1)
    ny = getindex.(pts.normal,2)
    #inplace modifications
    dX = nx .* dX 
    dY = ny .* dY
    #reuse B
    B = dX .+ dY
    #B is now normal derivative matrix (u function)
    T = (w_n .* B) #apply integration weights
    mul!(G,B',T)#norm matrix
    return F, G    
    
end

function solve(solver::DecompositionMethod,basis::Ba, pts::BoundaryPointsDM, k) where {Ba<:AbsBasis}
    F, G = construct_matrices(solver, basis, pts, k)
    mu = generalized_eigvals(Symmetric(F),Symmetric(G);eps=solver.eps)
    lam0 = mu[end]
    t = 1.0/lam0
    return  t
end

function solve(solver::DecompositionMethod,F,G)
    #F, G = construct_matrices(solver, basis, pts, k)
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
