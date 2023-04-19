using LinearAlgebra, StaticArrays, TimerOutputs

struct BoundaryIntegralMethod{T} <: SweepSolver where {T<:Real}
    dim_scaling_factor::T
    pts_scaling_factor::Vector{T}
    sampler::Vector
    eps::T
    min_dim::Int64
    min_pts::Int64
end

#=
function DecompositionMethod(dim_scaling_factor::T, pts_scaling_factor::Union{T,Vector{T}}; min_dim = 100, min_pts = 500) where T<:Real 
    d = dim_scaling_factor
    bs = typeof(pts_scaling_factor) == T ? [pts_scaling_factor] : pts_scaling_factor
    sampler = [GaussLegendreNodes()]
return DecompositionMethod(d, bs, sampler, eps(T), min_dim, min_pts)
end

function DecompositionMethod(dim_scaling_factor::T, pts_scaling_factor::Union{T,Vector{T}}, samplers::Vector{AbsSampler}; min_dim = 100, min_pts = 500) where {T<:Real} 
    d = dim_scaling_factor
    bs = typeof(pts_scaling_factor) == T ? [pts_scaling_factor] : pts_scaling_factor
    return DecompositionMethod(d, bs, samplers, eps(T), min_dim, min_pts)
end
=#

struct BoundaryPointsBIM{T} <: AbsPoints where {T<:Real}
    xy::Vector{SVector{2,T}}
    normal::Vector{SVector{2,T}} #normal vectors in points
    kappa::Vector{T}
    ds::Vector{T}
    cosphi::Matrix{T} # tension weights
    tau::Matrix{T} #normalization weights
end

function evaluate_points(solver::BoundaryIntegralMethod, billiard::Bi, k) where {Bi<:AbsBilliard}
    bs, samplers = adjust_scaling_and_samplers(solver, billiard)
    curves = billiard.fundamental_boundary
    type = eltype(solver.pts_scaling_factor)
    
    xy_all = Vector{SVector{2,type}}()
    normal_all = Vector{SVector{2,type}}()
    kappa_all::Vector{T}()
    w_all = Vector{type}()
   
    for i in eachindex(curves)
        crv = curves[i]
        if typeof(crv) <: AbsRealCurve
            L = crv.length
            N = max(solver.min_pts,round(Int, k*L*bs[i]/(2*pi)))
            sampler = samplers[i]
            t, dt = sample_points(sampler,N)
            ds = L*dt #modify this
            xy = curve(crv,t)
            normal = normal_vec(crv,t)
            kappa = curvature(curve,t)
            append!(xy_all, xy)
            append!(normal_all, normal)
            append!(kappa_all, kappa)
            append!(w_all, ds)
        end
    end

    N = length(xy_all)
    cosphi = zeros(type,(N,N))
    tau = similar(cosphi)

    for j in 1:N
        for i in 1:N
            pt = xy_all[i] .- xy_all[j]
            tt = hypot(pt[1],pt[2])
            tau[i,j] = tt
            cosphi[i,j] = dot(normal_all[i], pt) / tt 
        end
    end

    return BoundaryPointsBIM{type}(xy_all,normal_all,kappa_all,w_all,cosphi,tau)
end
#=
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
=#
function construct_matrices(solver::BoundaryIntegralMethod, basis::Ba, pts::BoundaryPointsBIM, k) where {Ba<:AbsBasis}
    #basis and gradient matrices
    
    symmetries = basis.symmetries
    if ~isnothing(symmetries)
        norm = (length(symmetries)+1.0)

    end
    xy = pts.xy
    ds = pts.ds
    N = basis.dim
    #basis matrix
    #B = basis_matrix(basis, k, xy)
    Q = -2.0 * dk_matrix(basis, k, xy)
    #type = eltype(G)
    #F = zeros(type,(N,N))
    #Fk = similar(F)
    T = (w .* B) #reused later
    mul!(F,B',T) #boundary norm matrix
    #reuse B
    #B = dk_matrix(basis,k, xy)
    mul!(Fk,B',T) #B is now derivative matrix
    #symmetrize matrix
    Fk = Fk + Fk' 
    return F, Fk    
    
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
