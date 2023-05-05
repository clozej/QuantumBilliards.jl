using LinearAlgebra
#=
function GESVDVALS(A, B; eps=1e-14)
    M = reduce(vcat, [A, B]) #concatenate columns
    Q, R =  qr(M)
    _ , sv_r, Vt_r = svd(R)
    mask = sv_r .> eps
    #println(mask)
    V1t = Vt_r[ : , mask]

    return svdvals((A * V1t), (B * V1t))
end
=#

function generalized_eigen(A,B;eps=1e-15)
    d, S = eigen(A)
    idx = (d/maximum(d)) .> eps
    q = 1.0 ./ sqrt.(d[idx])
    #println(length(q))
    C = q' .* S[:,idx] 
    D = B * C
    #println(size(D))
    E = Symmetric(C' * D)
    #println(size(E))
    mu, Z = eigen(E) #check eigenvectors
    return mu, Z, C
end

function generalized_eigvals(A,B;eps=1e-15)
    d, S = eigen(A)
    idx = (d/maximum(d)) .> eps
    q = 1.0 ./ sqrt.(d[idx])
    #println(length(q))
    C = q' .* S[:,idx] 
    #println(size(C))
    D = B * C
    #println(size(D))
    E = Symmetric(C' * D)
    #println(size(E))
    mu = eigvals(E)
    return mu
end

directsum(A,B) = [A zeros(size(A,1), size(B,2)); zeros(size(B,1), size(A,2)) B]

function adjust_scaling_and_samplers(solver::AbsSolver, billiard::AbsBilliard)
    bs = solver.pts_scaling_factor
    samplers = solver.sampler
    default = samplers[1]
    n_curves = length(billiard.fundamental_boundary)
    #=
    for crv in billiard.fundamental_boundary
        if typeof(crv) <: AbsRealCurve
            n_curves += 1
        end
    end
    =#
    b_min = minimum(bs)
    while length(bs)<n_curves
        push!(bs, b_min)
    end
    
    while length(samplers)<n_curves
        push!(samplers, default)
    end
    return bs, samplers
end
