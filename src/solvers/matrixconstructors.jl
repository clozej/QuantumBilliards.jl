include("../abstracttypes.jl")

function basis_matrix(basis::AbsBasis,k, x::Vector{T}, y::Vector{T}) where T<:Number
    M =  length(x)
    N = basis.dim
    B = zeros(T,M,N)  #basis matrix
    @inbounds Threads.@threads for i in axes(B,2)
        @inbounds @simd for j in axes(B,1)
            B[j,i] = basis_fun(basis, i, k, x[j], y[j])
        end
    end 
    return B
end

function basis_matrix(basis::AbsBasis,k, x::Vector{T}, y::Vector{T}, dim_ind) where T<:Number
    M =  length(x)
    N = basis.dim
    B = zeros(T,M,N)  #basis matrix
    @inbounds Threads.@threads for i in dim_ind
        @inbounds @simd for j in axes(B,1)
            B[j,i] = basis_fun( basis, i, k, x[j], y[j])
        end
    end 
    return B
end

function basis_matrix(basis::AbsBasis,k, x::Vector{T}, y::Vector{T}, dim_ind, pts_ind) where T<:Number
    M =  length(x)
    N = basis.dim
    B = zeros(T,M,N)  #basis matrix
    @inbounds Threads.@threads for i in dim_ind
        @inbounds @simd for j in pts_ind
            B[j,i] = basis_fun( basis, i, k, x[j], y[j])
        end
    end 
    return B
end

function U_matrix(basis::AbsBasis,k, x::Vector{T}, y::Vector{T},nx::Vector{T}, ny::Vector{T}) where T<:Number
    M =  length(x)
    N = basis.dim
    U = zeros(T,M,N)  #basis matrix
    @inbounds Threads.@threads for i in axes(U,2) #
        @inbounds @simd for j in axes(U,1)
            dx, dy = grad_fun(basis, i, k, x[j], y[j])
            U[j,i] = dx*nx[j] + dy*ny[j]
        end
    end 
    return U
end