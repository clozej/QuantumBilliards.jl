#include("../abstracttypes.jl")

function filter_matrix!(M; ϵ = eps(eltype(M)))
    type = eltype(M)
    #max = maximum(abs,M)
    k = 1
    for t in eachindex(M)
        if abs.(M[t]) <= ϵ
            M[t] = zero(type)
            k += 1
        end
    end
    return(M)
end

#this will be usefull for basis sets containing several functions (plane and evanscent waves etc.)
function basis_matrix(basis::Ba, k, pts::Vector{SVector{2,T}}; parallel_matrix = true) where {T<:Real, Ba<:AbsBasis}
    let dim = basis.dim
        B = basis_fun(basis,1:dim,k,pts; parallel_matrix)
        return filter_matrix!(B)
    end
end

function gradient_matrices(basis::Ba, k, pts::Vector{SVector{2,T}}; parallel_matrix = true) where {T<:Real, Ba<:AbsBasis}
    let dim = basis.dim
        dB_dx, dB_dy = gradient(basis,1:dim,k,pts; parallel_matrix)
        return filter_matrix!(dB_dx), filter_matrix!(dB_dy)
    end
end

function basis_and_gradient_matrices(basis::Ba, k, pts::Vector{SVector{2,T}}; parallel_matrix = true) where {T<:Real, Ba<:AbsBasis}
    let dim = basis.dim
        B, dB_dx, dB_dy = basis_and_gradient(basis,1:dim,k,pts; parallel_matrix)
        return filter_matrix!(B), filter_matrix!(dB_dx), filter_matrix!(dB_dy)
    end
end


function dk_matrix(basis::Ba, k, pts::Vector{SVector{2,T}}; parallel_matrix = true) where {T<:Real, Ba<:AbsBasis}
    let dim = basis.dim
        dB_dk = dk_fun(basis,1:dim,k,pts; parallel_matrix)
        return filter_matrix!(dB_dk)
    end
end

#=
#for now unnecesary
function basis_matrix(basis::Ba, k, pts::Vector{SVector{2,T}}, indices::AbstractArray) where {T<:Real, Ba<:AbsBasis}
    let N = basis.dim
        M =  length(pts)
        B = zeros(T,M,N) 
        B1 = basis_fun(basis,indices,k,pts)
        #println(size(B))
        #println(size(B1))
        for i in indices
            B[:,i] .= B1[:,i]
        end
        return filter_matrix!(B)
    end
end

function basis_and_gradient_matrices(basis::Ba, k, pts::Vector{SVector{2,T}}, indices::AbstractArray) where {T<:Real, Ba<:AbsBasis}
    let N = basis.dim
        M =  length(pts)
        dX = zeros(T,M,N)
        dY = zeros(T,M,N)
        B1 = zeros(T,M,N)
        B, dB_dx, dB_dy = basis_and_gradient(basis,indices,k,pts)
        filter_matrix!(B)
        filter_matrix!(dB_dx)
        filter_matrix!(dB_dy)
        #println(size(B))
        #println(size(B1))
        for i in indices
            B1[:,i] .= B[:,i]
            dX[:,i] .= dB_dx[:,i]
            dY[:,i] .= dB_dy[:,i]
        end
        return B1, dX, dY
    end
end



function gradient_matrices(basis::Ba, k, pts::Vector{SVector{2,T}}, indices::AbstractArray) where {T<:Real, Ba<:AbsBasis}
    let N = basis.dim
        M =  length(pts)
        dX = zeros(T,M,N)
        dY = zeros(T,M,N)  
        dB_dx, dB_dy = gradient(basis,indices,k,pts)
        filter_matrix!(dB_dx)
        filter_matrix!(dB_dy)
        #println(size(B))
        #println(size(B1))
        for i in indices
            dX[:,i] .= dB_dx[:,i]
            dY[:,i] .= dB_dy[:,i]
        end
        return dX, dY
    end
end

function dk_matrix(basis::Ba, k, pts::Vector{SVector{2,T}}, indices::AbstractArray) where {T<:Real, Ba<:AbsBasis}
    let N = basis.dim
        M =  length(pts)
        dB1 = zeros(T,M,N) 
        dB_dk = dk_fun(basis,indices,k,pts)
        filter_matrix!(dB_dk)
        #println(size(B))
        #println(size(B1))
        for i in indices
            dB1[:,i] .= dB_dk[:,i]
        end
        return dB1
    end
end
=#
#=
#rework these
function basis_matrix(basis::Ba, k, x_grid::AbstractArray, y_grid::AbstractArray) where {T<:Real, Ba<:AbsBasis}
    let dim = basis.dim
        pts = collect(SVector(x,y) for y in y_grid for x in x_grid) 
        return basis_fun(basis,1:dim,k,pts)
    end
end

function basis_matrix(basis::Ba, k, x_grid::AbstractArray, y_grid::AbstractArray, indices::AbstractArray) where {T<:Real, Ba<:AbsBasis}
    let N = basis.dim
        M =  length(x_grid)*length(y_grid)
        B = zeros(eltype(x_grid),M,N) 
        pts = collect(SVector(x,y) for y in y_grid for x in x_grid)
        B1 = basis_fun(basis,indices,k,pts)
        #println(size(B))
        #println(size(B1))
        for i in indices
            B[:,i] .= B1[:,i]
        end
        return B
    end
end

function basis_matrix(basis::Ba, k, x_grid::AbstractArray, y_grid::AbstractArray, indices::AbstractArray) where {T<:Real, Ba<:AbsBasis}
    let N = basis.dim
        M =  length(x_grid)*length(y_grid)
        B = zeros(eltype(x_grid),M,N)
        pts = collect(SVector(x,y) for y in y_grid for x in x_grid) 
        B1 = basis_fun(basis,indices,k,pts)
        #println(size(B))
        #println(size(B1))
        for i in indices
            B[:,i] .= B1[:,i]
        end
        return B
    end
end
=#