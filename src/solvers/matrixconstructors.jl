#include("../abstracttypes.jl")

"""
Set small values in a matrix to zero in-place. This enables sparse matrices representation

# Description
This function iterates over the elements of the matrix `M` and sets any value with an absolute magnitude smaller than or equal to the specified tolerance `ϵ` to zero.

# Arguments
- `M`: The matrix to be filtered.
- `ϵ`: The threshold below which elements are set to zero. Default is the machine epsilon for the element type of `M`.

# Returns
- The filtered matrix `M`, with small values set to zero.
"""
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
"""
Construct the basis matrix for a given set of points and a basis and filters it.

# Description
This function constructs the basis matrix `B` for a given set of points `pts` and a specified basis `basis`. The matrix is filtered to remove insignificant values. The indices of the matrix are constructed as `1:dim`, where dim is the dimension field of the `basis`.

# Arguments
- `basis::Ba`: The basis to be used, of type `AbsBasis`.
- `k`: The wavenumber for which the basis matrix is constructed.
- `pts::Vector{SVector{2,T}}`: A vector of 2D points where the basis functions are evaluated.

# Returns
- The filtered basis matrix `B`.
"""
function basis_matrix(basis::Ba, k, pts::Vector{SVector{2,T}}) where {T<:Real, Ba<:AbsBasis}
    let dim = basis.dim
        B = basis_fun(basis,1:dim,k,pts)
        return filter_matrix!(B)
    end
end

#perhaps unnecesary
"""
Construct the basis matrix for specific indices of the basis functions.

# Description
This function constructs the basis matrix `B` for specific `indices` of the basis functions, given a set of points `pts` and a basis `basis`. The matrix is filtered to remove insignificant values. This is different from the version where the indices are gotten from the the `dim`.

# Arguments
- `basis::Ba`: The basis to be used, of type `AbsBasis`.
- `k`: The wavenumber for which the basis matrix is constructed.
- `pts::Vector{SVector{2,T}}`: A vector of 2D points where the basis functions are evaluated.
- `indices::AbstractArray`: Indices specifying which basis functions to use.

# Returns
- The filtered basis matrix `B`.
"""
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

#TODO #rework these
"""
Construct the basis matrix using x and y grids.

# Description
This function constructs the basis matrix `B` using x and y coordinate grids, evaluating the basis functions at each grid point. The resulting matrix is filtered to remove insignificant values. It uses the indexes internally from the dimension of the basis from `basis`.

# Arguments
- `basis::Ba`: The basis to be used, of type `AbsBasis`.
- `k`: The wavenumber for which the basis matrix is constructed.
- `x_grid::AbstractArray`: The grid of x-coordinates.
- `y_grid::AbstractArray`: The grid of y-coordinates.

# Returns
- The filtered basis matrix `B`.
"""
function basis_matrix(basis::Ba, k, x_grid::AbstractArray, y_grid::AbstractArray) where {T<:Real, Ba<:AbsBasis}
    let dim = basis.dim
        pts = collect(SVector(x,y) for y in y_grid for x in x_grid) 
        return basis_fun(basis,1:dim,k,pts)
    end
end

"""
Construct the basis matrix for specific indices (not from the `basis` dimension) using x and y grids.

# Description
This function constructs the basis matrix `B` for specific indices of the basis functions, using x and y coordinate grids. The matrix is filtered to remove insignificant values.

# Arguments
- `basis::Ba`: The basis to be used, of type `AbsBasis`.
- `k`: The wavenumber for which the basis matrix is constructed.
- `x_grid::AbstractArray`: The grid of x-coordinates.
- `y_grid::AbstractArray`: The grid of y-coordinates.
- `indices::AbstractArray`: Indices specifying which basis functions to use.

# Returns
- The filtered basis matrix `B`.
"""
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

# TODO Check if same as previous
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

"""
Construct the gradient matrices for a basis at given points. The indices are constructed internally from the dimension of the `basis`.

# Description
This function constructs the gradient matrices `dB_dx` and `dB_dy` for a basis `basis` at a set of points `pts`. The matrices represent the derivatives of the basis functions with respect to x and y. The resulting matrices are filtered to remove insignificant values.

# Arguments
- `basis::Ba`: The basis to be used, of type `AbsBasis`.
- `k`: The wavenumber for which the gradient matrices are constructed.
- `pts::Vector{SVector{2,T}}`: A vector of 2D points where the basis functions are evaluated.

# Returns
- The filtered gradient matrices `dB_dx` and `dB_dy`.
"""
function gradient_matrices(basis::Ba, k, pts::Vector{SVector{2,T}}) where {T<:Real, Ba<:AbsBasis}
    let dim = basis.dim
        dB_dx, dB_dy = gradient(basis,1:dim,k,pts)
        return filter_matrix!(dB_dx), filter_matrix!(dB_dy)
    end
end

"""
Construct the gradient matrices for specific indices at given points.

# Description
This function constructs the gradient matrices `dB_dx` and `dB_dy` for specific indices of the basis functions at a set of points `pts`. The matrices represent the derivatives of the basis functions with respect to x and y. The resulting matrices are filtered to remove insignificant values.

# Arguments
- `basis::Ba`: The basis to be used, of type `AbsBasis`.
- `k`: The wavenumber for which the gradient matrices are constructed.
- `pts::Vector{SVector{2,T}}`: A vector of 2D points where the basis functions are evaluated.
- `indices::AbstractArray`: Indices specifying which basis functions to use.

# Returns
- The filtered gradient matrices `dB_dx` and `dB_dy`.
"""
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

"""
Construct the basis and gradient matrices for a basis at given points. The indices are constructed internally from the dimension of the `basis`.

# Description
This function constructs both the basis matrix `B` and the gradient matrices `dB_dx` and `dB_dy` for a basis `basis` at a set of points `pts`. The matrices are filtered to remove insignificant values.

# Arguments
- `basis::Ba`: The basis to be used, of type `AbsBasis`.
- `k`: The wavenumber for which the matrices are constructed.
- `pts::Vector{SVector{2,T}}`: A vector of 2D points where the basis functions are evaluated.

# Returns
- The filtered basis matrix `B` and the filtered gradient matrices `dB_dx` and `dB_dy`.
"""
function basis_and_gradient_matrices(basis::Ba, k, pts::Vector{SVector{2,T}}) where {T<:Real, Ba<:AbsBasis}
    let dim = basis.dim
        B, dB_dx, dB_dy = basis_and_gradient(basis,1:dim,k,pts)
        return filter_matrix!(B), filter_matrix!(dB_dx), filter_matrix!(dB_dy)
    end
end

"""
Construct the basis and gradient matrices for specific indices at given points.

# Description
This function constructs both the basis matrix `B` and the gradient matrices `dB_dx` and `dB_dy` for specific indices of the basis functions at a set of points `pts`. The matrices are filtered to remove insignificant values.

# Arguments
- `basis::Ba`: The basis to be used, of type `AbsBasis`.
- `k`: The wavenumber for which the matrices are constructed.
- `pts::Vector{SVector{2,T}}`: A vector of 2D points where the basis functions are evaluated.
- `indices::AbstractArray`: Indices specifying which basis functions to use.

# Returns
- The filtered basis matrix `B`, and the filtered gradient matrices `dB_dx` and `dB_dy`.
"""
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

"""
Construct the derivative matrix at the wavenumber `k` for the particular basis. The indices are constructed internally from the dimension of the `basis`.

# Description
This function constructs the derivative matrix `dB_dk` for a basis `basis` with respect to the wavenumber `k` at a set of points `pts`. The matrix is filtered to remove insignificant values.

# Arguments
- `basis::Ba`: The basis to be used, of type `AbsBasis`.
- `k`: The wavenumber for which the derivative matrix is constructed.
- `pts::Vector{SVector{2,T}}`: A vector of 2D points where the basis functions are evaluated.

# Returns
- The filtered derivative matrix `dB_dk`.
"""
function dk_matrix(basis::Ba, k, pts::Vector{SVector{2,T}}) where {T<:Real, Ba<:AbsBasis}
    let dim = basis.dim
        dB_dk = dk_fun(basis,1:dim,k,pts)
        return filter_matrix!(dB_dk)
    end
end

"""
Construct the derivative matrix with respect to the wavenumber for specific indices.

# Description
This function constructs the derivative matrix `dB_dk` for specific indices of the basis functions with respect to the wavenumber `k` at a set of points `pts`. The matrix is filtered to remove insignificant values.

# Arguments
- `basis::Ba`: The basis to be used, of type `AbsBasis`.
- `k`: The wavenumber for which the derivative matrix is constructed.
- `pts::Vector{SVector{2,T}}`: A vector of 2D points where the basis functions are evaluated.
- `indices::AbstractArray`: Indices specifying which basis functions to use.

# Returns
- The filtered derivative matrix dB_dk.
"""
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
