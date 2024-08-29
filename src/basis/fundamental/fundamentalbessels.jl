

using Bessels
using CoordinateTransformations, Rotations, StaticArrays

#aliases so one is able to change bessel library easily
Y0(x) = Bessels.bessely0(x)
dY0(x) = -Bessels.bessely1(x)

struct FundamentalBessel{T,Sy} <: AbsBasis where  {T<:Real, Sy<:Union{AbsSymmetry,Nothing}}
    #cs::PolarCS{T} #not fully implemented
    dim::Int64 #using concrete type
    dilation::T
    positions::Vector{SVector{2,T}}
    symmetries::Union{Vector{Sy},Nothing}
end

"""
This function creates a default `FundamentalBessel` basis with a dimension of 0, no dilation, a single position at the origin, and no symmetries.
"""
function FundamentalBessel()
    return FundamentalBessel{Float64,Nothing}(0, 0.0, [SVector(0.0,0.0)], nothing)
end

"""
This function creates a `FundamentalBessel` basis with a dimension of 0, no dilation, a single position at the origin, and the specified symmetries.

# Arguments
- `symmetries::Vector{Sy}`: A vector of symmetries to be applied to the basis functions, where `Sy` is a subtype of `AbsSymmetry`.

# Returns
- A `FundamentalBessel` object initialized with the specified symmetries.
"""
function FundamentalBessel(symmetries::Vector{Sy}) where  {Sy<:AbsSymmetry}
    return FundamentalBessel(0, 0.0, [SVector(0.0,0.0)], symmetries)
end

"""
This function resizes the `FundamentalBessel` basis to a specified dimension based on the given wavenumber `k` and the geometry of the billiard domain. The resized basis includes appropriately positioned points based on either the symmetries or the boundary of the billiard domain.

# Logic
- If the basis has symmetries, positions are determined using the `LinearNodes` sampler.
- If the basis has no symmetries, positions are determined using the `FourierNodes` sampler based on the full boundary of the billiard.
- The wavenumber `k` is adjusted to ensure it is above a minimum value of 25.0. This is to ensure we have matrices of enough size
- The function then creates a new `FundamentalBessel` basis with the resized dimension and the computed positions.

# Arguments
- `basis::Ba`: The `FundamentalBessel` basis to be resized.
- `billiard::Bi`: The billiard domain, where `Bi` is a subtype of `AbsBilliard`.
- `dim::Int`: The new dimension for the basis.
- `k`: The wavenumber.

# Returns
- A `FundamentalBessel` object resized to the specified dimension with updated positions based on the billiard domain and wavenumber.
"""
function resize_basis(basis::Ba, billiard::Bi, dim::Int, k) where {Ba<:FundamentalBessel,Bi<:AbsBilliard}
    k_min = 25.0
    k_basis = max(k_min,k)
    if ~isnothing(basis.symmetries)
        sampler = LinearNodes()
        pos = dilated_boundary_points(billiard, sampler, dim, k_basis)
    else
        boundary = billiard.full_boundary
        crv_lengths = [crv.length for crv in boundary]
        sampler = FourierNodes(nothing,crv_lengths)
        pos = dilated_boundary_points(billiard, sampler, dim, k_basis)
    end
    return FundamentalBessel{typeof(k),eltype(basis.symmetries)}(length(pos),2*pi/k_basis,pos,basis.symmetries)
end

"""
This function computes the Bessel basis function for a specified index `i` in the `FundamentalBessel` basis at a given wavenumber `k`. The function evaluates the Bessel function of the second kind `Y0` at each point in the input array, adjusted for the position of the basis function.

# Logic
- The points `pts` are translated by subtracting the position of the basis function indexed by `i`.
- The radial distance `hypot(pt[1], pt[2])` is computed for each translated point.
- The Bessel function of the second kind `Y0` is evaluated at each `k*r`.
- If symmetries are present, the symmetric images of the points are computed, and their contributions to the basis function are added.

# Returns
- A vector containing the values of the Bessel basis function evaluated at the input points.

# Arguments
- `basis`: The `FundamentalBessel` basis.
- `i`: The index of the basis function.
- `k`: The wavenumber.
- `pts`: An array of points where the basis function is to be evaluated.
"""
@inline function basis_fun(basis::FundamentalBessel{T}, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let symmetries = basis.symmetries
        pts_loc = [pt - basis.positions[i] for pt in pts]
        arg = [k*hypot(pt[1], pt[2]) for pt in pts_loc]
        b = Y0.(arg)
        if ~isnothing(symmetries)
            for sym in symmetries
                pts_sym = sym.sym_map.(pts)
                pts_loc = [pt - basis.positions[i] for pt in pts_sym]
                arg = [k*hypot(pt[1], pt[2]) for pt in pts_loc]
                b .+= sym.parity.*Y0.(arg)
            end
            b = b./(length(symmetries)+1.0)
        end #norm::T = one(T)/sqrt(basis.dim)
        return b
    end
end

#try using strided for this
"""
This function computes the Bessel basis functions for multiple specified indices in the `FundamentalBessel` basis at a given wavenumber `k`. The function evaluates the Bessel function of the second kind `Y0` at each point in the input array for each index, adjusted for the position of the corresponding basis function. 

# Logic
- For each index in `indices`:
- The points `pts` are translated by subtracting the position of the corresponding basis function.
- The radial distance `hypot(pt[1], pt[2])` is computed for each translated point.
- The Bessel function of the second kind `Y0` is evaluated at each `k*r`.
- If symmetries are present, the symmetric images of the points are computed, and their contributions to the basis function are added.

# Returns
- A matrix with dimensions `(number of points, number of indices)`, where each column contains the values of the Bessel basis function for a specific index evaluated at the input points.

# Arguments
- `basis`: The `FundamentalBessel` basis.
- `indices`: An array of indices for which the basis functions are to be computed.
- `k`: The wavenumber.
- `pts`: An array of points where the basis functions are to be evaluated.
"""
@inline function basis_fun(basis::FundamentalBessel{T}, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let symmetries = basis.symmetries
        M =  length(pts)
        N = length(indices)
        B = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            pts_loc = [pt - basis.positions[i] for pt in pts]
            arg = [k*hypot(pt[1], pt[2]) for pt in pts_loc]
            B[:,i] .= Y0.(arg)
            if ~isnothing(symmetries)
                for sym in symmetries
                    pts_sym = sym.sym_map.(pts)
                    pts_loc = [pt - basis.positions[i] for pt in pts_sym]
                    arg = [k*hypot(pt[1], pt[2]) for pt in pts_loc]
                    B[:,i] .+= sym.parity.*Y0.(arg)
                end
            end
        end
        if ~isnothing(symmetries)
            return B./(length(symmetries)+1.0)
        else
            return B
        end 
    end
end

"""
This function computes both the Bessel basis function and its gradient with respect to the Cartesian coordinates `x` and `y` for a specified index `i` in the `FundamentalBessel` basis. The function evaluates these quantities at the provided points, taking into account the position of the basis function.

# Logic
- The points `pts` are translated by subtracting the position of the basis function indexed by `i`.
- The radial distance `r = hypot(pt[1], pt[2])` is computed for each translated point.
- The Bessel function of the second kind `Y0(k*r)` is evaluated for each `k*r`.
- The gradient components `dx` and `dy` are computed as the partial derivatives of the basis function with respect to `x` and `y`, respectively.
- If symmetries are present, the symmetric images of the points are computed, and their contributions to the basis function and gradient are added.

# Returns
- `bf`: A vector containing the values of the Bessel basis function evaluated at the input points.
- `dx`: A vector containing the `x` component of the gradient of the basis function.
- `dy`: A vector containing the `y` component of the gradient of the basis function.

# Arguments
- `basis`: The `FundamentalBessel` basis.
- `i`: The index of the basis function.
- `k`: The wavenumber.
- `pts`: An array of points where the basis function and its gradient are to be evaluated.
"""
function basis_and_gradient(basis::FundamentalBessel{T}, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let symmetries = basis.symmetries
        pts_loc = [pt - basis.positions[i] for pt in pts]
        r = [hypot(pt[1], pt[2]) for pt in pts_loc]
        x = getindex.(pts_loc, 1)
        y = getindex.(pts_loc, 2) 
        arg = k.*r
        #norm::T = one(T)/sqrt(basis.dim)
        bf = Y0.(arg)
        dx = @. k*x/r*dY0(arg)
        dy = @. k*y/r*dY0(arg)
        if ~isnothing(symmetries)
            n = (length(symmetries)+1.0)
            for sym in symmetries
                pts_sym = sym.sym_map.(pts)
                pts_loc = [pt - basis.positions[i] for pt in pts_sym]
                r = [hypot(pt[1], pt[2]) for pt in pts_loc]
                x = getindex.(pts_loc, 1)
                y = getindex.(pts_loc, 2) 
                arg = k.*r
                bf .+= sym.parity.*Y0.(arg)
                dx .+= @. sym.parity*k*x/r*dY0(arg)
                dy .+= @. sym.parity*k*y/r*dY0(arg)
            end
            bf = bf./n
            dx = dx./n
            dy = dy./n
        end
        return bf, dx, dy
    end
end

"""
This function computes both the Bessel basis functions and their gradients with respect to the Cartesian coordinates `x` and `y` for multiple specified indices in the `FundamentalBessel` basis. The function evaluates these quantities at the provided points, considering the effect of the positions of the basis functions.

# Logic
- For each index in `indices`:
- The points `pts` are translated by subtracting the position of the corresponding basis function.
- The radial distance `r = hypot(pt[1], pt[2])` is computed for each translated point.
- The Bessel function of the second kind `Y0(k*r)` is evaluated for each `k*r`.
- The gradient components `dB_dx` and `dB_dy` are computed as the partial derivatives of the basis function with respect to `x` and `y`, respectively.
- If symmetries are present, the symmetric images of the points are computed, and their contributions to the basis function and gradient are added.

# Returns
- `B`: A matrix with dimensions `(number of points, number of indices)`, where each column contains the values of the Bessel basis function for a specific index evaluated at the input points.
- `dB_dx`: A matrix with dimensions `(number of points, number of indices)`, where each column contains the `x` component of the gradient of the basis function for a specific index.
- `dB_dy`: A matrix with dimensions `(number of points, number of indices)`, where each column contains the `y` component of the gradient of the basis function for a specific index.

# Arguments
- `basis`: The `FundamentalBessel` basis.
- `indices`: An array of indices for which the basis functions and gradients are to be computed.
- `k`: The wavenumber.
- `pts`: An array of points where the basis functions and gradients are to be evaluated.
"""
function basis_and_gradient(basis::FundamentalBessel{T}, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let symmetries = basis.symmetries
        M =  length(pts)
        N = length(indices)
        B = zeros(T,M,N)
        dB_dx = zeros(T,M,N)
        dB_dy = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            pts_loc = [pt - basis.positions[i] for pt in pts]
            r = [hypot(pt[1], pt[2]) for pt in pts_loc]
            x = getindex.(pts_loc, 1)
            y = getindex.(pts_loc, 2) 
            arg = k.*r
            #println(size(s))
            B[:,i] .= Y0.(arg)
            dB_dx[:,i] .= @. k*x/r*dY0(arg)
            dB_dy[:,i] .= @. k*y/r*dY0(arg)
            if ~isnothing(symmetries)
                for sym in symmetries
                    pts_sym = sym.sym_map.(pts)
                    pts_loc = [pt - basis.positions[i] for pt in pts_sym]
                    r = [hypot(pt[1], pt[2]) for pt in pts_loc]
                    x = getindex.(pts_loc, 1)
                    y = getindex.(pts_loc, 2) 
                    arg = k.*r
                    B[:,i] .+= sym.parity.*Y0.(arg)
                    dB_dx[:,i] .+= @. sym.parity*k*x/r*dY0(arg)
                    dB_dy[:,i] .+= @. sym.parity*k*y/r*dY0(arg)
                end
            end
        end
        if ~isnothing(symmetries)
            n = (length(symmetries)+1.0)
            return B./n, dB_dx./n, dB_dy./n
        else
        #println(size(s))
            return B, dB_dx, dB_dy
        end
    end
end

"""
This function computes the derivative of the Bessel basis function with respect to the wavenumber `k` for a specified index `i` in the `FundamentalBessel` basis. The derivative is evaluated at each point in the input array. If symmetries are present, the contributions from symmetric points are included.

# Logic
- The points `pts` are translated by subtracting the position of the basis function indexed by `i`.
- The radial distance `r = hypot(pt[1], pt[2])` is computed for each translated point.
- The derivative of the Bessel function of the second kind, `dY0(k*r)`, is evaluated for each `k*r`.
- The result is multiplied by the radial distance `r` to obtain the derivative of the basis function with respect to `k`.
- If symmetries are present, the symmetric images of the points are computed, and their contributions to the derivative are added.

# Returns
- A vector containing the derivative of the Bessel basis function with respect to `k` for the specified index `i`.

# Arguments
- `basis`: The `FundamentalBessel` basis.
- `i`: The index of the basis function.
- `k`: The wavenumber.
- `pts`: An array of points where the derivative is to be evaluated.
"""
@inline function dk_fun(basis::FundamentalBessel{T}, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let symmetries = basis.symmetries
        pts_loc = [pt - basis.positions[i] for pt in pts]
        r = [hypot(pt[1], pt[2]) for pt in pts_loc]
        dk = @. r*dY0(k*r)
        if ~isnothing(symmetries)
            n = (length(symmetries)+1.0)
            for sym in symmetries
                pts_sym = sym.sym_map.(pts)
                pts_loc = [pt - basis.positions[i] for pt in pts_sym]
                r = [hypot(pt[1], pt[2]) for pt in pts_loc]
                dk .+= @. sym.parity*r*dY0(k*r)
            end
            dk = dk./n
        end
        #norm::T = one(T)/sqrt(basis.dim)
        return dk
    end
end
    
"""
This function computes the derivatives of the Bessel basis functions with respect to the wavenumber `k` for multiple specified indices in the `FundamentalBessel` basis. The derivatives are evaluated at each point in the input array, considering the effect of the position of each basis function.

# Logic
- For each index in `indices`:
- The points `pts` are translated by subtracting the position of the corresponding basis function.
- The radial distance `r = hypot(pt[1], pt[2])` is computed for each translated point.
- The derivative of the Bessel function of the second kind, `dY0(k*r)`, is evaluated for each `k*r`.
- The result is multiplied by the radial distance `r` to obtain the derivative of the basis function with respect to `k`.
- If symmetries are present, the symmetric images of the points are computed, and their contributions to the derivative are added.

# Returns
- A matrix with dimensions `(number of points, number of indices)`, where each column contains the derivative of the Bessel basis function with respect to `k` for a specific index evaluated at the input points.

# Arguments
- `basis`: The `FundamentalBessel` basis.
- `indices`: An array of indices for which the derivatives are to be computed.
- `k`: The wavenumber.
- `pts`: An array of points where the derivatives are to be evaluated.
"""
@inline function dk_fun(basis::FundamentalBessel{T}, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let symmetries = basis.symmetries
        M =  length(pts)
        N = length(indices)
        dB_dk = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            pts_loc = [pt - basis.positions[i] for pt in pts]
            r = [hypot(pt[1], pt[2]) for pt in pts_loc]
            dB_dk[:,i] .= @. r*dY0(k*r)
            if ~isnothing(symmetries)
                for sym in symmetries
                    pts_sym = sym.sym_map.(pts)
                    pts_loc = [pt - basis.positions[i] for pt in pts_sym]
                    r = [hypot(pt[1], pt[2]) for pt in pts_loc]
                    dB_dk[:,i] .+= @. sym.parity*r*dY0(k*r)
                end
            end
        end
        if ~isnothing(symmetries)
            n = (length(symmetries)+1.0)
            return dB_dk./n
        else
            return dB_dk
        end
    end
end