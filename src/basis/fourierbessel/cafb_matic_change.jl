using StaticArrays
using Bessels

# new helpers for the bessels
Jv(nu, r) = Bessels.besselj(nu, r)

# Modified k derivative of bessel function
function Jvp(α::T, i::T, r::T) where {T<:Real}
    let
    j_minus = Jv(α*i-one(T),r)
    j_plus = Jv(α*i+one(T),r)
    return 0.5*(j_minus - j_plus)
    end
end

"""
Calculates all the radii for the given vector of points `pts`

# Arguments 
- `pts::Vector{SVector{T, 2}}` The points vector for which the radii will be calculated

# Returns
- `Vector{T}` The vector of all the radii of points
"""
function rs(pts::Vector{SVector{T, 2}}) where {T<:Real}
    return [sqrt(pt[i][1]^2 + pt[i][2]^2) for i in 1:length(pts)]
end

# Helper function to compute Bessel function values for all points
"""
Helper function that calculates all the Bessel functions for the given points. The Bessel function are indexed by `nu=α*i` where `α` is the `pi/corner_angle` and `i` the index of the basis function in the FB expansion

# Arguments
- `pts::Vector{SVector{T, 2}}` The points vector for which the Bessel functions will be calculated
- `α::T` The value of `α=pi/corner_angle`
- `i::Int` The index of the basis function
- `k::T` The wavenumber where we calculate the value of the Bessel function

# Returns
- `Vector{T}` The vector of all the Bessel function values at the given points. This is a vector of radial components of the basis
"""
function js(pts::Vector{SVector{T, 2}}, α::T, i::Int, k::T) where {T<:Real}
    # Get all the radii for points
    rs_values = rs(pts)
    return [Jv(α * i, k * r) for r in rs_values]
end

# Helper function to compute the derivative of Bessel function values for all points
"""
Computes the derivatives of the Bessel functions for the given points. The derivatives are indexed by `nu=α*i` where `α` is the `pi/corner_angle` and `i` the index of the basis function in the FB expansion

# Arguments
- `pts::Vector{SVector{T, 2}}` The points vector for which the derivative of Bessel functions will be calculated
- `α::T` The value of `α=pi/corner_angle`
- `i::Int` The index of the derivative of the  basis function
- `k::T` The wavenumber where we calculate the value of the derivative of the Bessel function

# Returns
- `Vector{T}` The vector of all the derivative of Bessel function values at the given points. This is a vector of radial components of the derivative of the basis
"""
function djs(pts::Vector{SVector{T, 2}}, α::T, i::Int, k::T) where {T<:Real}
    # Get all the radii for points
    rs_values = rs(pts)
    return [Jvp(α * i, k, r) for r in rs_values]
end

# Helper function for the sin_phi term that circumevants the discontinutity
"""
Computes the angular part of the Bessel function while circumventing the discontinutity of the atan2 function. It calculates the `sin(α*i*phi)` part of the basis.

# Arguments
- `pt::SVector{T, 2}` The point vector for which the angular part of the Bessel function will be calculated
- `i::Int` The index of the basis function
- `α::T` The value of `α=pi/corner_angle`

# Returns
- `Vector{T}` The values of the angular part of the FB basis
"""
function sin_phi(pt::SVector{T, 2}, i::Integer, α::T) where {T<:Real}
    return imag((pt[1]/sqrt(pt[1]^2 + pt[2]^2) + im*pt[2]/sqrt(pt[1]^2 + pt[2]^2))^(α*i))
end

"""
Helper function for calculating the `cos(α*i*phi)` that are useful in gradient calculations for the basis

# Arguments
- `pt::SVector{T, 2}` The point vector for which the `cos(α*i*phi)` will be calculated
- `i::Int` The index of the basis function
- `α::T` The value of `α=pi/corner_angle`

# Returns
- `Vector{T}` The values of the `cos(α*i*phi)` for the points
"""
function cos_phi(pt::SVector{T, 2}, i::Integer, α::T) where {T<:Real}
    return real((pt[1]/sqrt(pt[1]^2 + pt[2]^2) + im*pt[2]/sqrt(pt[1]^2 + pt[2]^2))^(α*i))
end

# Helper functions for the angular components to the ca_fb 
"""
Computes the `sin(α*i*phi)` part of the ca_fb basis while circumventing the discontinutity of the atan2 function. It does this for all the points in the input by evaluating each one with the single point function call of the same name

# Arguments
- `pts::Vector{SVector{T, 2}}` The points vector for which the angular parts of the Bessel function will be calculated
- `i::Int` The index of the basis function
- `α::T` The value of `α=pi/corner_angle`

# Returns
- `Vector{T}` The values of the angular part of the ca_fb basis for all the points
"""
function sin_phis(pts::Vector{SVector{T, 2}}, i::Integer, α::T) where {T<:Real}
    sin_phis_values = Vector{T}(length(pts))
    Threads.@threads for (j, pt) in enumerate(pts)
        sin_phis_vales[j] =sin_phi(pt, i, α)
    end
    return sin_phis_values
end

"""
Computes the `cos(α*i*phi)`'s for all the points given. It does this for all the points in the input by evaluating each one with the single point function call of the same name

# Arguments
- `pts::Vector{SVector{T, 2}}` The points vector for which the `cos(α*i*phi)`'s will be calculated
- `i::Int` The index of the basis function
- `α::T` The value of `α=pi/corner_angle`

# Returns
- `Vector{T}` The values of the `cos(α*i*phi)` for all the points
"""
function cos_phis(pts::Vector{SVector{T, 2}}, i::Integer, α::T) where {T<:Real}
    cos_phis_values = Vector{T}(length(pts))
    Threads.@threads for (j, pt) in enumerate(pts)
        cos_phis_vales[j] =cos_phi(pt, i, α)
    end
    return cos_phis_values
end

# Modified k derivative of the ca_fb function that use sin_phi directly
ca_fb(α, k, pt, i) = Jv(α*i, k*r)*sin_phi(pt, i, α)
ca_fb_dk(α, k, pt, i) = r*Jv(α*i, k*r)*sin_phi(pt, i, α)

"""
Struct representing the CornerAdaptedFourierBessel basis.

 Fields:
- `cs::PolarCS`: The coordinate system of the basis.
- `dim::Int64`: The dimension of the basis.
- `α::T`: The value of `α=pi/corner_angle`.
- `i::Integer`: The index of the basis
- `symmetries::Union{Vector{Sy}, Nothing}`: The symmetries of the billiard geoemtry
"""
struct CornerAdaptedFourierBessel{T,Sy} <: AbsBasis where  {T<:Real,Sy<:Union{AbsSymmetry,Nothing}}
    cs::PolarCS{T}
    dim::Int64 #using concrete type
    α::T # corner angle
    i::Integer #order constant intege -> nu = α * i
    symmetries::Union{Vector{Sy},Nothing}
end

"""
Constructor for the CornerAdaptedFourierBessel struct.

# Arguments
- `dim::Int64`: The dimension of the basis.
- `α::T`: The value of `α=pi/corner_angle`.
- `origin::SVector{T, 2}`: The origin of the coordinate system.
- `rot_angle::Real`: The rotation angle of the coordinate system.

 # Returns
 - `CornerAdaptedFourierBessel{Float64,Nothing}`: A CornerAdaptedFourierBessel struct with the given parameters.
"""
function CornerAdaptedFourierBessel(dim, α, origin, rot_angle)
    cs = PolarCS(origin, rot_angle)
    nu = pi/corner_angle
    return CornerAdaptedFourierBessel{Float64,Nothing}(cs, dim, corner_angle, nu, nothing)
end

"""
Constructor for CornerAdaptedFourierBessel struct.

# Arguments
- `dim::Int64`: The dimension of the basis.
- `α::T`: The value of `α=pi/corner_angle`.
- `cs::CoordinateSystem`: The coordinate system of the basis.

 # Returns
 - `CornerAdaptedFourierBessel{Float64, Nothing}`: A CornerAdaptedFourierBessel struct with the given parameters.
"""
function CornerAdaptedFourierBessel(dim, α, cs::CoordinateSystem)
    nu = pi/corner_angle
    return CornerAdaptedFourierBessel{Float64,Nothing}(cs, dim, corner_angle, nu, nothing)
end

# The regular mono index one
"""
Corner adapted Fourier Bessel basis function.

# Arguments
- `basis::CornerAdaptedFourierBessel`: The CornerAdaptedFourierBessel struct. This determines the potential rotation for the coordinate system
- `i::Int`: The index of the basis function.
- `k::T`: The wave number.

# Returns
- `Vector{T}`: The values of the basis function for the given points.
"""
@inline function basis_fun(basis::CornerAdaptedFourierBessel{T}, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map, pts=pts, α=basis.α
        pts_rotated_cart = collect([pm(pt) for pt in pts])
        return collect(ca_fb(α, k, pt, i) for pt in pts_rotated_cart)
    end
end

# The multi index one that calls the mono index one
"""
Corner Adapted Fourier Bessel basis construction for all the indices in the dimension `dim` of the `CornerAdaptedFourierBessel` struct.

# Arguments
- `basis::CornerAdaptedFourierBessel`: The CornerAdaptedFourierBessel struct. This determines the potential rotation for the coordinate system.
- `indices::AbstractArray`: The indices of the basis functions.
- `k::T`: The wave number.
- `pts::AbstractArray`: The points for which the basis functions will be evaluated.

# Returns
- `Matrix{T}`: The values of the basis functions for all the indices for the given points, that is "`B[:, i] for i in indices`"
"""
@inline function basis_fun(basis::CornerAdaptedFourierBessel{T}, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    α = basis.α
    # Initialize matrix for each idx storage
    M =  length(pts)
    N = length(indices)
    B = zeros(T,M,N)
    Threads.@threads for idx in eachindex(indices)
        B[:,idx] .= basis_fun(basis, idx, k, pts) # For each mono index we apply the affine transformation
    end
    return B 
end

"""
Computes the derivative of the Corner Adapted Fourier Bessel basis for the given index `i`.

# Arguments
- `basis::CornerAdaptedFourierBessel`: The CornerAdaptedFourierBessel struct. This determines the potential rotation for the coordinate system.
- `i::Int`: The index of the basis function.
- `k::T`: The wave number.
- `pts::AbstractArray`: The points for which the basis functions will be evaluated.

# Returns
- `Vector{T}`: The values of the derivative of the basis function for the given points.
"""
@inline function dk_fun(basis::CornerAdaptedFourierBessel{T}, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    #translation of coordiante origin
    let pm = basis.cs.local_map, pts=pts, α = basis.α
        pts_rotated_cart = collect([pm(pt) for pt in pts])
        # Prepare storage for sin_phi values
        ca_fb_dks = Vector{T}(length(pts))
        Threads.@threads for (j, pt) in enumerate(pts_rotated_cart)
            ca_fb_dks[j] = ca_fb_dk(α, k, pt, i)
        end
        return ca_fb_dks
    end
end

"""
Computes the derivative of the Corner Adapted Bessel basis for `indices`. This uses the single index `dk_fun` call under the hood.

# Arguments
- `basis::CornerAdaptedFourierBessel`: The CornerAdaptedFourierBessel struct. This determines the potential rotation for the coordinate system.
- `indices::AbstractArray`: The indices of the basis functions.
- `k::T`: The wave number.
- `pts::AbstractArray`: The points for which the basis functions will be evaluated.

# Returns
- `Matrix{T}`: The values of the derivative of the basis functions for all the indices for the given points, that is "`dB[:, i] for i in indices`"
"""
@inline function dk_fun(basis::CornerAdaptedFourierBessel{T}, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    α = basis.α
    # Initialize matrix for storaage for each index a dk_fun
    M =  length(pts)
    N = length(indices)
    dB_dk = zeros(T,M,N)
    Threads.@threads for idx in eachindex(indices)
        dB_dk[:,idx] = dk_fun(basis, idx, k, pts) # For each mono index we apply the affine transformation
    end
end

"""
Computes the gradient of the Corner Adapted Fourier Bessel function at the given index `i`

# Arguments
- `basis::CornerAdaptedFourierBessel`: The CornerAdaptedFourierBessel struct. This determines the potential rotation for the coordinate system.
- `i::Int`: The index of the basis function.
- `k::T`: The wave number.
- `pts::AbstractArray`: The points for which the basis functions will be evaluated.

# Returns
- `Tuple{Vector{T}, Vector{T}}`: The values of the gradient of the basis function at the given index for the given points. The first `Vector{T}` is the X-component of the grad and the second one is the Y-component
"""
function gradient(basis::CornerAdaptedFourierBessel, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map, nu=basis.α * i
        # Get the affine transform of cartesian to cartesian
        pt_xy = collect(pm(pt) for pt in pts)
        rs_values = rs(pt_xy) # Calculate rs
        # Initialize arrays to hold the gradient components
        N = length(pts)
        dx = Vector{T}(undef, N)
        dy = Vector{T}(undef, N)
        # Compute the Bessel function values and their derivatives
        j_values = js(pt_xy, basis.α, i, k)
        dj_values = djs(pt_xy, basis.α, i, k)
        # Precompute sin_phi and cos_phi values
        sin_phi_values = sin_phis(pt_xy, i, basis.α)
        cos_phi_values = cos_phis(pt_xy, i, basis.α)

        # multithread the dx & dy
        Threads.@threads for j in 1:N
            X, Y = pt_xy[j][1], pt_xy[j][2]
            r = rs_values[j]
            j_value = j_values[j]
            dj_value = dj_values[j]
            sin_phi_value = sin_phi_values[j]
            cos_phi_value = cos_phi_values[j]
            
            # Calculate the gradient components
            dx[j] = dj_value * k * (X / r) * sin_phi_value - nu * j_value * cos_phi_value * (Y / r^2)
            dy[j] = dj_value * k * (Y / r) * sin_phi_value + nu * j_value * cos_phi_value * (X / r^2)
        end
    return dx, dy
    end
end

"""
Computes the gradient of the Corner Adapted Fourier Bessel function for a set of `indices`. This allows to compute the basis for all the indices

# Arguments
- `basis::CornerAdaptedFourierBessel`: The CornerAdaptedFourierBessel struct. This determines the potential rotation for the coordinate system.
- `indices::AbstractArray`: The indices of the basis functions.
- `k::T`: The wave number.
- `pts::AbstractArray`: The points for which the basis functions will be evaluated.

# Returns
- `Tuple{Matrix{T}, Matrix{T}}`: The values of the gradients of the basis functions for all the indices with the given points, that is "`dB_dx[:, i] for i in indices`" and "`dB_dy[:, i] for i in indices`"

"""
function gradient(basis::CornerAdaptedFourierBessel, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map, α=basis.α
        # Initialize matrices for the columns gradients for each index
        M = length(pts)
        N = length(indices)
        dB_dx = zeros(T,M,N)
        dB_dy = zeros(T,M,N)
        Threads.@threads for idx in eachindex(indices)
            grad = gradient(basis, idx, k, pts) # For each mono index we apply the affine transformation
            dB_dx[:,i] .= grad[1]
            dB_dy[:,i] .= grad[2]
        end
    return dB_dx, dB_dy
    end
end

"""
Computes the Corner Adapted Fourier Bessel basis and it's gradient for the given index `i` and points `pts`

# Arguments
- `basis::CornerAdaptedFourierBessel`: The CornerAdaptedFourierBessel struct. This determines the potential rotation for the coordinate system.
- `i::Int`: The index of the basis function.
- `k::T`: The wave number.
- `pts::AbstractArray`: The points for which the basis functions will be evaluated.

# Returns
- `Tuple{Vector{T}, Vector{T}, Vector{T}}`: The values of the basis function, its X-component of the gradient, and its Y-component of the gradient for the given index for the given points, that is `bf[i], dx[i], dy[i] for i in indices`"
"""
function basis_and_gradient(basis::CornerAdaptedFourierBessel, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    bf = basis_fun(basis, i, k, pts)
    dx, dy = gradient(basis, i, k, pts)
    return bf, dx, dy
end

"""
Computes the Corner Adapted Fourier Bessel basis and it's gradient for a `Vector` of indexes `indices` and points `pts`

# Arguments
- `basis::CornerAdaptedFourierBessel`: The CornerAdaptedFourierBessel struct. This determines the potential rotation for the coordinate system.
- `indices::AbstractArray`: The indices of the basis functions.
- `k::T`: The wave number.
- `pts::AbstractArray`: The points for which the basis functions will be evaluated.

# Returns
- `Tuple{Matrix{T}, Matrix{T}, Matrix{T}}`: The values of the basis functions, their X-component of the gradients, and their Y-component of the gradients for all the indices with the given points, that is "`B[:, i], dB_dx[:, i], dB_dy[:, i] for i in indices`"
"""
function basis_and_gradient(basis::CornerAdaptedFourierBessel, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    M =  length(pts)
    N = length(indices)
    B = zeros(T,M,N)
    dB_dx = zeros(T,M,N)
    dB_dy = zeros(T,M,N)
    Threads.@threads for idx in eachindex(indices)
        bf, dx, dy = basis_and_gradient(basis, idx, k, pts) # For each mono index we apply the affine transformation
        B[:,idx] .= bf
        dB_dx[:,idx] .= dx
        dB_dy[:,idx] .= dy
    end
    return B, dB_dx, dB_dy
end