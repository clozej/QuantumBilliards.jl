
#include("../../abstracttypes.jl")
#include("../../utils/coordinatesystems.jl")
#include("../../utils/gridutils.jl")
#using BenchmarkTools
#using SpecialFunctions
using Bessels
using CoordinateTransformations, Rotations, StaticArrays

#Jv(nu, r) = SpecialFunctions.besselj(nu,r)
"""
Symbolic definition of the Bessel function of order `nu` with radius `r`
"""
Jv(nu, r) = Bessels.besselj(nu, r)

"""
This function calculates the cylindrical wave expansion term using the Bessel function of the first kind `Jv(nu, k*r)` and the sine of the angular component.

# Logic
- The Bessel function of the first kind `Jv(nu, k*r)` is computed for the given order `nu`, wavenumber `k`, and radial distance `r`.
- The result is multiplied by the sine of the product of `nu` and the azimuthal angle `phi`. This is the standard difference reccurence relation for the derivative

# Returns
- The value of the cylindrical wave expansion term for the given parameters.

# Arguments
- `nu`: The order of the Bessel function.
- `k`: The wavenumber.
- `r`: The radial coordinate.
- `phi`: The azimuthal angle.
"""
ca_fb(nu,k,r,phi) = Jv(nu, k*r)*sin(nu*phi)
function Jvp(nu, r::T) where {T<:Real}
    let
    j_minus = Jv(nu-one(T),r)
    j_plus = Jv(nu+one(T),r)
    return 0.5*( j_minus - j_plus)
    end
end

"""
Compute the derivative of the cylindrical wave expansion term with respect to the wavenumber `k`.

# Logic
- The derivative of the Bessel function `Jv(nu, r)` with respect to its argument `r` is computed using the `Jvp` function.
- The result is multiplied by the radial coordinate `r` (implicit differentiation) and the sine of the product of `nu` and the angle `phi`.

# Returns
- The derivative of the cylindrical wave expansion term with respect to `k`.

# Arguments
- `nu`: The order of the Bessel function.
- `k`: The wavenumber.
- `r`: The radial coordinate.
- `phi`: The angle.
"""
ca_fb_dk(nu,k,r,phi) = r * Jvp(nu, k*r) * sin(nu * phi)

#=
function _bessel_diff_formula(v, z, n, L, phase)
    # from AMS55.
    # L(v, z) = J(v, z), Y(v, z), H1(v, z), H2(v, z), phase = -1
    # L(v, z) = I(v, z) or exp(v*pi*i)K(v, z), phase = 1
    # For K, you can pull out the exp((v-k)*pi*i) into the caller
    p = 1.0
    s = L(v-n, z)
    for i in 1:n
        p = phase * (p * (n-i+1)) / i   # = choose(k, i)
        s .+= p*L(v-n + i*2, z)
    end
    return s ./ (2.0^n)
end
=#

struct CornerAdaptedFourierBessel{T,Sy} <: AbsBasis where  {T<:Real,Sy<:Union{AbsSymmetry,Nothing}}
    cs::PolarCS{T}
    dim::Int64 #using concrete type
    corner_angle::T
    nu::T #order constant, order=nu*i
    symmetries::Union{Vector{Sy},Nothing}
end

"""
This function constructs a `CornerAdaptedFourierBessel` basis, which is specifically designed to handle boundary value problems in domains with sharp corners. The basis functions are adapted to the corner angle, ensuring accurate representation of the solution in these regions.
For further reading please check Betcke's paper: Reviving the Method of Particular Solutions, Timo Betcke & Lloyd N. Trefethen

# Logic
- The function takes as input the dimension `dim`, the corner angle `corner_angle`, and the origin and rotation angle for defining a polar coordinate system.
- The order constant `nu` is computed as `pi/corner_angle`.
- A `PolarCS` coordinate system is created using the provided origin and rotation angle (depending on the implementation version chosen).
- The function returns a `CornerAdaptedFourierBessel` object initialized with these parameters.

# Returns
- A `CornerAdaptedFourierBessel` object, configured with the specified dimension, corner angle, and coordinate system.
"""
function CornerAdaptedFourierBessel(dim, corner_angle, origin, rot_angle)
    cs = PolarCS(origin, rot_angle)
    nu = pi/corner_angle
    return CornerAdaptedFourierBessel{Float64,Nothing}(cs, dim, corner_angle, nu, nothing)
end

function CornerAdaptedFourierBessel(dim, corner_angle, cs::CoordinateSystem)
    nu = pi/corner_angle
    return CornerAdaptedFourierBessel{Float64,Nothing}(cs, dim, corner_angle, nu, nothing)
end

"""
Convert a `CornerAdaptedFourierBessel` basis to use `Float32` precision.
# Returns
- A `CornerAdaptedFourierBessel` object with all relevant parameters converted to `Float32`.
"""
toFloat32(basis::CornerAdaptedFourierBessel) = CornerAdaptedFourierBessel(basis.dim, Float32(basis.corner_angle), Float32.(basis.cs.origin), Float32(basis.cs.rot_angle))

"""
This function resizes the `CornerAdaptedFourierBessel` basis to a new dimension, if necessary. It checks whether the current dimension matches the desired dimension and returns the resized basis if they differ.
- If the dimensions match, the original basis is returned.
- If the dimensions differ, a new `CornerAdaptedFourierBessel` object is created with the new dimension and the existing corner angle and coordinate system.

# Returns
- A `CornerAdaptedFourierBessel` object with the updated dimension, or the original basis if no resizing is needed.
"""
function resize_basis(basis::CornerAdaptedFourierBessel, billiard::Bi, dim::Int, k) where {Bi<:AbsBilliard}
    if basis.dim == dim
        return basis
    else
        return CornerAdaptedFourierBessel(dim, basis.corner_angle, basis.cs)
    end
end

#evaluation functions
#=
function local_coords(basis::CornerAdaptedFourierBessel{T}, pt::SVector{2,T}) where {T<:Number}
     #new_pt::SVector{2,T} = cartesian_to_polar(f(pt))
    return cartesian_to_polar(basis.cs.local_map(pt))
end  
=#

"""
This function computes the basis function for a specified index `i` in the `CornerAdaptedFourierBessel` basis at a given wavenumber `k`. The function maps the input points to local polar coordinates and evaluates the corresponding Fourier-Bessel function.

# Logic
- The points `pts` are mapped to local polar coordinates using the `local_map` of the basis's coordinate system.
- The polar coordinates (`r`, `phi`) are computed for each point.
- The Bessel function `ca_fb(nu*i, k, r, phi)` is evaluated for the given index `i` and wavenumber `k` for those points.

# Returns
- A vector containing the values of the basis function evaluated at the input points.

# Arguments
- `basis`: The `CornerAdaptedFourierBessel` basis.
- `i`: The index of the basis function.
- `k`: The wavenumber.
- `pts`: An array of points where the basis function is to be evaluated.
"""
@inline function basis_fun(basis::CornerAdaptedFourierBessel{T}, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map, nu=basis.nu, pts=pts
        pt_pol = (cartesian_to_polar(pm(pt)) for pt in pts)
        #norm::T = one(T)/sqrt(basis.dim)
        return collect(ca_fb(nu*i, k, pt[1], pt[2]) for pt in pt_pol)
    end
end

"""
This function computes the basis functions for multiple specified indices in the `CornerAdaptedFourierBessel` basis at a given wavenumber `k`. The function maps the input points to local polar coordinates and evaluates the corresponding Fourier-Bessel functions in parallel (using Threads).

# Logic
- The points `pts` are mapped to local polar coordinates using the `local_map` of the basis's coordinate system.
- The polar coordinates (`r`, `phi`) are computed for each point.
- For each index in `indices`, the Bessel function `ca_fb(nu*i, k, r, phi)` is evaluated in parallel using multithreading.

# Returns
- A matrix with dimensions `(number of points, number of indices)`, where each column contains the values of the basis function for a specific index evaluated at the input points.

# Arguments
- `basis`: The `CornerAdaptedFourierBessel` basis.
- `indices`: An array of indices for which the basis functions are to be computed.
- `k`: The wavenumber.
- `pts`: An array of points where the basis functions are to be evaluated.
"""
@inline function basis_fun(basis::CornerAdaptedFourierBessel{T}, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map, nu=basis.nu, pts=pts
        pt_pol = (cartesian_to_polar(pm(pt)) for pt in pts)
        #norm::T = one(T)/sqrt(basis.dim)
        M =  length(pts)
        N = length(indices)
        B = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            B[:,i] .= (ca_fb(nu*i, k, pt[1], pt[2]) for pt in pt_pol)
        end
        return B 
    end
end

"""
This function computes the derivative of the basis function with respect to the wavenumber `k` for a specified index `i` in the `CornerAdaptedFourierBessel` basis. The derivative is evaluated at the provided points in the local polar coordinates.

# Logic
- The points `pts` are mapped to local polar coordinates using the `local_map` of the basis's coordinate system.
- The radial coordinate `r` and the angular coordinate `phi` are extracted from the polar coordinates.
- The derivative of the Bessel function `Jv(nu*i, k*r)` with respect to its argument `k*r` is computed using the `Jvp` function.
- The result is multiplied by the radial coordinate `r` and the sine of the angular component `nu*i*phi` to obtain the derivative of the basis function with respect to `k`.

# Returns
- A vector containing the derivative of the basis function with respect to the wavenumber `k` for the specified index `i`.

# Arguments
- `basis`: The `CornerAdaptedFourierBessel` basis.
- `i`: The index of the basis function.
- `k`: The wavenumber.
- `pts`: An array of points where the derivative is to be evaluated.
"""
@inline function dk_fun(basis::CornerAdaptedFourierBessel{T}, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    #translation of coordiante origin
    let pm = basis.cs.local_map, nu=basis.nu, pts=pts
        pt_pol = [cartesian_to_polar(pm(pt)) for pt in pts]
        #norm::T = one(T)/sqrt(basis.dim)
        r = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2) 
        dj = @. Jvp(nu*i, k*r)
        s = @. sin(nu*i*phi)
        dk = @. r*dj*s
        return dk
    end
end
    
"""
This function computes the derivatives of the basis functions with respect to the wavenumber `k` for multiple specified indices in the `CornerAdaptedFourierBessel` basis. The derivatives are evaluated at the provided points in the local polar coordinates.

# Logic
- The points `pts` are mapped to local polar coordinates using the `local_map` of the basis's coordinate system.
- The radial coordinate `r` and the angular coordinate `phi` are extracted from the polar coordinates.
- For each index in `indices`, the derivative of the Bessel function `Jv(nu*i, k*r)` with respect to its argument `k*r` is computed using the `Jvp` function.
- The result for each index is multiplied by the radial coordinate `r` and the sine of the angular component `nu*i*phi`.
- The function returns a matrix where each column corresponds to the derivative of a basis function with respect to `k` for one of the indices, and each row corresponds to a point in `pts`.

# Returns
- A matrix with dimensions `(number of points, number of indices)`, where each column contains the derivative of the basis function with respect to `k` for a specific index evaluated at the input points.

# Arguments
- `basis`: The `CornerAdaptedFourierBessel` basis.
- `indices`: An array of indices for which the derivatives are to be computed.
- `k`: The wavenumber.
- `pts`: An array of points where the derivatives are to be evaluated.
"""
@inline function dk_fun(basis::CornerAdaptedFourierBessel{T}, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map, nu=basis.nu, pts=pts
        pt_pol = [cartesian_to_polar(pm(pt)) for pt in pts]
        #norm::T = one(T)/sqrt(basis.dim)
        r = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2)
        M =  length(pts)
        N = length(indices)
        dB_dk = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            dj = @. Jvp(nu*i, k*r)
            s = @. sin(nu*i*phi)
            dB_dk[:,i] .= @. r*dj*s
        end
        return dB_dk
    end
end

"""
This function computes the gradient of the basis function with respect to the Cartesian coordinates `x` and `y` for a specified index `i` in the `CornerAdaptedFourierBessel` basis. The gradient is evaluated at the provided points in local Cartesian coordinates.

# Logic
- The points `pts` are mapped to local Cartesian coordinates using the `local_map` of the basis's coordinate system.
- These local Cartesian coordinates are converted to polar coordinates (`r`, `phi`).
- The Bessel function `Jv(nu*i, k*r)` and its derivative with respect to `k*r` are computed.
- The gradient components `dx` and `dy` are calculated using the chain rule, which involves the Bessel function, its derivative, and the trigonometric functions of `phi`.
- The function returns the gradient as two vectors: one for the `x` component (`dx`) and one for the `y` component (`dy`).

# Returns
- Two vectors: `dx` and `dy`, representing the gradient of the basis function with respect to the Cartesian coordinates `x` and `y`.

# Arguments
- `basis`: The `CornerAdaptedFourierBessel` basis.
- `i`: The index of the basis function.
- `k`: The wavenumber.
- `pts`: An array of points where the gradient is to be evaluated.
"""
function gradient(basis::CornerAdaptedFourierBessel, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map, nu=basis.nu, pts=pts
        pt_xy = collect(pm(pt) for pt in pts)
        pt_pol = collect(cartesian_to_polar(pt) for pt in pt_xy) #local cartesian coords
        #norm::T = one(T)/sqrt(basis.dim)
        r = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2)
        X = getindex.(pt_xy,1)
        Y = getindex.(pt_xy,2)
        j = Jv.(nu*i, k*r)
        #println(size(j))
        dj = Jvp.(nu*i, k*r) 
        #println(size(dj))
        s = @. sin(nu*i*phi) 
        c = @. cos(nu*i*phi) 
        #println(size(s))
        dx = @. (dj*k*(X/r)*s-nu*i*j*c*Y/(r^2))
        dy = @. (dj*k*(Y/r)*s+nu*i*j*c*X/(r^2))
    return dx, dy
    end
end

"""
This function computes the gradients of the basis functions with respect to the Cartesian coordinates `x` and `y` for multiple specified indices in the `CornerAdaptedFourierBessel` basis. The gradients are evaluated at the provided points in local Cartesian coordinates.

# Logic
- The points `pts` are mapped to local Cartesian coordinates using the `local_map` of the basis's coordinate system.
- These local Cartesian coordinates are converted to polar coordinates (`r`, `phi`).
- For each index in `indices`, the Bessel function `Jv(nu*i, k*r)` and its derivative with respect to `k*r` are computed.
- The gradient components `dx` and `dy` are calculated for each index using the chain rule, which involves the Bessel function, its derivative, and the trigonometric functions of `phi`.
- The function returns two matrices: one for the `x` component (`dB_dx`) and one for the `y` component (`dB_dy`), where each column corresponds to an index and each row corresponds to a point.

# Returns
- Two matrices: `dB_dx` and `dB_dy`, as a `Tuple` representing the gradients of the basis functions with respect to the Cartesian coordinates `x` and `y`.

# Arguments
- `basis`: The `CornerAdaptedFourierBessel` basis.
- `indices`: An array of indices for which the gradients are to be computed.
- `k`: The wavenumber.
- `pts`: An array of points where the gradients are to be evaluated.
"""
function gradient(basis::CornerAdaptedFourierBessel, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map, nu=basis.nu, pts=pts
        #local cartesian coords
        pt_xy = collect(pm(pt) for pt in pts)
        pt_pol = collect(cartesian_to_polar(pt) for pt in pt_xy)
        #norm::T = one(T)/sqrt(basis.dim)
        r = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2)
        X = getindex.(pt_xy,1)
        Y = getindex.(pt_xy,2)
        M = length(pts)
        N = length(indices)
        dB_dx = zeros(T,M,N)
        dB_dy = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            j = Jv.(nu*i, k*r)
            #println(size(j))
            dj = Jvp.(nu*i, k*r) 
            #println(size(dj))
            s = @. sin(nu*i*phi) 
            c = @. cos(nu*i*phi) 
            #println(size(s))
            dB_dx[:,i] .= @. (dj*k*(X/r)*s-nu*i*j*c*Y/(r^2))
            dB_dy[:,i] .= @. (dj*k*(Y/r)*s+nu*i*j*c*X/(r^2))
        end
        #println(size(s))
    return dB_dx, dB_dy
    end
end

"""
This function computes both the basis function and its gradient with respect to the Cartesian coordinates `x` and `y` for a specified index `i` in the `CornerAdaptedFourierBessel` basis. The function evaluates these quantities at the provided points. This is a composite of `basis_fun` and `gradient`.

# Logic
- The points `pts` are mapped to local Cartesian coordinates using the `local_map` of the basis's coordinate system.
- These local Cartesian coordinates are converted to polar coordinates (`r`, `phi`).
- The Bessel function `Jv(nu*i, k*r)` and its derivative with respect to `k*r` are computed.
- The basis function `bf` is calculated as the product of the Bessel function and the sine of the angular coordinate `nu*i*phi`.
- The gradient components `dx` and `dy` are calculated using the chain rule, which involves the Bessel function, its derivative, and the trigonometric functions of `phi`.
- The function returns the basis function and the gradient as three separate vectors: `bf`, `dx`, and `dy`.

# Returns
- `bf`: A vector containing the values of the basis function evaluated at the input points.
- `dx`: A vector containing the `x` component of the gradient of the basis function.
- `dy`: A vector containing the `y` component of the gradient of the basis function.

# Arguments
- `basis`: The `CornerAdaptedFourierBessel` basis.
- `i`: The index of the basis function.
- `k`: The wavenumber.
- `pts`: An array of points where the basis function and its gradient are to be evaluated.
"""
function basis_and_gradient(basis::CornerAdaptedFourierBessel, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map, nu=basis.nu, pts=pts
        pt_xy = collect(pm(pt) for pt in pts)
        pt_pol = collect(cartesian_to_polar(pt) for pt in pt_xy) #local cartesian coords
        #norm::T = one(T)/sqrt(basis.dim)
        r = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2)
        X = getindex.(pt_xy,1)
        Y = getindex.(pt_xy,2)
        j = Jv.(nu*i, k*r)
        #println(size(j))
        dj = Jvp.(nu*i, k*r) 
        #println(size(dj))
        s = @. sin(nu*i*phi) 
        c = @. cos(nu*i*phi) 

        #println(size(s))
        bf = @. j*s
        dx = @. (dj*k*(X/r)*s-nu*i*j*c*Y/(r^2))
        dy = @. (dj*k*(Y/r)*s+nu*i*j*c*X/(r^2))
    return bf, dx, dy
    end
end

"""
This function computes both the basis functions and their gradients with respect to the Cartesian coordinates `x` and `y` for multiple specified indices in the `CornerAdaptedFourierBessel` basis. The function evaluates these quantities at the provided points. This is a composite of `basis_fun` and `gradient`.

# Logic
- The points `pts` are mapped to local Cartesian coordinates using the `local_map` of the basis's coordinate system.
- These local Cartesian coordinates are converted to polar coordinates (`r`, `phi`).
- For each index in `indices`, the Bessel function `Jv(nu*i, k*r)` and its derivative with respect to `k*r` are computed.
- The basis functions `B` are calculated as the product of the Bessel function and the sine of the angular coordinate `nu*i*phi`.
- The gradient components `dB_dx` and `dB_dy` are calculated for each index using the chain rule, which involves the Bessel function, its derivative, and the trigonometric functions of `phi`.
- The function returns the basis functions and their gradients as three matrices: `B`, `dB_dx`, and `dB_dy`.

# Returns
- `B`: A matrix with dimensions `(number of points, number of indices)`, where each column contains the values of the basis function for a specific index.
- `dB_dx`: A matrix with dimensions `(number of points, number of indices)`, where each column contains the `x` component of the gradient of the basis function for a specific index.
- `dB_dy`: A matrix with dimensions `(number of points, number of indices)`, where each column contains the `y` component of the gradient of the basis function for a specific index.

# Arguments
- `basis`: The `CornerAdaptedFourierBessel` basis.
- `indices`: An array of indices for which the basis functions and gradients are to be computed.
- `k`: The wavenumber.
- `pts`: An array of points where the basis functions and gradients are to be evaluated.
"""
function basis_and_gradient(basis::CornerAdaptedFourierBessel, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map, nu=basis.nu, pts=pts
        #local cartesian coords
        pt_xy = collect(pm(pt) for pt in pts)
        pt_pol = collect(cartesian_to_polar(pt) for pt in pt_xy)
         
        #norm::T = one(T)/sqrt(basis.dim)
        r = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2)
        X = getindex.(pt_xy,1)
        Y = getindex.(pt_xy,2)
        M =  length(pts)
        N = length(indices)
        B = zeros(T,M,N)
        dB_dx = zeros(T,M,N)
        dB_dy = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            j = Jv.(nu*i, k*r)
            #println(size(j))
            dj = Jvp.(nu*i, k*r) 
            #println(size(dj))
            s = @. sin(nu*i*phi) 
            c = @. cos(nu*i*phi) 
            #println(size(s))
            B[:,i] .= @. j*s
            dB_dx[:,i] .= @. (dj*k*(X/r)*s-nu*i*j*c*Y/(r^2))
            dB_dy[:,i] .= @. (dj*k*(Y/r)*s+nu*i*j*c*X/(r^2))
        end
        #println(size(s))
    return B, dB_dx, dB_dy
    end
end

