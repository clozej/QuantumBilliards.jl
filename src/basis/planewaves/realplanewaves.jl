"""
Basic struct containinng information for constructing the `RealPlaneWaves` basis

# Fields
- `dim::Int64`: The dimension of the basis.
- `symmetries::Union{Vector{Sy},Nothing}`: An optional vector of symmetries applied to the basis functions, where `Sy` is a subtype of `AbsSymmetry`.
- `angle_arc::T`: The angle used to define the scope of the plane waves.
- `angle_shift::T`: A shift applied to the angles.
- `angles::Vector{T}`: A vector of angles defining the directions of the plane waves.
- `parity_x::Vector{Int64}`: Parity values for the `x` direction, affecting the sine/cosine wave functions.
- `parity_y::Vector{Int64}`: Parity values for the `y` direction, affecting the sine/cosine wave functions.
- `sampler::Sa`: A sampler used to define the discretization (the basis indexes) of the angles, where `Sa` is a subtype of `AbsSampler`.
"""
struct RealPlaneWaves{T,Sy,Sa} <: AbsBasis where  {T<:Real, Sy<:Union{AbsSymmetry,Nothing}, Sa<:AbsSampler}
    #cs::PolarCS{T} #not fully implemented
    dim::Int64 #using concrete type
    symmetries::Union{Vector{Sy},Nothing}
    angle_arc::T
    angle_shift::T
    angles::Vector{T}
    parity_x::Vector{Int64}
    parity_y::Vector{Int64}
    sampler::Sa
end

"""
Provided (parity_x, parity_y) vectors that take into account all the symmetries of the system, depending on the axis chosen (determined internally by the :x_axis and :y_axis symbols). If not symmetries are present this method simply returns:  (parity_x = [1,1,-1,-1], parity_y = [1,-1,1,-1]). Otherwise it modifies parity_x if any :y_axis symmetry is present and modified the parity_y if any :x_axis symmetry is present and then returns them.
"""
function parity_pattern(symmetries)
    if isnothing(symmetries)
        parity_x = [1,1,-1,-1]
        parity_y = [1,-1,1,-1]
    else
        parity_x = [1,-1]
        parity_y = [1,-1]
        for sym in symmetries
            if sym.axis == :y_axis
                parity_x = [sym.parity, sym.parity]
            end
            if sym.axis == :x_axis
                parity_y = [sym.parity, sym.parity]   
            end
        end
    end            
     
    return parity_x, parity_y
end

"""
This function creates a `RealPlaneWaves` basis with a specified dimension, symmetries, and additional parameters for angle arc, angle shift, and sampler.

# Arguments
- `dim`: The dimension of the basis.
- `symmetries`: A vector of symmetries to be applied to the basis functions.
- `angle_arc`: The arc of angles used to define the spread of plane waves (default: `pi`).
- `angle_shift`: A shift applied to the angles (default: `0.0`).
- `sampler`: A sampler used to define the discretization of the angles (default: `LinearNodes()`).

# Returns
- A `RealPlaneWaves` object initialized with the specified parameters and symmetries.
"""
function RealPlaneWaves(dim, symmetries; angle_arc = pi, angle_shift=0.0, sampler=LinearNodes())
    par_x, par_y = parity_pattern(symmetries)
    pl = length(par_x)
    eff_dim = dim*pl
    t, dt = sample_points(sampler, dim)
    angles = @. t*angle_arc + angle_shift
    angles = repeat(angles, inner=pl)
    par_x = repeat(par_x, outer=dim)
    par_y = repeat(par_y, outer=dim)
    return RealPlaneWaves(eff_dim, symmetries, angle_arc, angle_shift, angles, par_x, par_y, sampler)
end

"""
This function creates a `RealPlaneWaves` basis with a specified dimension and additional parameters for angle arc, angle shift, and sampler. The basis functions are defined by angles adjusted by the angle shift. The parity patterns are set to default values as no symmetries are applied.

# Arguments
- `dim`: The dimension of the basis.
- `angle_arc`: The arc of angles used to define the plane waves (default: `pi`).
- `angle_shift`: A shift applied to the angles (default: `0.0`).
- `sampler`: A sampler used to define the discretization of the angles (default: `LinearNodes()`).

# Returns
- A `RealPlaneWaves` object initialized with the specified parameters.
"""
function RealPlaneWaves(dim; angle_arc = pi, angle_shift=0.0, sampler=LinearNodes())
    symmetries = nothing
    par_x, par_y = parity_pattern(symmetries)
    pl = length(par_x)
    eff_dim = dim*pl
    t, dt = sample_points(sampler, dim)
    angles = @. t*angle_arc + angle_shift
    angles = repeat(angles, inner=pl)
    par_x = repeat(par_x, outer=dim)
    par_y = repeat(par_y, outer=dim)
    return RealPlaneWaves{eltype(angles),Nothing,typeof(sampler)}(eff_dim, symmetries, angle_arc, angle_shift, angles, par_x, par_y, sampler)
end

"""
This function resizes the `RealPlaneWaves` basis to a specified dimension based on the given wavenumber `k` and retains the original angle arc, angle shift, and sampler settings. The resized basis includes updated angles and parities based on the new dimension.

# Arguments
- `basis::Ba`: The `RealPlaneWaves` basis to be resized.
- `billiard::Bi`: The billiard domain (not directly used here but maintained for consistency).
- `dim::Int`: The new dimension for the basis.
- `k`: The wavenumber (not directly used in this function).

# Returns
- A `RealPlaneWaves` object resized to the specified dimension with updated angles and parities.
"""
function resize_basis(basis::Ba, billiard::Bi, dim::Int, k) where {Ba<:RealPlaneWaves,Bi<:AbsBilliard}
    return RealPlaneWaves(dim, basis.symmetries; angle_arc = basis.angle_arc, angle_shift=basis.angle_shift, sampler=basis.sampler)
end

"""
This function computes the real plane wave component (either a cosine or sine function) based on the parity value. If the parity is `1`, it returns the cos of the argument. If the parity is `-1`, it returns the sin of the argument.

# Returns
- A vector of the computed real plane wave component (either cosine or sine) for each element in `arg`.

# Arguments
- `arg`: The argument to the sine or cosine function, typically a product of the wavenumber and a coordinate.
- `parity::Int64`: The parity value determining whether to use cosine (`1`) or sine (`-1`).
"""
@inline function rpw(arg, parity::Int64)
    if parity == 1
        return cos.(arg)
    else
        return sin.(arg)
    end 
end

"""
This function computes the derivative of the real plane wave component (either a sine or cosine function) based on the parity value. If the parity is `1`, it returns the negative sine of the argument (the derivative of cosine). If the parity is `-1`, it returns the cosine of the argument (the derivative of sine).

# Returns
- A vector of the computed derivative of the real plane wave component (either `-sin(arg)` or `cos(arg)`) for each element in `arg`.

# Arguments
- `arg`: The argument to the sine or cosine function, typically a product of the wavenumber and a coordinate.
- `parity::Int64`: The parity value determining whether to differentiate cosine (`1`) or sine (`-1`).
"""
@inline function d_rpw(arg, parity::Int64)
    if parity == 1
        return -sin.(arg)
    else
        return cos.(arg)
    end 
end

"""
This function computes the real plane wave basis function for a specified index `i` in the `RealPlaneWaves` basis at a given wavenumber `k`. The function evaluates the product of the real plane wave components (either sine or cosine) along the `x` and `y` directions based on the parity values.

# Logic
- The `x` and `y` components of the points `pts` are extracted.
- The direction cosines `vx` and `vy` are computed from the angle associated with the basis function at index `i`.
- The arguments for the sin or cos functions, `arg_x` and `arg_y`, are computed as the product of the wavenumber `k` and the corresponding direction cosines and coordinates (dot product).
- The basis function is computed as the product of the real plane wave components for `x` and `y`, determined by the parity values `par_x` and `par_y`.

# Returns
- A vector containing the values of the real plane wave basis function evaluated at the input points.

# Arguments
- `basis`: The `RealPlaneWaves` basis.
- `i`: The index of the basis function.
- `k`: The wavenumber.
- `pts`: An array of points where the basis function is to be evaluated.
"""
@inline function basis_fun(basis::RealPlaneWaves, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let par_x = basis.parity_x, par_y = basis.parity_y
        x = getindex.(pts,1)
        y = getindex.(pts,2)
        vx = cos(basis.angles[i])
        vy = sin(basis.angles[i])
        arg_x = k*vx.*x
        arg_y = k*vy.*y
        b = rpw(arg_x, par_x[i]).*rpw(arg_y, par_y[i])
        return b
    end
end

"""
This function computes the real plane wave basis functions for multiple specified indices in the `RealPlaneWaves` basis at a given wavenumber `k`. The function evaluates the product of the real plane wave components (either sine or cosine) along the `x` and `y` directions for each index based on the parity values.

# Logic
- The `x` and `y` components of the points `pts` are extracted.
- For each index in `indices`:
  - The direction cosines `vx` and `vy` are computed from the angle associated with the current basis function.
  - The arguments for the sin or cos functions, `arg_x` and `arg_y`, are computed as the product of the wavenumber `k` and the corresponding direction cosines and coordinates.
  - The basis function is computed as the product of the real plane wave components for `x` and `y`, determined by the parity values `par_x` and `par_y`.
- The function returns the basis functions as a matrix where each column corresponds to an index.

# Returns
- A matrix with dimensions `(number of points, number of indices)`, where each column contains the values of the real plane wave basis function for a specific index evaluated at the input points.

# Arguments
- `basis`: The `RealPlaneWaves` basis.
- `indices`: An array of indices for which the basis functions are to be computed.
- `k`: The wavenumber.
- `pts`: An array of points where the basis functions are to be evaluated.
"""
@inline function basis_fun(basis::RealPlaneWaves, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let par_x = basis.parity_x, par_y = basis.parity_y
        x = getindex.(pts,1)
        y = getindex.(pts,2)
        M =  length(pts)
        N = length(indices)
        B = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            vx = cos(basis.angles[i])
            vy = sin(basis.angles[i])
            arg_x = k*vx.*x
            arg_y = k*vy.*y
            B[:,i] .= rpw(arg_x, par_x[i]).*rpw(arg_y, par_y[i])
        end
        return B 
    end
end

"""
This function computes the gradient of the plane wave basis function with respect to the Cartesian coordinates `x` and `y` for a specified index `i` in the `RealPlaneWaves` basis. The function evaluates the derivatives of the sin or cos functions, depending on the parity, at each point in the input array.

# Arguments
- `basis`: The `RealPlaneWaves` basis.
- `i`: The index of the basis function.
- `k`: The wavenumber.
- `pts`: An array of points where the gradient is to be evaluated.

# Returns
- Two vectors: `dx` and `dy`, representing the gradient of the plane wave basis function with respect to the Cartesian coordinates `x` and `y`.
"""
function gradient(basis::RealPlaneWaves, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let par_x = basis.parity_x, par_y = basis.parity_y
        x = getindex.(pts,1)
        y = getindex.(pts,2)
        vx = cos(basis.angles[i])
        vy = sin(basis.angles[i])
        arg_x = k*vx.*x
        arg_y = k*vy.*y
        bx = rpw(arg_x, par_x[i])
        by = rpw(arg_y, par_y[i])
        dx = k*vx.*d_rpw(arg_x, par_x[i]).*by
        dy = bx.*k*vy.*d_rpw(arg_y, par_y[i])
        return dx, dy
    end
end

"""
This function computes the gradients of the plane wave basis functions with respect to the Cartesian coordinates `x` and `y` for multiple specified indices in the `RealPlaneWaves` basis. The function evaluates the derivatives of the sin or cos functions, depending on the parity, at each point in the input array for each index.

# Arguments
- `basis`: The `RealPlaneWaves` basis.
- `indices`: An array of indices for which the gradients are to be computed.
- `k`: The wavenumber.
- `pts`: An array of points where the gradients are to be evaluated.

# Returns
- Two matrices: `dB_dx` and `dB_dy`, representing the gradients of the plane wave basis functions with respect to the Cartesian coordinates `x` and `y`.
"""
function gradient(basis::RealPlaneWaves, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let par_x = basis.parity_x, par_y = basis.parity_y
        x = getindex.(pts,1)
        y = getindex.(pts,2)
        M =  length(pts)
        N = length(indices)
        dB_dx = zeros(T,M,N)
        dB_dy = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            vx = cos(basis.angles[i])
            vy = sin(basis.angles[i])
            arg_x = k*vx.*x
            arg_y = k*vy.*y
            bx = rpw(arg_x, par_x[i])
            by = rpw(arg_y, par_y[i])
            dB_dx[:,i] .= k*vx.*d_rpw(arg_x, par_x[i]).*by
            dB_dy[:,i] .= bx.*k*vy.*d_rpw(arg_y, par_y[i])
        end
        return dB_dx, dB_dy
    end
end

"""
This function computes both the plane wave basis function and its gradient with respect to the Cartesian coordinates `x` and `y` for a specified index `i` in the `RealPlaneWaves` basis. The function evaluates these quantities at the provided points, considering the sin or cos functions depending on the parity.

# Logic
- The points `pts` are separated into their `x` and `y` components.
- The direction cosines `vx` and `vy` are computed from the angle associated with the basis function at index `i`.
- The arguments for the sin or cos functions, `arg_x` and `arg_y`, are computed as the product of the wavenumber `k` and the corresponding direction cosines and coordinates.
- The basis function is computed as the product of the sine or cosine functions for `x` and `y`, depending on the parity values `par_x` and `par_y`.
- The gradient components `dx` and `dy` are computed by differentiating the sin or cos functions with respect to `x` and `y`.

# Returns
- `bf`: A vector containing the values of the plane wave basis function evaluated at the input points.
- `dx`: A vector containing the `x` component of the gradient of the basis function.
- `dy`: A vector containing the `y` component of the gradient of the basis function.

# Arguments
- `basis`: The `RealPlaneWaves` basis.
- `i`: The index of the basis function.
- `k`: The wavenumber.
- `pts`: An array of points where the basis function and its gradient are to be evaluated.
"""
function basis_and_gradient(basis::RealPlaneWaves, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let par_x = basis.parity_x, par_y = basis.parity_y
        x = getindex.(pts,1)
        y = getindex.(pts,2)
        vx = cos(basis.angles[i])
        vy = sin(basis.angles[i])
        arg_x = k*vx.*x
        arg_y = k*vy.*y
        bx = rpw(arg_x, par_x[i])
        by = rpw(arg_y, par_y[i])
        bf = bx.*by
        dx = k*vx.*d_rpw(arg_x, par_x[i]).*by
        dy = bx.*k*vy.*d_rpw(arg_y, par_y[i])
        return bf, dx, dy
    end
end

"""
This function computes both the plane wave basis functions and their gradients with respect to the Cartesian coordinates `x` and `y` for multiple specified indices in the `RealPlaneWaves` basis. The function evaluates these quantities at the provided points, considering the sin or cos functions depending on the parity.

# Logic
- The points `pts` are separated into their `x` and `y` components.
- For each index in `indices`:
- The direction cosines `vx` and `vy` are computed from the angle associated with the current basis function.
- The arguments for the sine or cosine functions, `arg_x` and `arg_y`, are computed as the product of the wavenumber `k` and the corresponding direction cosines and coordinates.
- The basis function is computed as the product of the sine or cosine functions for `x` and `y`, depending on the parity values `par_x` and `par_y`.
- The gradient components `dB_dx` and `dB_dy` are computed by differentiating the sine or cosine functions with respect to `x` and `y`, respectively.

# Returns
- `B`: A matrix with dimensions `(number of points, number of indices)`, where each column contains the values of the plane wave basis function for a specific index evaluated at the input points.
- `dB_dx`: A matrix with dimensions `(number of points, number of indices)`, where each column contains the `x` component of the gradient of the basis function for a specific index.
- `dB_dy`: A matrix with dimensions `(number of points, number of indices)`, where each column contains the `y` component of the gradient of the basis function for a specific index.

# Arguments
- `basis`: The `RealPlaneWaves` basis.
- `indices`: An array of indices for which the basis functions and gradients are to be computed.
- `k`: The wavenumber.
- `pts`: An array of points where the basis functions and gradients are to be evaluated.
"""
function basis_and_gradient(basis::RealPlaneWaves, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let par_x = basis.parity_x, par_y = basis.parity_y
        x = getindex.(pts,1)
        y = getindex.(pts,2)
        M =  length(pts)
        N = length(indices)
        B = zeros(T,M,N)
        dB_dx = zeros(T,M,N)
        dB_dy = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            vx = cos(basis.angles[i])
            vy = sin(basis.angles[i])
            arg_x = k*vx.*x
            arg_y = k*vy.*y
            bx = rpw(arg_x, par_x[i])
            by = rpw(arg_y, par_y[i])
            B[:,i] .= bx.*by
            dB_dx[:,i] .= k*vx.*d_rpw(arg_x, par_x[i]).*by
            dB_dy[:,i] .= bx.*k*vy.*d_rpw(arg_y, par_y[i])

        end
        return B, dB_dx, dB_dy
    end
end

"""
This function computes the derivative of the plane wave basis function with respect to the wavenumber `k` for a specified index `i` in the `RealPlaneWaves` basis. The function evaluates the product of the derivatives of the sin or cos functions, depending on the parity, at each point in the input array.

# Arguments
- `basis`: The `RealPlaneWaves` basis.
- `i`: The index of the basis function.
- `k`: The wavenumber.
- `pts`: An array of points where the derivative is to be evaluated.

# Returns
- A vector containing the derivative of the plane wave basis function with respect to the wavenumber `k`.
"""
@inline function dk_fun(basis::RealPlaneWaves, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let par_x = basis.parity_x, par_y = basis.parity_y
        x = getindex.(pts,1)
        y = getindex.(pts,2)
        vx = cos(basis.angles[i])
        vy = sin(basis.angles[i])
        arg_x = k*vx.*x
        arg_y = k*vy.*y
        bx = rpw(arg_x, par_x[i])
        by = rpw(arg_y, par_y[i])
        d_bx = d_rpw(arg_x, par_x[i])
        d_by = d_rpw(arg_y, par_y[i])
        dk = @. vx*x*d_bx*by + bx*vy*y*d_by
        return dk
    end
end
    
"""
This function computes the derivatives of the plane wave basis functions with respect to the wavenumber `k` for multiple specified indices in the `RealPlaneWaves` basis. The function evaluates the product of the derivatives of the sin or cos functions, depending on the parity, at each point in the input array for each index.

# Arguments
- `basis`: The `RealPlaneWaves` basis.
- `indices`: An array of indices for which the derivatives are to be computed.
- `k`: The wavenumber.
- `pts`: An array of points where the derivatives are to be evaluated.

# Returns
- A matrix `dB_dk` with dimensions `(number of points, number of indices)`, where each column contains the derivative of the plane wave basis function with respect to the wavenumber `k` for a specific index evaluated at the input points.
"""
@inline function dk_fun(basis::RealPlaneWaves, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let par_x = basis.parity_x, par_y = basis.parity_y
        x = getindex.(pts,1)
        y = getindex.(pts,2)
        M =  length(pts)
        N = length(indices)
        dB_dk = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            vx = cos(basis.angles[i])
            vy = sin(basis.angles[i])
            arg_x = k*vx.*x
            arg_y = k*vy.*y
            bx = rpw(arg_x, par_x[i])
            by = rpw(arg_y, par_y[i])
            d_bx = d_rpw(arg_x, par_x[i])
            d_by = d_rpw(arg_y, par_y[i])
            dB_dk[:,i] .=  @. vx*x*d_bx*by + bx*vy*y*d_by
        end
        return dB_dk
    end
end
