#include("../abstracttypes.jl")
#include("../utils/billiardutils.jl")
#include("../utils/gridutils.jl")
#include("../solvers/matrixconstructors.jl")
using FFTW

#this takes care of singular points
"""
Regularize a vector by replacing `NaN` values with the average of their neighboring elements. This will remove singular points.

# Arguments
- `u::Vector`: The vector to be regularized, with potential `NaN` values.
"""
function regularize!(u)
    idx = findall(isnan, u)
    for i in idx
        u[i] = (u[i+1] + u[i-1])/2.0
    end
end

"""
Compute the boundary function u for a given quantum state.

For further reading: https://users.flatironinstitute.org/~ahb/thesis_html/node73.html - Tension matrix in a scaling eigenfunction basis

# Logic
This function computes the boundary function `u` for a given quantum state represented by `state`. It calculates the gradient of the state, regularizes it (gets rid of potential NaN's), and computes the boundary norm.
- With the `sampler` in the `state` variable we distribute points (and the normals, arclengths adn the integration weights) via the `boundary_coords` method
- For a given basis and the discretization points (from the previous bullet point) constructs the gradient matrixes (which are the gradients of the `basis` calculated for the different (x,y) points on the boundary)
- From the normals (also gotten from `boundary_coords` method) we calculate the dot product between the normal vectors and the gradient matrices for every boundary point (x,y)
- The resulting (nx, ny) ̇(dX, dY) is stored in `U` (not-normalized) which we regularize to help get rid of NaN values
- Finally, it computes the boundary `norm` by integrating the absolute square of the boundary function divided by 2k^2 (for this check the hyperlink above, this division comes from the introduction of the scaling eigenfunction basis)

# Arguments
- `state::AbsState`: The quantum state for which the boundary function is computed.
- `b::Float64`: Scaling factor for the number of sampling points (default is 5.0).

# Returns
- `u::Vector{T}`: The computed boundary function values.
- `s::Vector{T}`: The arc length coordinates corresponding to the boundary.
- `norm::T`: The computed boundary norm.
"""
function boundary_function(state::S; b=5.0) where {S<:AbsState}
    let vec = state.vec, k = state.k, k_basis = state.k_basis, new_basis = state.basis, billiard=state.billiard
        type = eltype(vec)
        boundary = billiard.full_boundary
        crv_lengths = [crv.length for crv in boundary]
        sampler = FourierNodes([2,3,5],crv_lengths)
        L = billiard.length
        N = max(round(Int, k*L*b/(2*pi)), 512)
        pts = boundary_coords(billiard, sampler, N)
        dX, dY = gradient_matrices(new_basis, k_basis, pts.xy)
        nx = getindex.(pts.normal,1)
        ny = getindex.(pts.normal,2)
        dX = nx .* dX 
        dY = ny .* dY
        U::Array{type,2} = dX .+ dY
        u::Vector{type} = U * vec
        regularize!(u)
        #compute the boundary norm
        w = dot.(pts.normal, pts.xy) .* pts.ds
        integrand = abs2.(u) .* w
        norm = sum(integrand)/(2*k^2)
        #println(norm)
        return u, pts.s::Vector{type}, norm
    end
end

"""
Compute the boundary functions `u` for a bundle of quantum states.

For further reading: https://users.flatironinstitute.org/~ahb/thesis_html/node73.html - Tension matrix in a scaling eigenfunction basis

# Logic
This function computes the boundary functions `u` for each quantum state in the bundle represented by `state_bundle`. It calculates the gradients of the states, regularizes them (eliminating potential NaN values), and computes the boundary norms for each state.
- Using the sampler provided in the `state_bundle`, points on the boundary (along with normals, arc lengths, and integration weights) are distributed via the `boundary_coords` method.
- For the given basis and the discretization points (from the previous step), the gradient matrices are constructed. These matrices represent the gradients of the basis functions evaluated at various (x, y) points on the boundary.
- The dot product between the normal vectors and the gradient matrices is computed for every boundary point (x, y). The results are stored in `U`, which is then regularized to eliminate NaN values.
- Finally, the boundary norm is computed by integrating the absolute square of the boundary function divided by 2k^2. This normalization stems from the scaling eigenfunction basis, as described in the linked reference.

# Arguments
- `state_bundle::EigenstateBundle`: The bundle of quantum states for which the boundary functions are computed.
- `b::Float64`: Scaling factor for the number of sampling points (default is 5.0).

# Returns
- `us::Vector{Vector{T}}`: The computed boundary function values for each state.
- `s::Vector{T}`: The arc length coordinates corresponding to the boundary.
- `norms::Vector{T}`: The computed boundary norms for each state.
"""
function boundary_function(state_bundle::S; b=5.0) where {S<:EigenstateBundle}
    let X = state_bundle.X, k_basis = state_bundle.k_basis, ks = state_bundle.ks, new_basis = state_bundle.basis, billiard=state_bundle.billiard 
        type = eltype(X)
        boundary = billiard.full_boundary
        crv_lengths = [crv.length for crv in boundary]
        sampler = FourierNodes([2,3,5],crv_lengths)
        L = billiard.length
        N = max(round(Int, k_basis*L*b/(2*pi)), 512)
        pts = boundary_coords(billiard, sampler, N)
        dX, dY = gradient_matrices(new_basis, k_basis, pts.xy)
        nx = getindex.(pts.normal,1)
        ny = getindex.(pts.normal,2)
        dX = nx .* dX 
        dY = ny .* dY
        U::Array{type,2} = dX .+ dY
        u_bundle::Matrix{type} = U * X
        for u in eachcol(u_bundle)
            regularize!(u)
        end
        #compute the boundary norm
        w = dot.(pts.normal, pts.xy) .* pts.ds
        norms = [sum(abs2.(u_bundle[:,i]) .* w)/(2*ks[i]^2) for i in eachindex(ks)]
        #println(norm)
        us::Vector{Vector{type}} = [u for u in eachcol(u_bundle)]
        return us, pts.s, norms
    end
end

"""
Compute the momentum distribution function from the boundary function `u`.

# Logic
This function computes the momentum distribution of the boundary function `u` using the Fourier transform. The Fourier transform of the boundary function is computed via `rfft` function. The corresponding wavenumbers `ks` are calculated based on the arc length coordinates `s`.
- Also the sampling rate is calculated as the inverse of the first difference between the arcslengths. THis is then used to calculate the `ks`.

# Arguments
- `u::Vector{T}`: The boundary function values.
- `s::Vector{T}`: The arc length coordinates corresponding to the boundary.

# Returns
- `mf::Vector{T}`: The momentum distribution function, which is the squared magnitude of the Fourier coefficients.
- `ks::Vector{T}`: The corresponding wavenumbers.
"""
function momentum_function(u,s)
    fu = rfft(u)
    sr = 1.0/diff(s)[1]
    ks = rfftfreq(length(s),sr).*(2*pi)
    return abs2.(fu)/length(fu), ks
end

"""
Compute the momentum distribution function for a quantum state.

# Logic
This function first calculates the boundary function `u` for the given quantum state and then computes its momentum distribution using the Fourier transform. The `momentum_function(u, s)` is called internally to perform this calculation.

# Arguments
- `state::AbsState`: The quantum state for which the momentum distribution is computed.
- `b::Float64`: Scaling factor for the number of sampling points (default is 5.0).

# Returns
- `mf::Vector{T}`: The momentum distribution function.
- `ks::Vector{T}`: The corresponding wavenumbers.
"""
function momentum_function(state::S; b=5.0) where {S<:AbsState}
    u, s, norm = boundary_function(state; b=b)
    return momentum_function(u,s)
end

#this can be optimized by usinf FFTW plans
"""
Compute the momentum distribution functions for a bundle of quantum states.

# Logic
This function computes the boundary functions `us` for each quantum state in the bundle and then calculates their respective momentum distributions. It iteratively calls `momentum_function(u, s)` for each state in the bundle. The function could be optimized using FFTW plans to handle the Fourier transforms more efficiently.

# Arguments
- `state_bundle::EigenstateBundle`: The bundle of quantum states for which the momentum distributions are computed.
- `b::Float64`: Scaling factor for the number of sampling points (default is 5.0).

# Returns
- `mfs::Vector{Vector{T}}`: A vector of momentum distribution functions for each state in the bundle.
- `ks::Vector{T}`: The wavenumbers corresponding to the momentum distributions.
"""
function momentum_function(state_bundle::S; b=5.0) where {S<:EigenstateBundle}
    us, s, norms = boundary_function(state_bundle; b=b)
    mf, ks = momentum_function(us[1],s)
    type = eltype(mf)
    mfs::Vector{Vector{type}} = [mf]
    for i in 2:length(us)
        mf, ks = momentum_function(us[i],s)
        push!(mfs,mf)
    end
    return mfs, ks
end