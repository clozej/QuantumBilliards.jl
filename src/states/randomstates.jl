#include("../abstracttypes.jl")
#include("../utils/typeutils.jl")

using Random, Distributions

"""
Construct a `GaussianRandomState` with a random vector sampled from a normal distribution.

# Description
This function creates a `GaussianRandomState` with a random vector whose elements are sampled from a standard normal distribution. The precision `eps` is set based on the float type of the wavenumber `k`.

# Arguments
- `k`: The wavenumber for the state.
- `dim`: The dimension of the vector.

# Returns
- A `GaussianRandomState` object containing the random vector and associated parameters.
"""
struct GaussianRandomState{K,T} <: AbsState where {K<:Number, T<:Real}
    k::K
    k_basis::K
    vec::Vector{T}
    dim::Int64
    eps::T
    #basis type
end

function GaussianRandomState(k,dim)
    d = Distributions.Normal()
    vec = rand(d, N)
    eps = set_precision(k)
    #norm = sum(abs.(vec))
    return GaussianRandomState(k,k, vec, dim,eps)
end