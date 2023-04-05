#include("../abstracttypes.jl")
#include("../utils/typeutils.jl")

using Random, Distributions

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