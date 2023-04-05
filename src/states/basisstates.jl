#include("../abstracttypes.jl")
#include("../utils/billiardutils.jl")
#include("../utils/typeutils.jl")

struct BasisState{K,T,Ba} <: StationaryState 
    k::K
    k_basis::K
    vec::Vector{T}
    idx::Int64
    dim::Int64
    eps::T
    basis::Ba
end

function BasisState(basis, k, i)  
    dim = basis.dim
    typ = typeof(k)
    eps = set_precision(k)
    vec = zeros(typ,dim)
    vec[i] = one(typ)
    return BasisState(k,k, vec, i, dim, eps, basis)
end
