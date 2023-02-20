include("../abstracttypes.jl")
include("../utils/billiardutils.jl")
include("../utils/typeutils.jl")

struct BasisState{K,T} <: StationaryState where {K<:Number, T<:Real}
    k::K
    vec::Vector{T}
    dim::Int64
    eps::T
    #basis type
end

function BasisState(k, i, dim)  
    typ = typeof(k)
    eps = set_precision(k)
    vec = zeros(typ,dim)
    vec[i] = one(typ)
    return BasisState(k, vec, dim, eps)
end
