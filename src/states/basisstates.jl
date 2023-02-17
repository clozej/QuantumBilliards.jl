include("../abstracttypes.jl")
include("../billiards/billiard.jl")

function set_precision(a)
    #expand for other types of numbers
    t = typeof(a)
    return t == Float32 ? Float32(1e-8) : convert(t,1e-16) 
end

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
