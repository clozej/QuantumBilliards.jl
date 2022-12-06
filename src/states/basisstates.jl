include("../abstracttypes.jl")
include("../billiards/billiard.jl")

function set_precision(a)
    #expand for other types of numbers
    t = typeof(a)
    return t == Float32 ? Float32(1e-8) : convert(t,1e-16) 
end

struct Basisstate{K,T} <: StationaryState
    k::K
    vec::Vector{T}
    dim::Int
    eps::T
    #basis type
end

function Basisstate(k, i, dim)  
    typ = typeof(k)
    eps = set_precision(k)
    vec = zeros(typ,dim)
    vec[i] = one(typ)
    return Basisstate(k, vec, dim, eps)
end
