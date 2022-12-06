include("../abstracttypes.jl")
using Random, Distributions

struct GaussianRandomState <: StationaryState
    k::Float64
    vec::Vector{Float64}
    dim::Int
    eps::Float64
    #basis type
    function GaussianRandomState(k,dim)
        d = Distributions.Normal()
        vec = rand(d, N)
        #norm = sum(abs.(vec))
        return new(k, vec, dim)
    end
end