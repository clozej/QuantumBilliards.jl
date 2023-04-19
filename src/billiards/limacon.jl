using Elliptic
#include("curves.jl")

liamcon_radial(t, lam)  = (1.0 .+ lam .* cos.(t .* pi))

limacon_curve(t, lam, center::SVector{2,T}) = SVector(cos.(t .* pi) .* liamcon_radial(t, lam) .- center[1], sin.(t .* pi) .* liamcon_radial(t, lam) .- center[2])
limacon_arclength(t, lam) = 2*(1+lam) .* Elliptic.E.(t .* pi/2, 4*lam/(1+lam)^2)

struct LimaconCurve{T}  <: AbsRealCurve where T<:Real
    center::SVector{2,T}
    orientation::Int64
    length::T
end