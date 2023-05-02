using Bessels
using CoordinateTransformations, Rotations, StaticArrays

#aliases so one is able to change bessel library easily
H0(x) = Bessels.hankelh1(0, x) 
dH0(x) = -Bessels.hankelh1(1, x) 

struct FundamentalHankel{T,Sy} <: AbsFundamentalBasis where {T<:Real, Sy<:Union{AbsSymmetry,Nothing}}
    #cs::PolarCS{T} #not fully implemented
    dim::Int64 #using concrete type
    symmetries::Union{Vector{Sy},Nothing}
    return_type::Type
end

function FundamentalHankel(dim; type=Float64)
    return_type = Complex{type}
    return FundamentalHankel{type,Nothing}(dim, nothing, return_type)
end

function FundamentalHankel(dim, symmetries::Vector{Sy}; type=Float64) where {T<:Real, Sy<:Union{AbsSymmetry,Nothing}}
    #cs::PolarCS{T} #not fully implemented
    return_type = Complex{type}
    return FundamentalHankel{type,Sy}(dim, symmetries, return_type)
end

function resize_basis(basis::Ba, billiard::Bi, dim::Int, k) where {T<:Real, Sy<:Union{AbsSymmetry,Nothing}, Ba<:FundamentalHankel{T,Sy}, Bi<:AbsBilliard}
    return FundamentalHankel(dim, basis.symmetries)
end

@inline function basis_fun(basis::FundamentalHankel, arg::Union{T,Vector{T}}) where {T<:Real}
    h = H0.(arg)
    return -im*0.25.*h
end

@inline function derivative_fun(basis::FundamentalHankel, arg::Union{T,Vector{T}}) where {T<:Real}
    h = dH0.(arg)
    return -im*0.25*h
end

