include("../abstracttypes.jl")

using ForwardDiff, LinearAlgebra

function curve(t, x, y, p...; kw_params...)
    xx = [x(ti, p...; kw_params...) for ti in t]
    yy = [y(ti, p...; kw_params...) for ti in t]
    return xx, yy
end

function tangent(t, x, y, p...; kw_params...)
    tx = [ForwardDiff.derivative(t->x(t, p...; kw_params...), ti) for ti in t]
    ty = [ForwardDiff.derivative(t->y(t, p...; kw_params...), ti) for ti in t]
    return tx, ty
end

function tangent_length(t, x, y, p...; kw_params...)
    l = [norm([ForwardDiff.derivative(t->x(t, p...; kw_params...),ti), ForwardDiff.derivative(t->y(t, p...; kw_params...),ti)]) for ti in t]
    return l
end

function tangent_vec(t, x, y, p...; kw_params...)
    tx = [ForwardDiff.derivative(t->x(t, p...; kw_params...), ti) for ti in t]
    ty = [ForwardDiff.derivative(t->y(t, p...; kw_params...), ti) for ti in t]
    n = [norm([tx[i],ty[i]]) for i in 1:length(tx)]
    return tx ./ n, ty ./ n
end

function normal_vec(t, x, y, p...; kw_params...)
    tx, ty = tangent_vec(t, x, y, p...; kw_params...)
    return ty, -tx
end

#=
function normal(t, x, y, p...; kw_params...)
    tx, ty = tangent_vec(t, x, y, p...; kw_params...)
    return ty, -tx
end
=#
function curvature(t, x, y, p...; kw_params...)
    tx, ty = tangent_vec(t, x, y, p...; kw_params...)
    n = [norm([tx[i],ty[i]]) for i in 1:length(tx)]
    dtx = [ForwardDiff.derivative(t -> ForwardDiff.derivative(t->x(t, p...; kw_params...), t), ti) for ti in t]
    dty = [ForwardDiff.derivative(t -> ForwardDiff.derivative(t->y(t, p...; kw_params...), t), ti) for ti in t]
    n2 = [norm([dtx[i],dty[i]]) for i in 1:length(tx)]
    return n2 ./ n
end

struct BoundaryPoints{T} <: AbsPoints
    x :: Vector{T} 
    y :: Vector{T} 
    nx :: Vector{T} 
    ny :: Vector{T}
    s :: Vector{T}
    w :: Vector{T}
end

#abstract type AbsCurve end
#abstract type SymmetryAxis <: AbsCurve end
#=
struct BoundaryCurve  <: AbsCurve
    r::Function
    n::Function
    s::Function
    ds::Function
    orientation::Int
    length::Real
    #domain_function::Function
end

struct VirtualCurve  <: AbsCurve
    r::Function
    n::Function
    s::Function
    ds::Function
    orientation::Int
    length::Real
    #domain_function::Function
end
=#