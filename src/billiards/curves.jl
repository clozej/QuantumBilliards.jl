include("../abstracttypes.jl")
include("../billiards/coordinatesystems.jl")

using StaticArrays, ForwardDiff, LinearAlgebra

#curve convenience functions to transform coordinate systems
function curve(cs::C, f::Function, t::AbstractVector, p...; kw_params...) where C<:CoordinateSystem
    return [curve(cs,f,ti,p...; kw_params...) for ti in t]
end

function curve(cs::CartesianCS, f::Function, t::T, p...; kw_params...) where T<:Number
    pt = f(t,p...; kw_params...)
    return cs.affine_map(SVector(pt[1],pt[2]))
end

function curve(cs::PolarCS, f::Function, t::T, p...; kw_params...) where T<:Number
    pt = f(t,p...; kw_params...) #point in polar coordinates
    pt_xy = polar_to_cartesian(pt)
    return cs.affine_map(pt_xy)
end

#tangent vectors
function tangent(cs::C, f::Function, t::AbstractVector, p...; kw_params...) where C<:CoordinateSystem
    return [tangent(cs,f,ti,p...; kw_params...) for ti in t]
end

function tangent(cs::C, f::Function, t, p...; kw_params...) where C<:CoordinateSystem
    r = t -> curve(cs, f, t, p...; kw_params...)
    ta = ForwardDiff.derivative(r, t)
    return cs.affine_map(SVector(ta[1],ta[2]))
end

function tangent_vec(cs::C, f::Function, t::AbstractVector, p...; kw_params...) where C<:CoordinateSystem
    return [tangent_vec(cs,f,ti,p...; kw_params...) for ti in t]
end

function tangent_vec(cs::C, f::Function, t, p...; kw_params...) where C<:CoordinateSystem
    r = t -> curve(cs, f, t, p...; kw_params...)
    ta = ForwardDiff.derivative(r, t)
    norm = hypot(ta[1], ta[2])
    return SVector(ta[1]/norm,ta[2]/norm)
end


function normal_vec(cs::C, f::Function, t::AbstractVector, p...; kw_params...) where C<:CoordinateSystem
    return [normal_vec(cs,f,ti,p...; kw_params...) for ti in t]
end

function normal_vec(cs::C, f::Function, t, p...; kw_params...) where C<:CoordinateSystem
    ta = tangent_vec(cs, f, t, p...; kw_params...)
    return SVector(ta[2], -ta[1])
end


#=
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
=#


#=
function normal(t, x, y, p...; kw_params...)
    tx, ty = tangent_vec(t, x, y, p...; kw_params...)
    return ty, -tx
end
=#
#=
function curvature(t, x, y, p...; kw_params...)
    tx, ty = tangent_vec(t, x, y, p...; kw_params...)
    n = [norm([tx[i],ty[i]]) for i in 1:length(tx)]
    dtx = [ForwardDiff.derivative(t -> ForwardDiff.derivative(t->x(t, p...; kw_params...), t), ti) for ti in t]
    dty = [ForwardDiff.derivative(t -> ForwardDiff.derivative(t->y(t, p...; kw_params...), t), ti) for ti in t]
    n2 = [norm([dtx[i],dty[i]]) for i in 1:length(tx)]
    return n2 ./ n
end
=#
#=
struct BoundaryPoints{T} <: AbsPoints
    x :: Vector{T} 
    y :: Vector{T} 
    nx :: Vector{T} 
    ny :: Vector{T}
    s :: Vector{T}
    w :: Vector{T}
end
=#
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