include("../abstracttypes.jl")
include("../utils/coordinatesystems.jl")
#include("curves.jl")
using StaticArrays,LinearAlgebra, ForwardDiff

line_eq(pt0::SVector{2,T}, pt1::SVector{2,T}, t) where {T<:Number} = @. (pt1 - pt0) * t + pt0
line_domain(x0,y0,x1,y1,x,y) = ((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)

#line(t, x0, x1) = (x1-x0) * t + x0 
struct LineSegment{T}  <: AbsRealCurve where T<:Real
    cs::CartesianCS{T}
    pt0::SVector{2,T}
    pt1::SVector{2,T}
    orientation::Int64
    length::T
end

function LineSegment(pt0::SVector{2,T}, pt1::SVector{2,T}; origin=zero(pt0),rot_angle=zero(eltype(pt0)), orientation = 1) where T<:Real
    cs = CartesianCS(SVector(origin...),rot_angle)
    x, y = pt1 .- pt0        
    L = hypot(x,y)
    return LineSegment(cs,pt0,pt1,orientation,L)
end

struct VirtualLineSegment{T}  <: AbsVirtualCurve where T<:Real
    cs::CartesianCS{T}
    pt0::SVector{2,T}
    pt1::SVector{2,T}
    orientation::Int64
    length::T
end

function VirtualLineSegment(pt0::SVector{2,T}, pt1::SVector{2,T}; origin=zero(pt0),rot_angle=zero(eltype(pt0)), orientation = 1) where T<:Real
    cs = CartesianCS(SVector(origin...),rot_angle)
    x, y = pt1 .- pt0        
    L = hypot(x,y)
    return VirtualLineSegment(cs,pt0,pt1,orientation,L)
end

LineSegments{T} = Union{LineSegment{T},VirtualLineSegment{T} } where T<:Real


function curve(line::L, ts::AbstractArray{T,1}) where {T<:Real,L<:LineSegments{T}}
    let pt0 = line.cs.affine_map(line.pt0), pt1 = line.cs.affine_map(line.pt1)
    return collect(line_eq(pt0,pt1,t) for t in ts)
    end
end

function arc_length(line::L, ts::AbstractArray{T,1}) where {T<:Real,L<:LineSegments{T}}
    s::Vector{T} = line.length.*ts
    return s
end

function tangent(line::L, ts::AbstractArray{T,1}) where {T<:Real,L<:LineSegments{T}}
    let pt0 = line.cs.affine_map(line.pt0), pt1 = line.cs.affine_map(line.pt1), orient=line.orientation
        r(t) = line_eq(pt0,pt1,t)
        #ForwardDiff.derivative(r, t)
        return collect(orient*ForwardDiff.derivative(r, t) for t in ts)
    end
end


#domain fuction has node (zero) at the edge of the curve

#= worse version
function domain2(line::LineSegment{T}, pts::AbstractArray) where {T<:Real}
    let line = line
    return collect(domain(line,pt) for pt in pts)
    end
end
=#


function domain(line::L, pt::SVector{2,T}) where {T<:Real,L<:LineSegments{T}}
    let pt0 = line.cs.affine_map(line.pt0) 
        pt1 = line.cs.affine_map(line.pt1)
        orientation = line.orientation
    return line_domain(pt0[1],pt0[2],pt1[1],pt1[2],pt[1],pt[2])*orientation
    end
end

function domain(line::L, pts::AbstractArray) where {T<:Real,L<:LineSegments{T}}
    let pt0 = line.cs.affine_map(line.pt0) 
        pt1 = line.cs.affine_map(line.pt1)
        orientation = line.orientation
    return collect(line_domain(pt0[1],pt0[2],pt1[1],pt1[2],pt[1],pt[2])*orientation for pt in pts)
    end
end

function domain(line::L, x_grid::AbstractArray, y_grid::AbstractArray) where {T<:Real,L<:LineSegments{T}}
    let pt0 = line.cs.affine_map(line.pt0) 
        pt1 = line.cs.affine_map(line.pt1)
        orientation = line.orientation
    return collect(line_domain(pt0[1],pt0[2],pt1[1],pt1[2],x,y)*orientation for y in y_grid for x in x_grid)
    end
end

function is_inside(line::L, pts::AbstractArray) where {T<:Real,L<:LineSegments{T}}
    let pt0 = line.cs.affine_map(line.pt0) 
        pt1 = line.cs.affine_map(line.pt1)
        orientation = line.orientation
    return collect(line_domain(pt0[1],pt0[2],pt1[1],pt1[2],pt[1],pt[2])*orientation < zero(eltype(pt0)) for pt in pts)
    end
end

function is_inside(line::L, x_grid::AbstractArray, y_grid::AbstractArray) where {T<:Real,L<:LineSegments{T}}
    let pt0 = line.cs.affine_map(line.pt0) 
        pt1 = line.cs.affine_map(line.pt1)
        orientation = line.orientation
    return collect(line_domain(pt0[1],pt0[2],pt1[1],pt1[2],x,y)*orientation < zero(eltype(x_grid)) for y in y_grid for x in x_grid)
    end
end
#@.((y1-y0)*pt[1]-(x1-x0)*pt[2]+x1*y0-y1*x0)*orientation #negative inside

#function domain(line)
#typeof(range(0.0,1.0,10)) <: AbstractArray{Float64,1}
#line = LineSegment(SVector(0.0,0.0), SVector(1.0,1.0))
#@code_warntype r(line, range(0.0,1.0,10))





function circle_eq(R::T, arc_angle::T, shift_angle::T, center::SVector{2,T}, t) where T<:Real
    return SVector(R*cos(arc_angle*t+shift_angle) + center[1], R*sin(arc_angle*t+shift_angle)+center[2])
end


struct CircleSegment{T}  <: AbsRealCurve where T<:Real
    cs::PolarCS{T}
    radius::T
    arc_angle::T
    shift_angle::T
    center::SVector{2,T}
    orientation::Int64
    length::T
end

function CircleSegment(R, arc_angle, shift_angle, x0, y0; origin=(zero(x0),zero(x0)),rot_angle=zero(x0), orientation = 1)
    cs = PolarCS(SVector(origin...),rot_angle)
    center = SVector(x0,y0)
    L = R*arc_angle 
    return CircleSegment(cs,R,arc_angle,shift_angle,center,orientation,L)
end

struct VirtualCircleSegment{T}  <: AbsVirtualCurve where T<:Real
    cs::PolarCS{T}
    radius::T
    arc_angle::T
    shift_angle::T
    center::SVector{2,T}
    orientation::Int64
    length::T
end

function VirtualCircleSegment(R, arc_angle, shift_angle, x0, y0; origin=(zero(x0),zero(x0)),rot_angle=zero(x0), orientation = 1)
    cs = PolarCS(SVector(origin...),rot_angle)
    center = SVector(x0,y0)
    L = R*arc_angle 
    return VirtualCircleSegment(cs,R,arc_angle,shift_angle,center,orientation,L)
end

CircleSegments{T} = Union{CircleSegment{T},VirtualCircleSegment{T} } where T<:Real


function curve(circle::L, ts::AbstractArray{T,1}) where {T<:Real,L<:CircleSegments{T}}
    let affine_map = circle.cs.affine_map, R = circle.radius, c = circle.center, a=circle.arc_angle, s=circle.shift_angle 
        return collect(affine_map(circle_eq(R, a, s, c, t)) for t in ts)
    end
end

function curve(circle::L, t::T) where {T<:Real,L<:CircleSegments{T}}
    let affine_map = circle.cs.affine_map, R = circle.radius, c = circle.center, a=circle.arc_angle, s=circle.shift_angle 
        return affine_map(circle_eq(R, a, s, c, t))
    end
end

function arc_length(circle::L, ts::AbstractArray{T,1}) where {T<:Real,L<:CircleSegments{T}}    s = @. line.length * ts
    s::Vector{T} = circle.length.*ts 
    return s
end

function tangent(circle::L, ts::AbstractArray{T,1}) where {T<:Real,L<:CircleSegments{T}}
    let affine_map = circle.cs.affine_map, R = circle.radius, c = circle.center, a=circle.arc_angle, s=circle.shift_angle 
        orient = circle.orientation
        r(t) = affine_map(circle_eq(R, a, s, c, t))
        #ForwardDiff.derivative(r, t)
        return collect(orient*ForwardDiff.derivative(r, t) for t in ts)
    end
end

circle_domain(R, center, x, y) = @. hypot(y-center[2],x-center[1]) - R

function circle_segment_domain(pt0, pt1, R, center, orient, x, y)
    let cd = orient*circle_domain(R, center, x, y)
        ld = orient*line_domain(pt0[1],pt0[2],pt1[1],pt1[2],x,y)
        return ld < zero(ld) ? ld : cd
    end
end
    
function domain(circle::L, pt::SVector{2,T}) where {T<:Real,L<:CircleSegments{T}}
    let pt0 = curve(circle,zero(T)), pt1 = curve(circle,one(T)), R=circle.radius, orient = circle.orientation
        center = circle.cs.affine_map(circle.center) #move center to correct position
        return circle_segment_domain(pt0, pt1, R, center, orient, pt[1], pt[2])*orient
    end
end

function domain(circle::L, pts::AbstractArray) where {T<:Real,L<:CircleSegments{T}}
    let pt0 = curve(circle,zero(T)), pt1 = curve(circle,one(T)), R=circle.radius, orient = circle.orientation
        center = circle.cs.affine_map(circle.center) #move center to correct position
        return collect(circle_segment_domain(pt0, pt1, R, center, orient, pt[1], pt[2]) for pt in pts)
    end
end

function domain(circle::L, x_grid::AbstractArray, y_grid::AbstractArray) where {T<:Real,L<:CircleSegments{T}}
    let pt0 = curve(circle,zero(T)), pt1 = curve(circle,one(T)), R=circle.radius, orient = circle.orientation
        center = circle.cs.affine_map(circle.center) #move center to correct position
        return collect(circle_segment_domain(pt0, pt1, R, center, orient,x,y) for y in y_grid for x in x_grid)
    end
end

function is_inside(circle::L, pts::AbstractArray) where {T<:Real,L<:CircleSegments{T}}
    let pt0 = curve(circle,zero(T)), pt1 = curve(circle,one(T)), R=circle.radius, orient = circle.orientation
        center = circle.cs.affine_map(circle.center) #move center to correct position
        return collect(circle_segment_domain(pt0, pt1, R, center, orient,x,y) < zero(eltype(pt0)) for pt in pts)
    end
end

function is_inside(circle::L, x_grid::AbstractArray, y_grid::AbstractArray) where {T<:Real,L<:CircleSegments{T}}
    let pt0 = curve(circle,zero(T)), pt1 = curve(circle,one(T)), R=circle.radius, orient = circle.orientation
        center = circle.cs.affine_map(circle.center) #move center to correct position
        return collect(circle_segment_domain(pt0, pt1, R, center, orient,x,y) < zero(eltype(x_grid)) for y in y_grid for x in x_grid)
    end
end


#this part is general
function tangent_vec(curve::L, ts::AbstractArray{T,1}) where {T<:Real,L<:AbsCurve}
    ta = tangent(curve, ts)
    return collect(ti/norm(ti) for ti in ta)
end

function normal_vec(curve::L, ts::AbstractArray{T,1}) where {T<:Real,L<:AbsCurve}
    ta = tangent_vec(curve, ts)
    return [SVector(ti[2], -ti[1]) for ti in ta]
end
#=
function d(pt)
    circ =  d2(pt[1],pt[2])
    lin = d1(pt[1],pt[2])
    return circ .* lin .* sign(circ)
end
s(t) = L .* t
ds(t) = tangent_length(cs, f, t)
c_x0, c_y0 = r(0.0)
c_x1, c_y1 = r(1.0)
# line domain used to cut off circle and extend domain

=#

function make_polygon(corners, curve_types; origin=(zero(corners[1][1]),zero(corners[1][1])),rot_angle=zero(corners[1][1]))
    N = length(corners)
    boundary = []
    circular_idx(i) = mod1(i,N)
    #d(x, y, x0, y0, x1, y1) = @.((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)
    for i in 1:N
        idx0 = circular_idx(i)
        c0 = corners[idx0]
        c1 = corners[circular_idx(i+1)]
        line = (curve_types[idx0] == :Real) ? LineSegment(c0, c1;origin=origin,rot_angle=rot_angle) : VirtualLineSegment(c0,c1;origin=origin,rot_angle=rot_angle)
        push!(boundary, line)
    end
    return boundary
end

