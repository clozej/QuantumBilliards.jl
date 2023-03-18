
#include("../abstracttypes.jl")
#include("../basis/fourierbessel/corneradapted.jl")
#include("geometry.jl")

function make_quarter_stadium(half_width;radius=one(half_width),x0=zero(half_width),y0=zero(half_width),rot_angle=zero(half_width))
    #d(x, y, x0, y0, x1, y1) = @.((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)
    origin = SVector(x0,y0)
    type = typeof(half_width)
    circle = CircleSegment(radius,pi/2,zero(type), half_width, zero(type); origin=origin, rot_angle = rot_angle)
    corners = [SVector(half_width, radius), SVector(zero(type), radius), SVector(zero(type), zero(type)), SVector(half_width + radius, zero(type))]
    
    line1 = LineSegment(corners[1],corners[2];origin=origin,rot_angle=rot_angle)
    line2 = VirtualLineSegment(corners[1],corners[2];origin=origin,rot_angle=rot_angle)
    line3 = VirtualLineSegment(corners[1],corners[2];origin=origin,rot_angle=rot_angle)
    boundary = [circle, line1, line2, line3]
    return boundary, corners
end

function make_full_stadium(half_width;radius=one(half_width),x0=zero(half_width),y0=zero(half_width),rot_angle=zero(half_width))
    #d(x, y, x0, y0, x1, y1) = @.((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)
    origin = SVector(x0,y0)
    type = typeof(half_width)
    
    corners = [SVector(half_width, radius), SVector(-half_width, radius), SVector(-half_width, -radius), SVector(half_width, -radius)]
    circle1 = CircleSegment(radius,1.0*pi, -pi*0.5, half_width, zero(type); origin=origin, rot_angle = rot_angle)
    line1 = LineSegment(corners[1],corners[2];origin=origin,rot_angle=rot_angle)
    circle2 = CircleSegment(radius,1.0*pi, pi*0.5, -half_width, zero(type); origin=origin, rot_angle = rot_angle)
    line2 = LineSegment(corners[3],corners[4];origin=origin,rot_angle=rot_angle)
    boundary = [circle1, line1, circle2, line2]
    return boundary, corners
end

struct Stadium{T}  <: AbsBilliard where {T<:Real}
    boundary::Vector
    length::T
    area::T
    half_width::T
    radius::T
    corners::Vector{SVector{2,T}}
end

function Stadium(half_width;full_domain=false, radius=1.0,x0=0.0,y0=0.0)
    if full_domain
        boundary, corners = make_full_stadium(half_width;radius=radius,x0=x0,y0=y0)
        area = 4.0*half_width*radius + (pi*radius^2)
    else
        boundary, corners = make_quarter_stadium(half_width;radius=radius,x0=x0,y0=y0)
        area = half_width*radius + (pi*radius^2)/4.0 
    end
    length = sum([crv.length for crv in boundary])
    #PolygonOps.area(collect(zip(x,y)))
    return Stadium(boundary,length,area,half_width,radius,corners)
end 


