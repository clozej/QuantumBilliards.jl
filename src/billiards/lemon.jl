
function make_quarter_lemon(half_separation;radius=one(half_separation),x0=zero(half_separation),y0=zero(half_separation),rot_angle=zero(half_separation))
    #d(x, y, x0, y0, x1, y1) = @.((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)
    origin = SVector(x0,y0)
    type = typeof(half_separation)
    theta = acos(half_separation/radius)
    shift = pi/2 -theta
    circle = CircleSegment(radius,theta, shift, zero(type), -half_separation; origin=origin, rot_angle = rot_angle)
    x_0 = radius*sin(theta)
    y_1 = radius-half_separation  
    corners = [SVector(x_0, zero(type)), SVector(zero(type), y_1), SVector(zero(type), zero(type))]
    
    line1 = VirtualLineSegment(corners[2],corners[3];origin=origin,rot_angle=rot_angle)
    line2 = VirtualLineSegment(corners[3],corners[1];origin=origin,rot_angle=rot_angle)
    boundary = [circle, line1, line2]
    return boundary, corners
end

function make_full_lemon(half_separation;radius=one(half_separation),x0=zero(half_separation),y0=zero(half_separation),rot_angle=zero(half_separation))
    #d(x, y, x0, y0, x1, y1) = @.((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)
    origin = SVector(x0,y0)
    type = typeof(half_separation)
    theta = acos(half_separation/radius)
    shift = pi/2 -theta
    circle1 = CircleSegment(radius,2*theta, shift, zero(type), -half_separation; origin=origin, rot_angle = rot_angle)
    circle2 = CircleSegment(radius,2*theta, pi+shift, zero(type), half_separation; origin=origin, rot_angle = rot_angle)

    x_0 = radius*sin(theta)
    corners = [SVector(x_0, zero(type)), SVector(-x_0, zero(type))]
    
    boundary = [circle1, circle2]
    return boundary, corners
end

struct Lemon{T}  <: AbsBilliard where {T<:Real}
    fundamental_boundary::Vector
    full_boundary::Vector
    length::T
    area::T
    half_separation::T
    radius::T
    corners::Vector{SVector{2,T}}
end

function Lemon(half_separation; radius=1.0,x0=0.0,y0=0.0)
    theta = acos(half_separation/radius)
    area = radius^2*(2*theta - sin(2*theta))
    full_boundary, _ = make_full_lemon(half_separation;radius=radius,x0=x0,y0=y0)
    fundamental_boundary, corners = make_quarter_lemon(half_separation;radius=radius,x0=x0,y0=y0)
    length = sum([crv.length for crv in full_boundary])
    #PolygonOps.area(collect(zip(x,y)))
    return Lemon(fundamental_boundary,full_boundary,length,area,half_separation,radius,corners)
end 