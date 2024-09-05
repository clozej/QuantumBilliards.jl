
#include("../abstracttypes.jl")
#include("../basis/fourierbessel/corneradapted.jl")
#include("geometry.jl")

"""
Constructs a quarter stadium billiard, which consists of a quarter circle connected to a straight line.

# Arguments
- `half_width::T`: The half-width of the stadium's straight section.
- `radius::T=one(half_width)`: The radius of the quarter circle section.
- `x0::T=zero(half_width)`: The x-coordinate of the origin.
- `y0::T=zero(half_width)`: The y-coordinate of the origin.
- `rot_angle::T=zero(half_width)`: The rotation angle of the billiard table.

# Returns
- A tuple containing:
- `boundary::Vector{Union{LineSegment{T}, CircleSegment{T}, VirtualLineSegment{T}}}`: The boundary segments of the quarter stadium.
- `corners::Vector{SVector{2,T}}`: The corner points of the quarter stadium.

# Logic
- The function starts by defining the circular segment as a quarter circle with the given `radius`, centered at the end of the straight section defined by `half_width`.
- The corners of the quarter stadium are calculated based on the input parameters, defining the a line attached to the quarter circle.
- Three line segments are created: one real segment and two virtual segments connecting the union of real curves to the origin.
"""
function make_quarter_stadium(half_width;radius=one(half_width),x0=zero(half_width),y0=zero(half_width),rot_angle=zero(half_width))
    #d(x, y, x0, y0, x1, y1) = @.((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)
    origin = SVector(x0,y0)
    type = typeof(half_width)
    circle = CircleSegment(radius,pi/2,zero(type), half_width, zero(type); origin=origin, rot_angle = rot_angle)
    corners = [SVector(half_width, radius), SVector(zero(type), radius), SVector(zero(type), zero(type)), SVector(half_width + radius, zero(type))]
    
    line1 = LineSegment(corners[1],corners[2];origin=origin,rot_angle=rot_angle)
    line2 = VirtualLineSegment(corners[2],corners[3];origin=origin,rot_angle=rot_angle)
    line3 = VirtualLineSegment(corners[3],corners[4];origin=origin,rot_angle=rot_angle)
    boundary = [circle, line1, line2, line3]
    return boundary, corners
end

"""
Constructs a full stadium billiard, which consists of two semicircles connected by two straight lines.

# Logic
- The function starts by defining two semicircular segments with the given `radius`, positioned at the ends of the straight section defined by `half_width`.
- The corners of the full stadium are calculated based on the input parameters, defining the 2 straight lines that connect the two semicircles.
- The function returns the boundary segments and the corner points in a tuple.

# Arguments
- `half_width::T`: The half-width of the stadium's straight section.
- `radius::T=one(half_width)`: The radius of the semicircular sections.
- `x0::T=zero(half_width)`: The x-coordinate of the origin.
- `y0::T=zero(half_width)`: The y-coordinate of the origin.
- `rot_angle::T=zero(half_width)`: The rotation angle of the billiard table.

# Returns
- A tuple containing:
  - `boundary::Vector{Union{LineSegment{T}, CircleSegment{T}}}`: The boundary segments of the full stadium.
  - `corners::Vector{SVector{2,T}}`: The corner points of the full stadium.
"""
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
    fundamental_boundary::Vector
    full_boundary::Vector
    length::T
    area::T
    half_width::T
    radius::T
    corners::Vector{SVector{2,T}}
end

"""
Constructs a stadium billiard with semicircular ends and a straight edge middle section.

# Arguments
- `half_width::T`: The half-width of the stadium's straight section.
- `radius::T=1.0`: The radius of the semicircular sections.
- `x0::T=0.0`: The x-coordinate of the origin.
- `y0::T=0.0`: The y-coordinate of the origin.

# Returns
- An instance of the `Stadium` struct representing the stadium billiard.

# Logic
- The constructor first creates the full stadium boundary and the corners using the `make_full_stadium` function.
- The `area` is calculated as the sum of the area of the rectangular section and the area of the two semicircles.
- The `fundamental_boundary` is constructed using the `make_quarter_stadium` function, representing the fundamental domain of the billiard.
- The `length` is calculated as the sum of the lengths of the segments in the `full_boundary`.
"""
function Stadium(half_width;radius=1.0,x0=0.0,y0=0.0)
    full_boundary, corners = make_full_stadium(half_width;radius=radius,x0=x0,y0=y0)
    area = 4.0*half_width*radius + (pi*radius^2)
    fundamental_boundary, _ = make_quarter_stadium(half_width;radius=radius,x0=x0,y0=y0)
    length = sum([crv.length for crv in full_boundary])
    #PolygonOps.area(collect(zip(x,y)))
    return Stadium(fundamental_boundary,full_boundary,length,area,half_width,radius,corners)
end 


