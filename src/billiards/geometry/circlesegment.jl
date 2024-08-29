"""
Compute the coordinates of a point on a circle segment.

# Arguments
- `R::T`: The radius of the circle.
- `arc_angle::T`: The angular span of the circle segment.
- `shift_angle::T`: The angle by which the segment/point is shifted.
- `center::SVector{2,T}`: The center of the circle.
- `t`: Parameter `t` which varies between `0` and `1` to represent points along the segment.

# Returns
- `SVector{2,T}`: The 2D coordinates of the point on the circle segment.
"""
function circle_eq(R::T, arc_angle::T, shift_angle::T, center::SVector{2,T}, t) where T<:Real
    return SVector(R*cos(arc_angle*t+shift_angle) + center[1], R*sin(arc_angle*t+shift_angle)+center[2])
end

"""
Struct representing a segment of a circle.

# Fields
- `cs::PolarCS{T}`: The polar coordinate system associated with the segment.
- `radius::T`: The radius of the circle.
- `arc_angle::T`: The angular span of the segment.
- `shift_angle::T`: The angle by which the segment is shifted.
- `center::SVector{2,T}`: The center of the circle.
- `orientation::Int64`: The orientation of the segment (1 for clockwise, -1 for counterclockwise).
- `length::T`: The arc length of the segment.
"""
struct CircleSegment{T}  <: AbsRealCurve where T<:Real
    cs::PolarCS{T}
    radius::T
    arc_angle::T
    shift_angle::T
    center::SVector{2,T}
    orientation::Int64
    length::T
end

"""
Default constructor for the CircleSegment struct.

# Arguments
- `R::T`: The radius of the circle.
- `arc_angle::T`: The angular span of the circle segment.
- `shift_angle::T`: The angle by which the segment is shifted.
- `x0::T`: The x-coordinate of the center of the circle.
- `y0::T`: The y-coordinate of the center of the circle.
- `origin::SVector{2,T}`: The origin of the coordinate system (default is x=zero(x0),y=zero(x0))
- `rot_angle::T`: The rotation angle of the coordinate system (default is zero(x0)).
- `orientation::Int64` The orientation of the circle. By default it is clockwise (+1)

# Returns
- `CircleSegment{T}`: A CircleSegment object with the provided parameters.
"""
function CircleSegment(R, arc_angle, shift_angle, x0, y0; origin=(zero(x0),zero(x0)),rot_angle=zero(x0), orientation = 1)
    cs = PolarCS(SVector(origin...),rot_angle)
    center = SVector(x0,y0)
    L = R*arc_angle 
    return CircleSegment(cs,R,arc_angle,shift_angle,center,orientation,L)
end

"""
Struct representing a virtual segment of a circle, used in cases where the segment is part of a virtual boundary.

# Fields
- `cs::PolarCS{T}`: The polar coordinate system associated with the segment.
- `radius::T`: The radius of the circle.
- `arc_angle::T`: The angular span of the segment.
- `shift_angle::T`: The angle by which the segment is shifted.
- `center::SVector{2,T}`: The center of the circle.
- `orientation::Int64`: The orientation of the segment (1 for clockwise, -1 for counterclockwise).
- `length::T`: The arc length of the segment.
"""
struct VirtualCircleSegment{T}  <: AbsVirtualCurve where T<:Real
    cs::PolarCS{T}
    radius::T
    arc_angle::T
    shift_angle::T
    center::SVector{2,T}
    orientation::Int64
    length::T
end

"""
Default constructor for the VirtualCircleSegment struct.

# Arguments
- `R::T`: The radius of the circle.
- `arc_angle::T`: The angular span of the circle segment.
- `shift_angle::T`: The angle by which the segment is shifted.
- `x0::T`: The x-coordinate of the center of the circle.
- `y0::T`: The y-coordinate of the center of the circle.
- `origin::SVector{2,T}`: The origin of the coordinate system (default is x=zero(x0),y=zero(x0))
- `rot_angle::T`: The rotation angle of the coordinate system (default is zero(x0)).
- `orientation::Int64` The orientation of the circle. By default it is clockwise (+1)

# Returns
- `VirtualCircleSegment{T}`: A VirtualCircleSegment object with the provided parameters.
"""
function VirtualCircleSegment(R, arc_angle, shift_angle, x0, y0; origin=(zero(x0),zero(x0)),rot_angle=zero(x0), orientation = 1)
    cs = PolarCS(SVector(origin...),rot_angle)
    center = SVector(x0,y0)
    L = R*arc_angle 
    return VirtualCircleSegment(cs,R,arc_angle,shift_angle,center,orientation,L)
end

"""
Type that contains all `CircleSegment`s, be it real or virtual
"""
CircleSegments{T} = Union{CircleSegment{T},VirtualCircleSegment{T} } where T<:Real

"""
Calculate the coordinates of parameters `ts` for a circle segment.

# Arguments
- `circle::CircleSegments{T}`: The circle segment or virtual circle segment.
- `ts::AbstractArray{T,1}`: An array of parameter values between `0` and `1` representing points along the segment.

# Returns
- `Vector{SVector{2,T}}`: A vector of 2D coordinates of points on the circle segment.
"""
function curve(circle::L, ts::AbstractArray{T,1}) where {T<:Real,L<:CircleSegments{T}}
    let affine_map = circle.cs.affine_map, R = circle.radius, c = circle.center, a=circle.arc_angle, s=circle.shift_angle 
        return collect(affine_map(circle_eq(R, a, s, c, t)) for t in ts)
    end
end

"""
Calculate the coordinates of a single parameter `t` on a circle segment.

# Arguments
- `circle::CircleSegments{T}`: The circle segment or virtual circle segment.
- `t::T`: A parameter value between `0` and `1` representing a point along the segment.

# Returns
- `SVector{2,T}`: The 2D coordinates of the point on the circle segment.
"""
function curve(circle::L, t::T) where {T<:Real,L<:CircleSegments{T}}
    let affine_map = circle.cs.affine_map, R = circle.radius, c = circle.center, a=circle.arc_angle, s=circle.shift_angle 
        return affine_map(circle_eq(R, a, s, c, t))
    end
end

"""
Compute the arc length of points along a circle segment.

# Arguments
- `circle::CircleSegments{T}`: The circle segment or virtual circle segment.
- `ts::AbstractArray{T,1}`: An array of parameter values between `0` and `1` representing points along the segment.

# Returns
- `Vector{T}`: A vector of arc lengths corresponding to the parameter values.
"""
function arc_length(circle::L, ts::AbstractArray{T,1}) where {T<:Real,L<:CircleSegments{T}}
    s::Vector{T} = circle.length.*ts 
    return s
end

"""
Compute the tangent vectors at points along a circle segment. One needs to be careful with orientation in this case due to how the ForwardDiff operator is multiplied.

# Arguments
- `circle::CircleSegments{T}`: The circle segment or virtual circle segment.
- `ts::AbstractArray{T,1}`: An array of parameter values between `0` and `1` representing points along the segment.

# Returns
- `Vector{SVector{2,T}}`: A vector of tangent vectors at the specified points on the circle segment.
"""
function tangent(circle::L, ts::AbstractArray{T,1}) where {T<:Real,L<:CircleSegments{T}}
    let affine_map = circle.cs.affine_map, R = circle.radius, c = circle.center, a=circle.arc_angle, s=circle.shift_angle 
        orient = circle.orientation
        r(t) = affine_map(circle_eq(R, a, s, c, t))
        #ForwardDiff.derivative(r, t)
        return collect(orient*ForwardDiff.derivative(r, t) for t in ts)
    end
end

"""
Compute the signed distance from the circumference of a circle to a point or a set of points.

# Logic
- The function computes the Euclidean distance between each point `(x, y)` and the `center` of the circle.
- It then subtracts the radius `R` from this distance to obtain the signed distance to the circumference.
- If `x` and `y` are arrays, the operation is applied element-wise to each corresponding point using broadcasting (`@.` macro).

# Arguments
- `R::T`: Radius of the circle.
- `center::SVector{2,T}`: The center of the circle, given as a 2D point.
- `x`: The x-coordinate(s) of the point(s) to evaluate. Can be a scalar or an array.
- `y`: The y-coordinate(s) of the point(s) to evaluate. Can be a scalar or an array.

# Returns
- The signed distance(s) to the circle's circumference. Positive values indicate points outside the circle, and negative values indicate points inside the circle.
"""
circle_domain(R, center, x, y) = @. hypot(y-center[2],x-center[1]) - R

"""
Compute the signed distance from a circular segment to a point or set of points.

# Logic
- The function first computes the distance to the circular domain and the linear segment domain.
- Depending on the orientation, it returns the appropriate distance to the circle segment.
- If `x` and `y` are arrays, the operation is applied element-wise to each corresponding point using broadcasting.

# Arguments
- `pt0::SVector{2,T}`: Start point of the circular segment.
- `pt1::SVector{2,T}`: End point of the circular segment.
- `R::T`: Radius of the circle.
- `center::SVector{2,T}`: Center of the circle.
- `orient::Int64`: Orientation of the circular segment (1 for counterclockwise, -1 for clockwise).
- `x`: The x-coordinate(s) of the point(s) to evaluate. Can be a scalar or an array.
- `y`: The y-coordinate(s) of the point(s) to evaluate. Can be a scalar or an array.

# Returns
- The signed distance(s) to the circle segment. Positive values indicate points outside the segment, and negative values indicate points inside.
"""
function circle_segment_domain(pt0, pt1, R, center, orient, x, y)
    let cd = orient*circle_domain(R, center, x, y)
        ld = orient*line_domain(pt0[1],pt0[2],pt1[1],pt1[2],x,y)
        if orient > 0
            return ld < zero(ld) ? ld : cd
        else
            return ld > zero(ld) ? ld : cd
        end
    end
end
#=    
function domain(circle::L, pt::SVector{2,T}) where {T<:Real,L<:CircleSegments{T}}
    let pt0 = curve(circle,zero(T)), pt1 = curve(circle,one(T)), R=circle.radius, orient = circle.orientation
        center = circle.cs.affine_map(circle.center) #move center to correct position
        return circle_segment_domain(pt0, pt1, R, center, orient, pt[1], pt[2])*orient
    end
end
=#

"""
Compute the signed distance from a circular segment to multiple points.

# Logic
- The function calculates the distance of each point in `pts` from the circle segment defined by the `circle` object.
- It handles multiple points by using the `@.` macro to apply the calculations element-wise to each point in `pts`.

# Arguments
- `circle::L`: A `CircleSegment` or `VirtualCircleSegment` object defining the circular segment.
- `pts::AbstractArray`: An array of points, where each point is an `SVector{2,T}` representing its coordinates.

# Returns
- A vector of signed distances for each point. Positive values indicate points outside the segment, and negative values indicate points inside.
"""
function domain(circle::L, pts::AbstractArray) where {T<:Real,L<:CircleSegments{T}}
    let pt0 = curve(circle,zero(T)), pt1 = curve(circle,one(T)), R=circle.radius, orient = circle.orientation
        center = circle.cs.affine_map(circle.center) #move center to correct position
        return collect(circle_segment_domain(pt0, pt1, R, center, orient, pt[1], pt[2]) for pt in pts)
    end
end

"""
Determine if points are inside a circular segment.

# Logic
- The function computes the signed distance of each point in `pts` from the circle segment defined by the `circle` object.
- It returns `true` for points inside the segment and `false` for points outside.
- The function handles multiple points by using the `@.` macro to apply the calculations element-wise.

# Arguments
- `circle::L`: A `CircleSegment` or `VirtualCircleSegment` object defining the circular segment.
- `pts::AbstractArray`: An array of points, where each point is an `SVector{2,T}` representing its coordinates.

# Returns
- A vector of Boolean values indicating whether each point is inside the circle segment.
"""
function is_inside(circle::L, pts::AbstractArray) where {T<:Real,L<:CircleSegments{T}}
    let pt0 = curve(circle,zero(T)), pt1 = curve(circle,one(T)), R=circle.radius, orient = circle.orientation
        center = circle.cs.affine_map(circle.center) #move center to correct position
        return collect(circle_segment_domain(pt0, pt1, R, center, orient,pt[1],pt[2]) < zero(eltype(pt0)) for pt in pts)
    end
end
