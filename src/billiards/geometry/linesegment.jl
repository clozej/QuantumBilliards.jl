"""
Compute the equation of a line segment at a parameter `t`.

# Logic
- This function computes the position along a line segment defined by two points `pt0` and `pt1` using the parameter `t`, where `t` ranges from 0 to 1.
- The result is a point along the line segment corresponding to the value of `t`.

# Arguments
- `pt0::SVector{2,T}`: The starting point of the line segment.
- `pt1::SVector{2,T}`: The ending point of the line segment.
- `t`: The parameter value (scalar or array) representing the position along the line segment.

# Returns
- The position on the line segment at parameter `t`.
"""
line_eq(pt0::SVector{2,T}, pt1::SVector{2,T}, t) where {T<:Number} = @. (pt1 - pt0) * t + pt0

"""
Compute the signed distance from a point to a line defined by two points.

# Logic
- This function calculates the signed distance from a point `(x, y)` to the line passing through `(x0, y0)` and `(x1, y1)`.
- The result is the value of the line equation evaluated at the point `(x, y)`.

# Arguments
- `x0`: The x-coordinate of the first point defining the line.
- `y0`: The y-coordinate of the first point defining the line.
- `x1`: The x-coordinate of the second point defining the line.
- `y1`: The y-coordinate of the second point defining the line.
- `x`: The x-coordinate(s) of the point(s) to evaluate. Can be a scalar or an array.
- `y`: The y-coordinate(s) of the point(s) to evaluate. Can be a scalar or an array.

# Returns
- The signed distance(s) to the line. Positive values indicate points on one side of the line, and negative values indicate points on the other side.
"""
line_domain(x0,y0,x1,y1,x,y) = ((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)

#line(t, x0, x1) = (x1-x0) * t + x0 
struct LineSegment{T}  <: AbsRealCurve where T<:Real
    cs::CartesianCS{T}
    pt0::SVector{2,T}
    pt1::SVector{2,T}
    orientation::Int64
    length::T
end

"""
Construct a `LineSegment` object representing a line segment in a Cartesian coordinate system.

# Logic
- This function initializes a `LineSegment` object with specified endpoints `pt0` and `pt1`, and optional origin, rotation angle, and orientation.
- It calculates the length of the line segment using the Euclidean distance between the endpoints.

# Arguments
- `pt0::SVector{2,T}`: The starting point of the line segment.
- `pt1::SVector{2,T}`: The ending point of the line segment.
- `origin`: The origin of the Cartesian coordinate system (default is the zero vector).
- `rot_angle`: The rotation angle of the Cartesian coordinate system (default is 0).
- `orientation`: The orientation of the line segment (default is 1).

# Returns
- A `LineSegment{T}` object representing the line segment.
"""
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

"""
Construct a `VirtualLineSegment` object representing a virtual line segment in a Cartesian coordinate system. Used for desymmetrization.

# Logic
- This function initializes a `VirtualLineSegment` object with specified endpoints `pt0` and `pt1`, and optional origin, rotation angle, and orientation.
- It calculates the length of the line segment using the Euclidean distance between the endpoints.

# Arguments
- `pt0::SVector{2,T}`: The starting point of the line segment.
- `pt1::SVector{2,T}`: The ending point of the line segment.
- `origin`: The origin of the Cartesian coordinate system (default is the zero vector).
- `rot_angle`: The rotation angle of the Cartesian coordinate system (default is 0).
- `orientation`: The orientation of the line segment (default is 1).

# Returns
- A `VirtualLineSegment{T}` object representing the virtual line segment.
"""
function VirtualLineSegment(pt0::SVector{2,T}, pt1::SVector{2,T}; origin=zero(pt0),rot_angle=zero(eltype(pt0)), orientation = 1) where T<:Real
    cs = CartesianCS(SVector(origin...),rot_angle)
    x, y = pt1 .- pt0        
    L = hypot(x,y)
    return VirtualLineSegment(cs,pt0,pt1,orientation,L)
end

"""
Type that contains all `LineSegment`s, be it real or virtual
"""
LineSegments{T} = Union{LineSegment{T},VirtualLineSegment{T} } where T<:Real

"""
Compute the points along a line segment for given parameter values.

# Logic
- This function computes the positions along the line segment defined by the `line` object for each value in `ts`.
- The positions are calculated by applying the affine transformation to the line segment's endpoints and then using the `line_eq` function.

# Arguments
- `line::L`: A `LineSegment` or `VirtualLineSegment` object defining the line segment.
- `ts::AbstractArray{T,1}`: An array of parameter values, typically ranging from 0 to 1.

# Returns
- An array of points on the line segment corresponding to the parameter values in `ts`.
"""
function curve(line::L, ts::AbstractArray{T,1}) where {T<:Real,L<:LineSegments{T}}
    let pt0 = line.cs.affine_map(line.pt0), pt1 = line.cs.affine_map(line.pt1)
    return collect(line_eq(pt0,pt1,t) for t in ts)
    end
end

"""
Compute the arc length along a line segment for given parameter values.

# Logic
- This function calculates the arc length along the line segment defined by the `line` object for each value in `ts`.
- The arc length is simply the length of the line segment scaled by the parameter value.

# Arguments
- `line::L`: A `LineSegment` or `VirtualLineSegment` object defining the line segment.
- `ts::AbstractArray{T,1}`: An array of parameter values, typically ranging from 0 to 1.

# Returns
- A vector of arc lengths corresponding to the parameter values in `ts`.
"""
function arc_length(line::L, ts::AbstractArray{T,1}) where {T<:Real,L<:LineSegments{T}}
    s::Vector{T} = line.length.*ts
    return s
end

"""
Compute the tangent vector along a line segment for given parameter values.

# Logic
- This function calculates the tangent vector along the line segment defined by the `line` object for each value in `ts`.
- The tangent vector is computed as the derivative of the line equation with respect to `t`, scaled by the orientation of the line segment.

# Arguments
- `line::L`: A `LineSegment` or `VirtualLineSegment` object defining the line segment.
- `ts::AbstractArray{T,1}`: An array of parameter values, typically ranging from 0 to 1.

# Returns
- A vector of tangent vectors corresponding to the parameter values in `ts`.
"""
function tangent(line::L, ts::AbstractArray{T,1}) where {T<:Real,L<:LineSegments{T}}
    let pt0 = line.cs.affine_map(line.pt0), pt1 = line.cs.affine_map(line.pt1), orient=line.orientation
        r(t) = line_eq(pt0,pt1,t)
        #ForwardDiff.derivative(r, t)
        return collect(orient*ForwardDiff.derivative(r, t) for t in ts)
    end
end

"""
Compute the signed distance from a line segment to multiple points.

# Logic
- This function calculates the signed distance of each point in `pts` from the line segment defined by the `line` object.
- The signed distance is computed using the `line_domain` function, and the result is scaled by the orientation of the line segment.

# Arguments
- `line::L`: A `LineSegment` or `VirtualLineSegment` object defining the line segment.
- `pts::AbstractArray`: An array of points, where each point is an `SVector{2,T}` representing its coordinates.

# Returns
- A vector of signed distances for each point. Positive values indicate points on one side of the line, and negative values indicate points on the other side.
"""
function domain(line::L, pts::AbstractArray) where {T<:Real,L<:LineSegments{T}}
    let pt0 = line.cs.affine_map(line.pt0) 
        pt1 = line.cs.affine_map(line.pt1)
        orientation = line.orientation
    return collect(line_domain(pt0[1],pt0[2],pt1[1],pt1[2],pt[1],pt[2])*orientation for pt in pts)
    end
end

"""
Determine if points are on a specific side of a line segment.

# Logic
- This function computes the signed distance of each point in `pts` from the line segment defined by the `line` object.
- It returns `true` for points on one side of the line segment and `false` for points on the other side, depending on the orientation of the line segment.

# Arguments
- `line::L`: A `LineSegment` or `VirtualLineSegment` object defining the line segment.
- `pts::AbstractArray`: An array of points, where each point is an `SVector{2,T}` representing its coordinates.

# Returns
- A vector of Boolean values indicating whether each point is on the specified side of the line segment.
"""
function is_inside(line::L, pts::AbstractArray) where {T<:Real,L<:LineSegments{T}}
    let pt0 = line.cs.affine_map(line.pt0) 
        pt1 = line.cs.affine_map(line.pt1)
        orientation = line.orientation
    return collect(line_domain(pt0[1],pt0[2],pt1[1],pt1[2],pt[1],pt[2])*orientation < zero(eltype(pt0)) for pt in pts)
    end
end

#=
function domain(line::L, pt::SVector{2,T}) where {T<:Real,L<:LineSegments{T}}
    let pt0 = line.cs.affine_map(line.pt0) 
        pt1 = line.cs.affine_map(line.pt1)
        orientation = line.orientation
    return line_domain(pt0[1],pt0[2],pt1[1],pt1[2],pt[1],pt[2])*orientation
    end
end
=#
