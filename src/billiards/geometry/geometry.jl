
#order of inclusions is important because linesegments need to be first
include("linesegment.jl")
include("circlesegment.jl")
include("dispersingcirclesegment.jl")
using StaticArrays,LinearAlgebra, ForwardDiff


#this part is general
"""
Compute the normalized tangent vectors for a curve at specified parameter values. The underlying `tangent` function is called based on the subtype of curve.

# Logic
- This function calculates the tangent vectors for a given curve at the parameter values specified in `ts` using the `tangent` function.
- It then normalizes each tangent vector to unit length.

# Arguments
- `curve::L`: An object of type `AbsCurve` representing the curve.
- `ts::AbstractArray{T,1}`: An array of parameter values where the tangent vectors are to be computed.

# Returns
- A vector of normalized tangent vectors corresponding to the parameter values in `ts`.
"""
function tangent_vec(curve::L, ts::AbstractArray{T,1}) where {T<:Real,L<:AbsCurve}
    ta = tangent(curve, ts)
    return collect(ti/norm(ti) for ti in ta)
end

"""
Compute the normal vectors for a curve at specified parameter values.

# Logic
- This function calculates the normalized tangent vectors for the given curve using the `tangent_vec` function.
- It then computes the normal vector at each point by rotating the tangent vector 90 degrees counterclockwise.

# Arguments
- `curve::L`: An object of type `AbsCurve` representing the curve.
- `ts::AbstractArray{T,1}`: An array of parameter values where the normal vectors are to be computed.

# Returns
- A vector of normal vectors corresponding to the parameter values in `ts`.
"""
function normal_vec(curve::L, ts::AbstractArray{T,1}) where {T<:Real,L<:AbsCurve}
    ta = tangent_vec(curve, ts)
    return [SVector(ti[2], -ti[1]) for ti in ta]
end

"""
Creates a polygon billiard from LineSegments or VirtualLineSegments. The later are included when the user wants to construct a desymmetrized a polygon

# Arguments
- `corners::Vector{SVector{2,T}}`: A vector of corner points representing the vertices of the polygon.
- `curve_types::Vector{Symbol}`: A vector represeting the either :Real or :Imag line segments that construct the polygon billiard
- `origin::Tuple{T,T}`: The origin of the coordinate system, defaults to the (0,0) of the type T
- `rot_angle::T`: The rotation angle in radians, defaults to 0 of the type T

# Returns
- A vector of `LineSegment` or `VirtualLineSegment` representing the boundary of the polygon billiard.
"""
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



#domain fuction has node (zero) at the edge of the curve

#= worse version
function domain2(line::LineSegment{T}, pts::AbstractArray) where {T<:Real}
    let line = line
    return collect(domain(line,pt) for pt in pts)
    end
end
=#

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


function domain(line::L, x_grid::AbstractArray, y_grid::AbstractArray) where {T<:Real,L<:LineSegments{T}}
    let pt0 = line.cs.affine_map(line.pt0) 
        pt1 = line.cs.affine_map(line.pt1)
        orientation = line.orientation
    return collect(line_domain(pt0[1],pt0[2],pt1[1],pt1[2],x,y)*orientation for y in y_grid for x in x_grid)
    end
end

function is_inside(line::L, x_grid::AbstractArray, y_grid::AbstractArray) where {T<:Real,L<:LineSegments{T}}
    let pt0 = line.cs.affine_map(line.pt0) 
        pt1 = line.cs.affine_map(line.pt1)
        orientation = line.orientation
    return collect(line_domain(pt0[1],pt0[2],pt1[1],pt1[2],x,y)*orientation < zero(eltype(x_grid)) for y in y_grid for x in x_grid)
    end
end

function domain(circle::L, x_grid::AbstractArray, y_grid::AbstractArray) where {T<:Real,L<:CircleSegments{T}}
    let pt0 = curve(circle,zero(T)), pt1 = curve(circle,one(T)), R=circle.radius, orient = circle.orientation
        center = circle.cs.affine_map(circle.center) #move center to correct position
        return collect(circle_segment_domain(pt0, pt1, R, center, orient,x,y) for y in y_grid for x in x_grid)
    end
end

function is_inside(circle::L, x_grid::AbstractArray, y_grid::AbstractArray) where {T<:Real,L<:CircleSegments{T}}
    let pt0 = curve(circle,zero(T)), pt1 = curve(circle,one(T)), R=circle.radius, orient = circle.orientation
        center = circle.cs.affine_map(circle.center) #move center to correct position
        return collect(circle_segment_domain(pt0, pt1, R, center, orient,x,y) < zero(eltype(x_grid)) for y in y_grid for x in x_grid)
    end
end
=#