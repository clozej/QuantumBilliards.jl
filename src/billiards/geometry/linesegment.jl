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

function domain(line::L, pts::AbstractArray) where {T<:Real,L<:LineSegments{T}}
    let pt0 = line.cs.affine_map(line.pt0) 
        pt1 = line.cs.affine_map(line.pt1)
        orientation = line.orientation
    return collect(line_domain(pt0[1],pt0[2],pt1[1],pt1[2],pt[1],pt[2])*orientation for pt in pts)
    end
end



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
