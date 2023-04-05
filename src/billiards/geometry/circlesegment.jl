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

function arc_length(circle::L, ts::AbstractArray{T,1}) where {T<:Real,L<:CircleSegments{T}}
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

function domain(circle::L, pts::AbstractArray) where {T<:Real,L<:CircleSegments{T}}
    let pt0 = curve(circle,zero(T)), pt1 = curve(circle,one(T)), R=circle.radius, orient = circle.orientation
        center = circle.cs.affine_map(circle.center) #move center to correct position
        return collect(circle_segment_domain(pt0, pt1, R, center, orient, pt[1], pt[2]) for pt in pts)
    end
end


function is_inside(circle::L, pts::AbstractArray) where {T<:Real,L<:CircleSegments{T}}
    let pt0 = curve(circle,zero(T)), pt1 = curve(circle,one(T)), R=circle.radius, orient = circle.orientation
        center = circle.cs.affine_map(circle.center) #move center to correct position
        return collect(circle_segment_domain(pt0, pt1, R, center, orient,pt[1],pt[2]) < zero(eltype(pt0)) for pt in pts)
    end
end
