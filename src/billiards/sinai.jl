function make_quarter_sinai(angle1, angle2; x0=zero(angle1),y0=zero(angle1),rot_angle=zero(angle1))
    origin = SVector(x0,y0)
    type = typeof(angle1)
    pos_x = cot(angle1)
    R1 = 1.0/sin(angle1)
    pos_y = cot(angle2)
    R2 = 1.0/sin(angle2)

    circle1 = DispersingCircleSegment(R1,angle1, zero(type), 1.0+pos_x, zero(type);origin=origin,rot_angle=rot_angle)
    circle2 = DispersingCircleSegment(R2,angle2, pi/2+angle2, zero(angle1), 1.0+pos_y;origin=origin,rot_angle=rot_angle)    
    corners = [SVector(one(type)+pos_x-R1, zero(type)), SVector(one(type), one(type)),SVector(zero(type),one(type)+pos_y-R2), SVector(zero(type), zero(type)) ]
    
    line1 = VirtualLineSegment(corners[3],corners[4];origin=origin,rot_angle=rot_angle)
    line2 = VirtualLineSegment(corners[4],corners[1];origin=origin,rot_angle=rot_angle)
    boundary = [circle1, circle2, line1, line2]
    return boundary, corners
end

function make_full_sinai(angle1, angle2; x0=zero(angle1),y0=zero(angle1),rot_angle=zero(angle1))
    origin = SVector(x0,y0)
    type = typeof(angle1)
    pos_x = cot(angle1)
    R1 = 1.0/sin(angle1)
    pos_y = cot(angle2)
    R2 = 1.0/sin(angle2)

    circle1 = DispersingCircleSegment(R1,2.0*angle1, angle1, 1.0+pos_x, zero(type);origin=origin,rot_angle=rot_angle)
    circle2 = DispersingCircleSegment(R2,2.0*angle2, pi/2+angle2, zero(angle1), 1.0+pos_y;origin=origin,rot_angle=rot_angle)    
    circle3 = DispersingCircleSegment(R1,2.0*angle1, pi+angle1, -(1.0+pos_x), zero(type);origin=origin,rot_angle=rot_angle)
    circle4 = DispersingCircleSegment(R2,2.0*angle2, 3*pi/2+angle2, zero(angle1), -(1.0+pos_y);origin=origin,rot_angle=rot_angle)    
   
    corners = [SVector(one(type), one(type)), SVector(-one(type), one(type)), SVector(-one(type), -one(type)), SVector(one(type), -one(type))]
    boundary = [circle1, circle2, circle3, circle4]
    return boundary, corners
end

struct Sinai{T}  <: AbsBilliard where {T<:Real}
    fundamental_boundary::Vector
    full_boundary::Vector
    length::T
    area::T
    angle1::T
    angle2::T
    corners::Vector{SVector{2,T}}
end

function Sinai(angle1,angle2; full_boundary=false, x0=0.0,y0=0.0,rot_angle=0.0)
    circ_seg_area1 = 1.0/(sin(angle1))^2*(2*angle1-sin(2.0*angle1))
    circ_seg_area2 = 1.0/(sin(angle2))^2*(2*angle2-sin(2.0*angle2))

    full_boundary, corners = make_full_sinai(angle1,angle2;x0=x0,y0=y0,rot_angle=rot_angle)
    area = 4.0 - (circ_seg_area1 + circ_seg_area2)
    boundary, _ = make_quarter_sinai(angle1,angle2;x0=x0,y0=y0,rot_angle=rot_angle)

    length = sum([crv.length for crv in full_boundary])
    #PolygonOps.area(collect(zip(x,y)))
    if full_boundary
        Sinai(full_boundary,full_boundary,length,area,angle1,angle2,corners)
    else
        return Sinai(boundary,full_boundary,length,area,angle1,angle2,corners)
    end
end