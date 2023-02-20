
include("../abstracttypes.jl")
include("../basis/fourierbessel/corneradapted.jl")
include("geometry.jl")

function make_stadium(width;radius=one(width),x0=zero(width),y0=zero(width), curve_types=[:Real,:Real,:Virtual,:Virtual],rot_angle=zero(width))
    #d(x, y, x0, y0, x1, y1) = @.((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)
    origin = SVector(x0,y0)
    circle = CircleSegment(radius,pi/4,zero(width), width, zero(width); origin=origin, rot_angle = rot_angle)
    line = (curve_types[idx0] == :Real) ? LineSegment(c0, c1;origin=origin,rot_angle=rot_angle) : VirtualLineSegment(c0,c1;origin=origin,rot_angle=rot_angle)

    line1 = (curve_types[2] == :Real) ? LineSegment(SVector(width, radius), SVector(zero(width), radius+y0);origin=origin,rot_angle=rot_angle) : VirtualLineSegment(SVector(width, radius), SVector(zero(width), radius+y0);origin=origin,rot_angle=rot_angle)
    line2 = (curve_types[3] == :Real) ? LineSegment(SVector(zero(width), radius), SVector(zero(width), zero(width));origin=origin,rot_angle=rot_angle) : VirtualLineSegment(SVector(zero(width), radius), SVector(zero(width), zero(width));origin=origin,rot_angle=rot_angle)
    line3 = (curve_types[4] == :Real) ? LineSegment(SVector(zero(width), zero(width)), SVector(width + radius, zero(width));origin=origin,rot_angle=rot_angle) : VirtualLineSegment(SVector(zero(width), zero(width)), SVector(width + radius, zero(width));origin=origin,rot_angle=rot_angle)
    boundary = [circle, line1, line2, line3]
    dm = [crv.domain for crv in boundary]
    function domain(x,y) #make better
        res = [true for i in 1:length(x)] 
        for f in dm
            res .= res .&& f(x,y)
        end
        return res
    end
    return boundary
end

struct Stadium{T}  <: AbsBilliard where {T<:Number}
    boundary :: Vector{Any}
    length:: T
    area :: T
    parameter::T
    radius::T
         
end

function Stadium(width;radius=1.0,x0=0.0,y0=0.0, curve_types=[:Real,:Real,:Virtual,:Virtual])
    boundary = make_stadium(width;radius=radius,x0=x0,y0=y0, curve_types=curve_types)
    length = sum([crv.length for crv in boundary])
    area = width*radius + (pi*radius^2)/4.0 #PolygonOps.area(collect(zip(x,y)))
    return Stadium
end 

function make_stadium_and_basis(width; curve_types = [:Real,:Real,:Virtual,:Virtual],radius=1.0,x0=zero(width),y0=zero(width), rot_angle=zero(width))
    billiard = Stadium(width;radius=radius,x0=x0,y0=y0, curve_types=curve_types)
    basis = CornerAdaptedFourierBessel(1, pi/2.0, SVector(x0,y0),rot_angle) 
    return billiard, basis 
end
