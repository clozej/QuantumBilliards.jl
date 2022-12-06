
include("../abstracttypes.jl")
include("../basis/fourierbessel.jl")
include("geometry.jl")

function make_stadium(eps;R=1,x0=0.0,y0=0.0, curve_types=[:Real,:Real,:Virtual,:Virtual])
    #d(x, y, x0, y0, x1, y1) = @.((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)

    circle = CircleSegment(R, x0+eps, y0, pi/2.0, 0.0)

    line1 = (curve_types[2] == :Real) ? LineSegment(eps+x0, R+y0, x0, R+y0) : VirtualLineSegment(eps+x0, R+y0, x0, R+y0)
    line2 = (curve_types[3] == :Real) ? LineSegment(x0, R+y0, x0, y0) : VirtualLineSegment(x0, R+y0, x0, y0)
    line3 = (curve_types[4] == :Real) ? LineSegment(x0, y0, eps + R + x0, y0) : VirtualLineSegment(x0, y0, eps + R + x0, y0)
    boundary = [circle, line1, line2, line3]
    dm = [crv.domain for crv in boundary]
    function domain(x,y) #make better
        res = [true for i in 1:length(x)] 
        for f in dm
            res .= res .&& f(x,y)
        end
        return res
    end
    return boundary, domain
end

struct Stadium  <: AbsBilliard
    boundary :: Vector{Any}
    length:: Float64
    area :: Float64
    domain :: Function 
    function Stadium(eps;R=1.0,x0=0.0,y0=0.0, curve_types=[:Real,:Real,:Virtual,:Virtual])
        boundary, domain = make_stadium(eps;R=1,x0=x0,y0=y0, curve_types=curve_types)
        length = sum([crv.length for crv in boundary])
        area = eps*R + (pi*R^2)/4.0 #PolygonOps.area(collect(zip(x,y)))
        return new(boundary,length,area,domain)
    end      
end

function make_stadium_and_basis(eps; curve_types = [:Real,:Real,:Virtual,:Virtual],R=1.0,x0=0.0,y0=0.0)
    billiard = Stadium(eps;R=1,x0=x0,y0=y0, curve_types=curve_types)
    basis = CornerAdaptedFourierBessel(1, pi/2.0, 0.0, x0, y0) 
    return billiard, basis 
end
