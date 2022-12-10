include("../abstracttypes.jl")
include("curves.jl")

line(t, x0, x1) = (x1-x0) * t + x0 

circle_x(t, R, x0, y0, angle, shift) = x0 .+ (R .* cos.(t .* angle .+ shift))
circle_y(t, R, x0, y0, angle, shift) = y0 .+ (R .* sin.(t .* angle .+ shift))


struct LineSegment  <: AbsRealCurve
    r::Function
    n::Function
    s::Function
    ds::Function
    orientation::Int
    length::Real
    domain::Function
    function LineSegment(x0, y0, x1, y1; orientation = 1)
        line_x(t; x0 = x0, x1= x1, kwargs...) = line(t, x0, x1)
        line_y(t; y0 = y0, y1= y1, kwargs...) = line(t, y0, y1)
        r = t -> curve(t, line_x, line_y; x0 = x0, x1= x1, y0 = y0, y1= y1 )
        n = t -> orientation .* normal_vec(t, line_x, line_y; x0 = x0, x1= x1, y0 = y0, y1= y1 )
        L = sqrt((x1 - x0)^2 + (y1 -y0)^2)
        s = t -> L .* t
        ds = t -> tangent_length(t, line_x, line_y; x0 = x0, x1= x1, y0 = y0, y1= y1 )
        d(x, y) = @.((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)*orientation .<0.0
        return new(r,n,s,ds,orientation,L,d)
    end
end

struct VirtualLineSegment  <: AbsVirtualCurve
    r::Function
    n::Function
    s::Function
    ds::Function
    orientation::Int
    length::Real
    domain::Function
    function VirtualLineSegment(x0, y0, x1, y1; orientation = 1)
        line_x(t; x0 = x0, x1= x1, kwargs...) = line(t, x0, x1)
        line_y(t; y0 = y0, y1= y1, kwargs...) = line(t, y0, y1)
        r = t -> curve(t, line_x, line_y; x0 = x0, x1= x1, y0 = y0, y1= y1 )
        n = t -> orientation .* normal_vec(t, line_x, line_y; x0 = x0, x1= x1, y0 = y0, y1= y1 )
        L = sqrt((x1 - x0)^2 + (y1 -y0)^2)
        s = t -> L .* t
        ds = t -> tangent_length(t, line_x, line_y; x0 = x0, x1= x1, y0 = y0, y1= y1 )
        d(x, y) = @.((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)*orientation .<0.0
        return new(r,n,s,ds,orientation,L,d)
    end
end

struct CircleSegment  <: AbsRealCurve
    r::Function
    n::Function
    s::Function
    ds::Function
    orientation::Int
    length::Real
    domain::Function
    function CircleSegment(R, x0, y0, angle, shift; orientation = 1)
        r = t -> curve(t, circle_x, circle_y,  R, x0, y0, angle, shift)
        n = t -> orientation .* normal_vec(t, circle_x, circle_y,  R, x0, y0, angle, shift)
        L = R*angle 
        s = t -> L .* t
        ds = t -> tangent_length(t, circle_x, circle_y,  R, x0, y0, angle, shift)
        c_x0, c_y0 = r(0.0)
        c_x1, c_y1 = r(1.0)
        # line domain used to cut off circle and extend domain
        d1(x, y) = @. ((c_y1-c_y0)*x-(c_x1-c_x0)*y+c_x1*c_y0-c_y1*c_x0)*orientation < 0.0
        d2(x, y) = @. hypot(y-y0,x-x0)*orientation < R
        
        d(x,y) = @. d1(x,y) || d2(x,y)
        return new(r,n,s,ds,orientation,L,d)
    end
end


function make_polygon(x, y, curve_types)
    N = length(x)
    boundary = []
    dm = []
    circular_idx(i) = mod1(i,N)
    #d(x, y, x0, y0, x1, y1) = @.((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)

    for i in 1:N
        idx0 = circular_idx(i)
        idx1 = circular_idx(i+1)
        x0, y0, x1, y1 = x[idx0], y[idx0], x[idx1], y[idx1]
        line = (curve_types[idx0] == :Real) ? LineSegment(x0, y0, x1, y1) : VirtualLineSegment(x0, y0, x1, y1)
        push!(boundary, line)
        push!(dm, line.domain) 
    end
    
    function domain(x,y) #make better
        res = [true for i in 1:length(x)] 
        for f in dm
            res .= res .&& f(x,y)
        end
        return res
    end
    
    return boundary, domain
end

#=
x0 = 0.0
x1 = 2.0
y0 = 0.0
y1 = 0.0

line1 = LineSegment(x0, y0, x1, y1)

t = collect(LinRange(0,1.0,10000))
x,y = line1.r(t)
ds = line1.s(t)
=#