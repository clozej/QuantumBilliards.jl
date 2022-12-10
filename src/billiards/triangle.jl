include("../abstracttypes.jl")
include("../basis/fourierbessel.jl")
include("geometry.jl")
using StaticArrays #PolygonOps


function triangle_corners(angles, x0, y0, h) #x0, y0 position of gamma corner
    alpha, beta, gamma = angles
    xB, yB = (-h/tan(beta+alpha))-x0, h-y0
    xA, yA = h/tan(alpha)+xB, 0-y0
    xC, yC =  0-x0, 0-y0
    x = [xA, xB, xC]
    y = [yA, yB, yC]
    return x, y
end


struct Triangle  <: AbsBilliard
    boundary :: Vector{Any}
    length:: Float64
    area :: Float64
    corners :: Vector{SVector{2, Float64}}
    angles :: Vector{Float64}
    domain :: Function 
    function Triangle(gamma::Float64, chi::Float64; curve_types = [:Real, :Virtual, :Virtual] , x0=0.0, y0=0.0, h = 1.0)
        alpha = (pi-gamma)/(1+chi)
        beta = alpha*chi
        angles = [alpha, beta, gamma]
        #println("α=$alpha, β=$beta, γ=$gamma")
        x, y = triangle_corners(angles, x0, y0, h)
        corners = [SVector(x[i],y[i]) for i in eachindex(x)]
        boundary, domain = make_polygon(x, y, curve_types)
        length = sum([crv.length for crv in boundary])
        area = 0.5*h*abs(x[1]-x[3])#PolygonOps.area(collect(zip(x,y)))
        return new(boundary,length,area,corners,angles,domain)
    end
      
end

function VeechRightTriangle(kind; curve_types = [:Real, :Virtual, :Virtual] , x0=0.0, y0=0.0, h = 1.0)
    n = 4 + kind
    gamma = pi/2
    chi = (n-2)/2
    return Triangle(gamma,chi; curve_types = curve_types,  x0 = x0, y0 = y0, h=h)
end

function adapt_basis(triangle::Triangle,i,N)
    c = triangle.corners
    #println("corners $(c[3])")
    #x_axis = SVector(1.0,0.0)
    i0 = mod1(i,N)
    i1 = mod1(i+1,N)
    
    a = c[i1] - c[i0]
    phi0 = atan(a[2],a[1])#angle(x_axis, a)
    x0, y0 = c[i0]
    return triangle.angles[i0], phi0, x0, y0
end

function make_triangle_and_basis(gamma,chi; edge_i=1)
    cor = Triangle(gamma,chi).corners
    x0,y0 = cor[mod1(edge_i+2,3)]
    re = [:Virtua, :Virtual, :Virtual]
    re[edge_i] = :Real 
    tr = Triangle(gamma,chi; curve_types = re, x0 = x0, y0 =y0)
    #dim = round(Int, tr.boundary[edge_i].length*k*solver.dim_scaling_factor/(2*pi))
    basis = CornerAdaptedFourierBessel(1, adapt_basis(tr,edge_i+2,3)...) 
    return tr, basis 
end