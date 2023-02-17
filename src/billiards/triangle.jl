include("../abstracttypes.jl")
include("../basis/fourierbessel/corneradapted.jl")
include("geometry.jl")
using StaticArrays #PolygonOps


function triangle_corners(angles, x0, y0, h) #x0, y0 position of gamma corner
    alpha, beta, gamma = angles
    B = SVector((-h/tan(beta+alpha))-x0, h-y0)
    A = SVector(h/tan(alpha)+B[1], -y0)
    C =  SVector(-x0, -y0)
    return SVector(A,B,C)
end


struct Triangle{T}  <: AbsBilliard where {T<:Number}
    boundary :: Vector{Any}
    length:: T
    area :: T
    corners :: SVector{3, SVector{2, T}}
    angles :: SVector{3, T}
    #domain :: Function 
end

function Triangle(gamma, chi; curve_types = [:Real, :Virtual, :Virtual] , x0=zero(gamma), y0=zero(gamma), h = one(gamma))
    alpha = (pi-gamma)/(1+chi)
    beta = alpha*chi
    angles = SVector(alpha, beta, gamma)
    #println("α=$alpha, β=$beta, γ=$gamma")
    corners = triangle_corners(angles, x0, y0, h)
    boundary = make_polygon(corners, curve_types)
    length = sum([crv.length for crv in boundary])
    area = 0.5*h*abs(corners[1][1]-corners[3][1])#PolygonOps.area(collect(zip(x,y)))
    return Triangle(boundary,length,area,corners,angles)
end


function VeechRightTriangle(kind; curve_types = [:Real, :Virtual, :Virtual] , x0=0.0, y0=0.0, h = 1.0)
    n = 4 + kind
    gamma = pi/2
    chi = (n-2)/2
    return Triangle(gamma,chi; curve_types = curve_types,  x0 = x0, y0 = y0, h=h)
end

function adapt_basis(triangle::Triangle,i)
    N = 3
    c = triangle.corners
    #println("corners $(c[3])")
    #x_axis = SVector(1.0,0.0)
    i0 = mod1(i,N)
    i1 = mod1(i+1,N)
    
    a = c[i1] - c[i0]
    rot_angle = atan(a[2],a[1])#angle(x_axis, a)
    origin = c[i0]
    cs = PolarCS(origin,rot_angle)
    return triangle.angles[i0], cs
end

#convenience function
function make_triangle_and_basis(gamma,chi; edge_i=1)
    cor = Triangle(gamma,chi).corners
    x0,y0 = cor[mod1(edge_i+2,3)]
    re = [:Virtua, :Virtual, :Virtual]
    re[edge_i] = :Real 
    tr = Triangle(gamma,chi; curve_types = re, x0 = x0, y0 =y0)
    #dim = round(Int, tr.boundary[edge_i].length*k*solver.dim_scaling_factor/(2*pi))
    basis = CornerAdaptedFourierBessel(1, adapt_basis(tr,edge_i+2)...) 
    return tr, basis 
end