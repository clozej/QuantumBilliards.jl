#include("../abstracttypes.jl")
#include("../basis/fourierbessel/corneradapted.jl")
#include("geometry.jl")
using StaticArrays #PolygonOps


function triangle_corners(angles, x0, y0, h) #x0, y0 position of gamma corner
    alpha, beta, gamma = angles
    B = SVector((-h/tan(beta+alpha))-x0, h-y0)
    A = SVector(h/tan(alpha)+B[1], -y0)
    C =  SVector(-x0, -y0)
    return SVector(A,B,C)
end

struct Triangle{T}  <: AbsBilliard where {T<:Real}
    fundamental_boundary :: Vector
    full_boundary::Vector
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
    full_boundary = make_polygon(corners, [:Real, :Real, :Real])
    length = sum([crv.length for crv in full_boundary])
    area = 0.5*h*abs(corners[1][1]-corners[3][1])#PolygonOps.area(collect(zip(x,y)))
    return Triangle(boundary,full_boundary,length,area,corners,angles)
end

function adapt_basis(triangle::T,i) where {T<:Triangle}
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




#=
function VeechRightTriangle(kind; curve_types = [:Real, :Virtual, :Virtual] , x0=0.0, y0=0.0, h = 1.0)
    n = 4 + kind
    gamma = pi/2
    chi = (n-2)/2
    return Triangle(gamma,chi; curve_types = curve_types,  x0 = x0, y0 = y0, h=h)
end
=#