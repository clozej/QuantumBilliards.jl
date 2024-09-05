#include("../abstracttypes.jl")
#include("../basis/fourierbessel/corneradapted.jl")
#include("geometry.jl")
using StaticArrays 

"""
Calculate the corner points of a triangle given its angles, base height, and the position of one corner.

# Logic
- The function starts by calculating the coordinates of point B, using trigonometry to project from the given gamma angle and the base height `h`.
- It then calculates the coordinates of point A using the previously calculated coordinates of point B.
- Finally, the coordinates of point C are taken as the given origin point `(-x0, -y0)`.
- The corners are returned in the order A, B, C, which corresponds to the order of the angles in the input vector `angles`.

# Arguments
- `angles::SVector{3, T}`: A static vector containing the three interior angles of the triangle.
- `x0::T`: The x-coordinate of the position of the gamma corner (third vertex).
- `y0::T`: The y-coordinate of the position of the gamma corner (third vertex).
- `h::T`: The height of the triangle from the gamma corner to the base.

# Returns
- `SVector{3, SVector{2, T}}`: A static vector containing the positions of the three corners of the triangle in the order: A, B, C.
"""
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

"""
Constructs a triangular billiard table given one interior angle, the ratio of the other two angles, and optional position and height.

# Logic
- The constructor first calculates the angles `alpha` and `beta` using the provided `gamma` angle and the `chi` ratio.
- These angles are used to calculate the positions of the corners of the triangle using the `triangle_corners` function.
- The `make_polygon` function is used to create the boundary segments of the triangle. The `curve_types` argument specifies which segments are real and which are virtual.
- The `full_boundary` is constructed using by using all the segments as real.
- The `length` is calculated as the sum of the lengths of the segments in the `full_boundary`.
- The `area` is calculated as half the product of the height and the absolute difference between the x-coordinates of the corners A and C. The standard k1*k2/2 rule

# Arguments
- `gamma::T`: The third interior angle of the triangle.
- `chi::T`: The ratio of the angles.
- `curve_types::Vector{Symbol}=:[:Real, :Virtual, :Virtual]`: A vector specifying the type of each boundary segment. 
- `x0::T=zero(gamma)`: The x-coordinate of the gamma corner.
- `y0::T=zero(gamma)`: The y-coordinate of the gamma corner.
- `h::T=one(gamma)`: The height of the triangle from the gamma corner to the base.

# Returns
- An instance of the `Triangle` struct.
"""
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

"""
Adapts the basis for a given corner of the triangle, returning the associated angle and coordinate system.

# Logic
- The function identifies the two corners that define the edge starting at the `i`-th corner.
- It calculates the vector `a` representing the edge and uses it to determine the rotation angle relative to the x-axis.
- The `origin` for the polar coordinate system is set as the `i`-th corner.
- The function then returns the angle at the `i`-th corner and the corresponding `PolarCS` coordinate system, which is aligned with the edge starting at the corner.


# Arguments
- `triangle::Triangle{T}`: The `Triangle` instance for which the basis needs to be adapted.
- `i::Int`: The index of the corner for which to adapt the basis. This is a corner index

# Returns
- `Tuple{T, PolarCS{T}}`: A tuple containing the angle at the corner and the corresponding polar coordinate system.
"""
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