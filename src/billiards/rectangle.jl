#include("../abstracttypes.jl")
#include("../basis/fourierbessel.jl")
#include("geometry.jl")
using StaticArrays

function rectangle_corners(x0, y0, l; h = 1.0) #x0, y0 position of gamma corner
    x = [0.0, l, l, 0.0] .+ x0
    y = [0.0, 0.0, h, h] .+ y0
    return x, y
end

struct Rectangle  <: AbsBilliard
    boundary :: Vector{Any}
    length:: Float64
    area :: Float64
    corners :: Vector{SVector{2, Float64}}
    angles :: Vector{Float64}
    domain :: Function 
    function Rectangle(chi::Float64; curve_types = [:Virtual,:Real, :Real, :Virtual] , x0=0.0, y0=0.0, h=1.0)
        angles = [pi/2.0 for i in 1:4]
        #println("α=$alpha, β=$beta, γ=$gamma")
        l = chi * h
        x, y = rectangle_corners(x0, y0, l; h = h)
        corners = [SVector(x[i],y[i]) for i in eachindex(x)]
        boundary, domain = make_polygon(x, y, curve_types)
        length = sum([crv.length for crv in boundary])
        area = h*l#PolygonOps.area(collect(zip(x,y)))
        return new(boundary,length,area,corners,angles,domain)
    end
      
end

function make_rectangle_and_basis(chi; curve_types = [:Virtual,:Real, :Real, :Virtual] , x0=0.0, y0=0.0, h=1.0)
    billiard = Rectangle(chi; curve_types = curve_types , x0=x0, y0=y0, h=h)
    basis = CornerAdaptedFourierBessel(1, pi/2.0, 0.0, x0, y0) 
    return billiard, basis 
end

