using  CoordinateTransformations, Rotations, StaticArrays
using GLMakie
using Revise
include("../src/billiards/coordinatesystems.jl")
include("../src/billiards/curves.jl")

fun(t, a) = SVector(2*a*(one(t)+cos(t*(2*pi))), t*(2*pi))
cs = PolarCS(SVector(1.0,1.0),0.3) 

a = 1.0
ts = collect(range(0.0,1.0,200))
r = curve(cs,fun,ts,a)
ta = tangent_vec(cs,fun,ts, a)
na = normal_vec(cs,fun,ts, a)

f = Figure(resolution = (1000,1000));
#axis = [Axis(f[1,i]) for i in 1:3]
ax = Axis(f[1,1],xlabel=L"x", ylabel=L"y")
#lines!(ax,  r_new)
lines!(ax,  r)
arrows!(ax,getindex.(r,1),getindex.(r,2), getindex.(ta,1),getindex.(ta,2), color = :black, lengthscale = 0.1)
ax.aspect=DataAspect()
display(f)


f = Figure(resolution = (1000,1000));
#axis = [Axis(f[1,i]) for i in 1:3]
ax = Axis(f[1,1],xlabel=L"x", ylabel=L"y")
#lines!(ax,  r_new)
lines!(ax,  r)
arrows!(ax,getindex.(r,1),getindex.(r,2), getindex.(na,1),getindex.(na,2), color = :black, lengthscale = 0.1)
ax.aspect=DataAspect()
display(f)

zero(BigFloat(1.0))
#=
R = LinearMap(Angle2d(0))
T = Translation(0.5,0.5)
A = compose(T, R)
=#
#A(r[1])
#r_rotated = map(pt -> R_2d * pt, r)
