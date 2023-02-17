#using  CoordinateTransformations, Rotations, StaticArrays
using GLMakie
#using Revise
#include("../src/billiards/coordinatesystems.jl")
#include("../src/billiards/geometry.jl")
#include("../src/billiards/triangle.jl")
include("../src/basis/fourierbessel/corneradapted.jl")
include("../src/plotting/plottingmakie.jl")
include("../src/utils/gridutils.jl")
#testing all curve types from geometry.jl
using StaticArrays

dim::Int = 100
corner_angle = Float32(sqrt(2)/4*pi)
origin = SVector(Float32.([-0.5,0.5])...) 
rot_angle = Float32(0.0)# -pi/4.0
basis = CornerAdaptedFourierBessel(dim, corner_angle, origin, rot_angle)

x_plot, y_plot, gen = make_grid((1000,1000),Float32.([-1.5,1.5]), Float32.([-1.5,1.5]))
Z = [basis_fun(basis, 1, Float32(1000.0), pt) for pt in gen]

f = Figure(resolution = (1000,1000));
#axis = [Axis(f[1,i]) for i in 1:3]
plot_heatmap_balaced!(f,x_plot, y_plot, Z)
display(f)
