#using  CoordinateTransformations, Rotations, StaticArrays
using GLMakie
#using Revise
#include("../src/billiards/coordinatesystems.jl")
#include("../src/billiards/geometry.jl")
include("../src/billiards/triangle.jl")
include("../src/solvers/acceleratedmethods/acceleratedmethods.jl")
include("../src/solvers/acceleratedmethods/scalingmethod.jl")
include("../src/plotting/plottingmakie.jl")
#include("../src/utils/gridutils.jl")
#testing all curve types from geometry.jl

gamma = sqrt(2)/2 * pi
chi  = 2.0

tri, basis = make_triangle_and_basis(gamma, chi)

d = 5.0
b = 5.0
acc_solver = ScalingMethod(d,b)

k0 = 300.8
dk = 0.05
k, ten = solve_wavenumber(acc_solver,basis,tri,k0,dk)
ks, tens = solve_spectrum(acc_solver,basis,tri,k0,dk)

