#using  CoordinateTransformations, Rotations, StaticArrays
using GLMakie
using Traceur, BenchmarkTools
#using Revise
#include("../src/billiards/coordinatesystems.jl")
include("../src/billiards/geometry.jl")
include("../src/billiards/triangle.jl")
include("../src/plotting/plottingmakie.jl")
include("../src/utils/gridutils.jl")
#testing all curve types from geometry.jl

line1 = LineSegment(SVector(1.0,0.0),SVector(0.0,1.0))
origin2 = (1.0,0.0)
rot_angle2 = -0.2
line2 = LineSegment(SVector(1.0,0.0),SVector(0.0,1.0);origin=origin2, rot_angle=rot_angle2, orientation=-1)

virt_line1 = VirtualLineSegment(SVector(1.0,0.0),SVector(0.0,1.0);origin=origin2, rot_angle=rot_angle2, orientation=1)
circle1 = CircleSegment(1.0,pi/3, pi/2, 0.0, 0.0;orientation = 1)
ts = range(0.0,1.0,1000000)
pts = curve(line2,ts)
Base.format_bytes(sizeof(pts))
@code_warntype curve(circle1,ts)
@benchmark curve(circle1,ts)
@code_warntype curve(line2,ts)
@btime curve(line2,ts)
@btime curve(line1,ts)
@benchmark curve(line2,ts)

@code_warntype normal_vec(line2,ts)

@btime tangent_vec(line2,ts)


@code_warntype is_inside(line2,pts)
@btime domain(line1,pts)
@btime is_inside(line2,pts)
x_grid = range((-1.0,1.0)... ,1000)
y_grid = range((-1.0,1.0)... ,1000)

@btime is_inside(line2, x_grid, y_grid)
@btime is_inside(circle1, x_grid, y_grid)

f = Figure(resolution = (1000,1000));
#axis = [Axis(f[1,i]) for i in 1:3]
ax = Axis(f[1,1],xlabel=L"x", ylabel=L"y")
plot_curve!(ax, line1)
plot_curve!(ax, line2)
plot_curve!(ax, circle1)
display(f)

f = Figure(resolution = (1000,1000));
#axis = [Axis(f[1,i]) for i in 1:3]
ax = Axis(f[1,1],xlabel=L"x", ylabel=L"y")
plot_domain!(ax, line1)
display(f)

f = Figure(resolution = (1000,1000));
#axis = [Axis(f[1,i]) for i in 1:3]
#ax = Axis(f[1,1],xlabel=L"x", ylabel=L"y")
ax,hmap = plot_domain_fun!(f,circle1)
plot_curve!(ax, circle1)
display(f)

f = Figure(resolution = (1000,1000));
#axis = [Axis(f[1,i]) for i in 1:3]
ax = Axis(f[1,1],xlabel=L"x", ylabel=L"y")
plot_domain!(ax, circle1)
plot_curve!(ax, circle1)
display(f)

gamma = sqrt(2)/2 * pi
chi  = 2.0
tri = Triangle(gamma, chi)

f = Figure(resolution = (1000,1000));
#axis = [Axis(f[1,i]) for i in 1:3]
ax = Axis(f[1,1],xlabel=L"x", ylabel=L"y")
plot_domain!(ax,tri;dens=200.0)
plot_boundary!(ax,tri)
#plot_lattice!(ax,tri; dens= 25.0)
display(f)

x_plot, y_plot, gen = interior_grid(tri,(100,100),(-1.0,1.0),(-1.0,1.0))

X = [pt.xy  for pt in gen if pt.inside]
pts = [i for i in gen]
X = [pt.inside for pt in gen]
collect(range(0.0,1.0;step=0.01))

[true,true,true] .& [false,false,true]

