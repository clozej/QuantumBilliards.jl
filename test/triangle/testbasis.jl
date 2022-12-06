include("../../src/billiards/triangle.jl")
include("../../src/billiards/billiard.jl")
include("../../src/basis/fourierbessel.jl")
include("../../src/plotting/plottingmakie.jl")
include("../../src/plotting/matrixplotting.jl")
include("../../src/solvers/particularsolutionsmethod.jl")
include("../../src/solvers/samplers.jl")
include("../../src/solvers/decompositions.jl")
include("../../src/utils/angleutils.jl")
include("../../src/utils/benchmarkutils.jl")
using GLMakie
using Revise, BenchmarkTools
using LinearAlgebra, StaticArrays, Optim
using Statistics

x0 = 0.0 #bug in basisi for x0,y0 isnot 0.0
y0 = 0.0
h = 1.0
gamma = sqrt(2)/2 * pi
chi  = 2.0
idx = 1 #select edge 
curve_types = [:Virtua, :Virtual, :Virtual]
curve_types[idx] = :Real #set selected edge to real

triangle = Triangle(gamma,chi; curve_types = curve_types, x0 = x0, y0 = y0)
alpha, beta, _ = triangle.angles
dim = 30
#tc = triangle_corners([alpha,beta,gamma], x0, y0, h)


fb_basis = [CornerAdaptedFourierBessel(dim, adapt_basis(triangle,i+2,3)...) for i in 1:3]
#basis = CornerAdaptedFourierBessel(dim, gamma, 0.5, x0, y0)
limits = [(-2,3),(-1,2)] #[extrema(row) for row in eachrow(reduce(hcat, triangle.corners))]
k = 10.0
f = Figure(resolution = (1000,1000));
ax = Axis(f[1,1],xlabel=L"x", ylabel=L"y")
plot_basis_function!(ax,fb_basis[idx], 1, k; xlim=limits[1], ylim=limits[2])
plot_billiard!(ax, triangle)
display(f)
#c = triangle.corners
#an = triangle.angles
#angle(c[1]-c[2], c[3]-c[2])


axs = []
fig = Figure()
for (i,idx2) in enumerate(1:5)
    ax = Axis(fig[i,1])
    push!(axs, ax)
    plot_basis_function!(ax, fb_basis[idx], idx2, k, triangle.boundary[idx], gauss_legendre_nodes)
    #axislegend()
end
display(fig)


#Veech triangles
kind = 1 # Veech triangle index
idx = 1 #select edge 
curve_types = [:Virtua, :Virtual, :Virtual]
curve_types[idx] = :Real #set selected edge to real

triangle = VeechRightTriangle(kind)
alpha, beta, _ = triangle.angles
dim = 30

fb_basis = [CornerAdaptedFourierBessel(dim, adapt_basis(triangle,i+2,3)...) for i in 1:3]
#basis = CornerAdaptedFourierBessel(dim, gamma, 0.5, x0, y0)
limits = [(-2,3),(-1,2)] #[extrema(row) for row in eachrow(reduce(hcat, triangle.corners))]
k = 10.0
f = Figure(resolution = (1000,1000));
ax = Axis(f[1,1],xlabel=L"x", ylabel=L"y")
plot_basis_function!(ax,fb_basis[idx], 1, k; xlim=limits[1], ylim=limits[2])
plot_billiard!(ax, triangle)
display(f)