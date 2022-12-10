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

#Strong-Mixing
gamma = pi/2*(1+sqrt(5))/2
chi  = (1+sqrt(2))/2

#Weak-Mixing
gamma = sqrt(2)/2 * pi
chi  = 2.0

#Veech
kind = 5
n = 4 + kind
gamma = pi/2
chi = (n-2)/2


function plot_angle!(ax, alpha, phi0, x0, y0)
    r = 0.2
    x1 = cos(phi0) + x0
    x2 = cos(alpha + phi0) + x0
    y1 = sin(phi0) + y0
    y2 = sin(alpha + phi0) + y0
    lines!(ax, [x1,x0,x2],[y1,y0,y2], color=:Red)
end

f = Figure(resolution = (1000,2000));
triangle = Triangle(gamma,chi)
fb_basis = [CornerAdaptedFourierBessel(10, adapt_basis(triangle,i+2,3)...) for i in 1:3]

k = 10.0
limits = [(-2,3),(-1,2)]
axis = [Axis(f[i,1]) for i in 1:3]
for idx in 1:3
    re = [:Virtua, :Virtual, :Virtual]
    re[idx] = :Real 
    tr = Triangle(gamma,chi; curve_types = re)
    alpha,phi0,x0,y0 = adapt_basis(tr,idx+2,3)
    println("angle=$alpha,  phi=$phi0")
    basis = CornerAdaptedFourierBessel(10, alpha,phi0,x0,y0)#fb_basis[idx]
    #println(tr.angles)
    #println(pi/basis.nu)
    println(basis)
    plot_basis_function!(axis[idx],basis, 1, k; xlim=limits[1], ylim=limits[2])
    plot_angle!(axis[idx], adapt_basis(triangle,idx+2,3)...)
    plot_billiard!(axis[idx], tr)
end
display(f)


function plot_triangle_sv_test(gamma,chi,ks;b = 40.0,b_int=40.0)
    f = Figure(resolution = (1000,500));
    axis = [Axis(f[1,i]) for i in 1:3]
    ax = Axis(f[2:3,:],xlabel=L"k", ylabel=L"sv(k)")

    k_min, k_max = extrema(ks)
    triangle = Triangle(gamma,chi; curve_types = [:Real, :Real, :Real])
    k = k_max
    dk = 0.5
    solver = ParticularSolutionsMethod(b)
    solver_dist = ParticularSolutionsMethod([b for i in 1:3],b_int) #distributed solver
    pts_vec = evaluate_points(solver_dist, triangle, gauss_legendre_nodes, k)
    basis_vec =[] 
    #println(basis_vec)
    pad = (-0.2,0.2)
    limits = [extrema(row) for row in eachrow(reduce(hcat, triangle.corners))]
    xlim, ylim = limits
    
    #limits = [(-2,3),(-1,2)]
    for i in 1:3
        println("")
        println("Computing single edge $i")
        re = [:Virtua, :Virtual, :Virtual]
        re[i] = :Real 
        tr = Triangle(gamma,chi; curve_types = re)
        dim = round(Int, triangle.boundary[i].length*k*solver_dist.dim_scaling_factor[i]/(2*pi))
        al,phi0,x0,y0 = adapt_basis(tr,i+2,3)
        #println("angle=$alpha,  phi=$phi0")
        basis = CornerAdaptedFourierBessel(dim, al,phi0,x0,y0) 
        push!(basis_vec,basis)
        #println(basis)
        plot_basis_function!(axis[i],basis, 1, k; xlim=xlim.+ pad, ylim=ylim .+ pad)
        plot_angle!(axis[i], al,phi0,x0,y0)
        plot_billiard!(axis[i], tr)
        benchmark_solver(solver, basis, tr, gauss_legendre_nodes, k, dk)
        #println(basis)
        res = k_sweep(solver,basis, pts_vec[i], ks)
        lines!(ax, ks, log10.(res), label = "edge $i")
        

    end
    
    println("")
    println("Computing distributed method")
    benchmark_solver(solver_dist, basis_vec, triangle, gauss_legendre_nodes, k, dk)
    res = k_sweep(solver_dist,basis_vec, pts_vec, ks)
    lines!(ax, ks, log10.(res), label = "distributed")
    alpha, beta, gamma = triangle.angles
    Label(f[:,:, Top()], L"$\alpha=%$(alpha/pi) \pi$, $\beta=%$(beta/pi) \pi$, $\gamma=%$(gamma/pi)\pi$ , $\chi=%$chi$", valign = :center)
    leg = Legend(f[2:3, :],ax, tellheight=false, tellwidth=false,
    margin = (10, 10, 10, 10),
    halign = :right, valign = :bottom)
    return f
end


#Strong-Mixing
gamma = pi/2*(1+sqrt(5))/2
chi  = (1+sqrt(2))/2

#Weak-Mixing
gamma = sqrt(2)/2 * pi
chi  = 2.0

#Veech
kind = 7
n = 4 + kind
gamma = pi/2
chi = (n-2)/2


ks = collect(LinRange(6.0,15.0,1000))
#ks = collect(LinRange(100.0,100.1,500))
f = plot_triangle_sv_test(gamma,chi,ks;b = 40.0,b_int=40.0)
display(f)

save("testsweep2.png",f,pt_per_unit = 2)