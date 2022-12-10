include("../src/billiards/triangle.jl")
include("../src/billiards/billiard.jl")
include("../src/basis/fourierbessel.jl")
include("../src/plotting/plottingmakie.jl")
include("../src/plotting/matrixplotting.jl")
include("../src/solvers/particularsolutionsmethod.jl")
include("../src/solvers/samplers.jl")
include("../src/solvers/decompositions.jl")
include("../src/utils/angleutils.jl")
include("../src/utils/benchmarkutils.jl")
using GLMakie
using Revise, BenchmarkTools
using LinearAlgebra, StaticArrays, Optim
using Statistics

gamma = sqrt(2)/2 * pi
chi  = 2.0
triangle = Triangle(gamma,chi; curve_types = [:Real, :Real, :Real])
alpha, beta, _ = triangle.angles
dim = 30

function plot_triangle_sv!(ax,gamma,chi,ks;b = 40.0,b_int=40.0)
    k_min, k_max = extrema(ks)
    triangle = Triangle(gamma,chi; curve_types = [:Real, :Real, :Real])
    k = k_max
    dk = 0.5
    solver = ParticularSolutionsMethod(b)
    solver_dist = ParticularSolutionsMethod([b for i in 1:3],b_int) #distributed solver
    pts_vec = evaluate_points(solver_dist, triangle, gauss_legendre_nodes, k)
    basis_vec = [CornerAdaptedFourierBessel(round(Int, triangle.boundary[i].length*k*solver_dist.dim_scaling_factor[i]/(2*pi)), adapt_basis(triangle,i+2,3)...) for i in 1:3]

    for i in 1:3
        println("Computing single edge $i")
        re = [:Virtua, :Virtual, :Virtual]
        re[i] = :Real 
        tr = Triangle(gamma,chi; curve_types = re)
        benchmark_solver(solver, basis_vec[i], tr, gauss_legendre_nodes, k, dk)
        res = k_sweep(solver,basis_vec[i], pts_vec[i], ks)
        lines!(ax, ks, log10.(res), label = "edge $i")
    end
    println("Computing distributed method")
    benchmark_solver(solver_dist, basis_vec, triangle, gauss_legendre_nodes, k, dk)
    res = k_sweep(solver_dist,basis_vec, pts_vec, ks)
    lines!(ax, ks, log10.(res), label = "distributed")
end

#ks = collect(LinRange(6.96,6.97,500))
#ks = collect(LinRange(6.5,15.0,500))
ks = collect(LinRange(100.0,100.1,500))
f = Figure(resolution = (1000,1000));
ax = Axis(f[1,1],xlabel=L"k", ylabel=L"sv(k)")

plot_triangle_sv!(ax,gamma,chi,ks)

Label(f[1,1, Top()], "γ=$(gamma/pi) π, χ=$chi", valign = :bottom)
leg = Legend(f[1, 1],ax, tellheight=false, tellwidth=false,
margin = (10, 10, 10, 10),
halign = :right, valign = :bottom)
display(f)

struct TestResult 
    k::Float64
    t::Float64
    dec_time::Float64
    mat_time::Float64
end

function compute_tests(gamma, chi, b_dim, b_pts, b_int, sampler, k, dk)
    triangle = Triangle(gamma,chi; curve_types = [:Real, :Real, :Real])
    solver = ParticularSolutionsMethod(b_dim, b_pts, b_int)
    solver_dist = DistributedParticularSolutionsMethod([b_dim for i in 1:3],
    [b_pts for i in 1:3], b_int, 1e-15) #distributed solver
    #pts_vec = evaluate_points(solver_dist, triangle, gauss_legendre_nodes, k)
    basis_vec = [CornerAdaptedFourierBessel(round(Int, triangle.boundary[i].length*k*solver_dist.dim_scaling_factor[i]/(2*pi)), adapt_basis(triangle,i+2,3)...) for i in 1:3]
    res = []
    for i in 1:3
        #println("Computing single edge $i")
        re = [:Virtua, :Virtual, :Virtual]
        re[i] = :Real 
        tr = Triangle(gamma,chi; curve_types = re)
        k0, t0, decomp_time, mat_time = benchmark_solver(solver, basis_vec[i], tr, sampler, k, dk; print_info=false)
        r = TestResult(k0, t0, decomp_time, mat_time)
        push!(res,r)
    end
    #println("Computing distributed method")
    k0, t0, decomp_time, mat_time = benchmark_solver(solver_dist, basis_vec, triangle, sampler, k, dk; print_info=false)
    r = TestResult(k0, t0, decomp_time, mat_time)
    push!(res,r)
    return res
end


function plot_tests(gamma, chi, sampler, k, dk)
    

    bs = collect(Float64,1:0.5:10)
    k_vals = Array{Float64}(undef,4,length(bs))
    sv = Array{Float64}(undef,4,length(bs))
    tim_mat = Array{Float64}(undef,4,length(bs))
    tim_decomp = Array{Float64}(undef,4,length(bs))
    
    for (i,b_dim) in enumerate(bs)
        println(i)
        b_pts = 3*b_dim
        b_int = 20.0
        res = compute_tests(gamma, chi, b_dim, b_pts, b_int, sampler, k, dk)
        println(res[1].k)
        println(res[1].t)
        k_vals[:,i] = [r.k for r in res]
        sv[:,i] = [r.t for r in res]
        tim_mat[:,i] = [r.mat_time for r in res]
        tim_decomp[:,i] = [r.dec_time for r in res]
        
    end

    f = Figure(resolution = (1000,1000));
    ax1 = Axis(f[1,1],xlabel=L"b", ylabel=L"sv")
    ax2 = Axis(f[1,2],xlabel=L"b", ylabel=L"k")
    ax3 = Axis(f[2,1],xlabel=L"\log b", ylabel=L"matrix time")
    ax4 = Axis(f[2,2],xlabel=L"\log b", ylabel=L"decomp time")
    for i in 1:4
        lines!(ax1, bs, log10.(sv[i,:]),label = "edge $i")
        lines!(ax2, bs, k_vals[i,:],label = "edge $i")
        lines!(ax3, log10.(bs), log10.(tim_mat[i,:]),label = "edge $i")
        lines!(ax4, log10.(bs), log10.(tim_decomp[i,:]),label = "edge $i")
    end
  
    Label(f[1,1, Top()], "γ=$(gamma/pi) π, χ=$chi", valign = :bottom)
    #leg = Legend(f[1, 1],ax1, tellheight=false, tellwidth=false,
    #margin = (10, 10, 10, 10),
    #halign = :right, valign = :top)
    
    display(f)
end
        
k, dk = 6.95, 0.05
plot_tests(gamma, chi, gauss_legendre_nodes, k, dk)
#=
function plot_triangle_sv!(ax,gamma,chi,solver,ks)
    k_min, k_max = extrema(ks)
    triangle = Triangle(gamma,chi; curve_types = [:Real, :Real, :Real])
    k = k_max
    L = real_length(triangle)
    M = round(Int, L*k*solver.int_pts_scaling_factor/(2*pi))
    #f = Figure(resolution = (3000,1000))
    x_int, y_int = random_interior_points(triangle,M)
6.9654290
6.9654290
    function solve(k)
        B3 = []
        B3_int = []
        res = Float64[]
        for (i,curve) in enumerate(triangle.boundary)
            L = curve.length
            dim = round(Int, L*k*solver.dim_scaling_factor/(2*pi))
            fb_basis = CornerAdaptedFourierBessel(dim, adapt_basis(triangle,i+2,3)...)
            x, y = evaluate_points(solver, curve, gauss_legendre_nodes, k) 
            pts = PSMPoints(x,y,x_int,y_int)
            B, B_int = construct_matrices(solver,fb_basis, pts, k)
            push!(B3, B)
            push!(B3_int, B_int)
            solution = svdvals(B, B_int)
            push!(res, minimum(solution))
        end
        #B = directsum(directsum(B3[1],B3[2]),B3[3])
        B = reduce(directsum, B3)
        B_int = reduce(hcat, B3_int)
        solution = svdvals(B, B_int)
        push!(res, minimum(solution))
        return res
    end
    result = Array{Float64}(undef,4,length(ks))
    #println(size(result))
    for (i,k) in enumerate(ks)
        #sol = solve(k)
        println(k)
        #println(sol)
        result[:,i] .= solve(k)
        #println(size(result))
    end
    for i in 1:4
        lines!(ax, ks, log10.(result[i,:]))
    end
end
=#


ks1 = collect(LinRange(6.96,6.97,1000))
ks2 = collect(LinRange(6.0,8.0,1000))
ks = sort(append!(ks1,ks2))

#ks = collect(LinRange(6.9654285,6.9654295,100))
#res1 = k_sweep(solver,triangle,fb_basis,chebyshev_nodes, ks)
solver_configs = [ParticularSolutionsMethod(b, b, b)  for b in [40.0]]
solver_distributed = ParticularSolutionsMethod([40.0 for i in 1:3], 40.0)

k = 6.965
dk = 0.01
benchmark_solver(solver_configs[1], fb_basis[1], triangle, gauss_legendre_nodes, k, dk)
benchmark_solver(solver_distributed, fb_basis, triangle, gauss_legendre_nodes, k, dk)




function tests(solvers, sampler)
    constr_times = Float64[]
    GESVD_times = Float64[]
    GESVD_sol = Float64[]
    #GESVD_k = Float64[]
    svd_times = Float64[]
    svd_sol = Float64[]

    gamma = sqrt(2)/2 * pi
    chi  = 2.0
    k_min, k_max = 6.96542, 6.96543
    triangle = Triangle(gamma,chi)
    k = k_max
    L = real_length(triangle)

    for solver in solvers
        dim = round(Int, L*solver.dim_scaling_factor/(2*pi))
        println(dim)
        basis =  CornerAdaptedFourierBessel(dim, gamma, 0.0, 0.0, 0.0)
        pts = evaluate_points(solver, triangle, sampler, k) 
        B, B_int = construct_matrices(solver,basis, pts, k)
        push!(constr_times, mean([@elapsed B, B_int for i in 1:10]))
        push!(GESVD_times, mean([@elapsed GESVDVALS(B, B_int; eps=solver.eps) for i in 1:10]))
        push!(svd_times, mean([@elapsed svdvals(B, B_int) for i in 1:10]))

        function f1(k)
            B, B_int = construct_matrices(solver,basis, pts, k)
            solution = GESVDVALS(B, B_int; eps=solver.eps)
            return minimum(solution)
        end

        function f2(k)
            B, B_int = construct_matrices(solver,basis, pts, k)
            solution = svdvals(B, B_int)
            return minimum(solution)
        end
        sol1 = optimize(f1, k_min, k_max)
        sol2 = optimize(f2, k_min, k_max)
        push!(GESVD_sol,sol1.minimum)
        push!(svd_sol,sol2.minimum)
    end
    return [constr_times GESVD_times svd_times GESVD_sol svd_sol]
end

solver_configs = [ParticularSolutionsMethod(b,b,b, 1e-8)  for b in 1:2:100]

tests(solver_configs[1:10], gauss_legendre_nodes)

function plot_benchmarks(bs)
    solver_configs = [ParticularSolutionsMethod(b,2*b,2*b, 1e-8)  for b in bs]
    f = Figure(resolution = (1000,1000));
    ax1 = Axis(f[1,1],xlabel=L"b", ylabel=L"sv")
    #ax2 = Axis(f[2,1],xlabel=L"b", ylabel=L"time")
    
    test = tests(solver_configs, gauss_legendre_nodes)
    lines!(ax1, bs, log10.(test[:,4]), color=:blue)
    lines!(ax1, bs, log10.(test[:,5]), color=:red)

    test = tests(solver_configs, chebyshev_nodes)
    lines!(ax1, bs, log10.(test[:,4]), color=:blue, linestyle = :dash)
    lines!(ax1, bs, log10.(test[:,5]), color=:red, linestyle = :dash)

    test = tests(solver_configs, linear_nodes)
    lines!(ax1, bs, log10.(test[:,4]), color=:blue, linestyle = :dot)
    lines!(ax1, bs, log10.(test[:,5]), color=:red, linestyle = :dot)
    display(f)
    return f
end
bs = collect(Float64,1:2:100)
plot_benchmarks(bs)


#x_int, y_int =  random_interior_points(triangle, 1000; grd = 1000)
#scatter!(x_int, y_int)

#=
f = Figure(resolution = (1000,1000));
grid = (200,200)
ax = Axis(f[1,1],xlabel=L"x", ylabel=L"y")
xlim, ylim = limits
x_plot = LinRange(xlim... , grid[1])
y_plot = LinRange(ylim... , grid[2])
x = repeat(x_plot , outer = length(y_plot))
y = repeat(y_plot , inner = length(x_plot))
phi = triangle.domain(x,y) #triangle.domain(x, y) 
Z = reshape(phi, grid)
heatmap!(ax, x_plot,y_plot,Z, colormap = :balance, colorrange=(-1,1))
ax.aspect=DataAspect()
display(f)
=#

k = 12.0
f = Figure(resolution = (1000,1000));
ax = Axis(f[1,1],xlabel=L"x", ylabel=L"y")

plot_basis_function!(ax, fb_basis[1], 1, k)
plot_billiard!(ax, triangle)
solver = ParticularSolutionsMethod(10.0,10.0,10.0) 
pts1 = evaluate_points(solver, triangle, gauss_legendre_nodes, k)
pts2 = evaluate_points(solver, triangle, linear_nodes, k)

pts1.x
pts1.x_int

scatter!(ax,pts1.x,pts1.y, ms=1.0, color= :black)
scatter!(ax,pts1.x_int,pts1.y_int, ms=10.0, marker=:x, color= :black)
display(f)
solver
B, B_int = construct_matrices(solver, fb_basis[1], pts1, k)

B_int

M = reduce(vcat, [B, B])
@btime V = reduce(directsum, [B, B, B])
@btime V1 = directsum(B,directsum(B,B))
V - V1
 #concatenate columns
Q, R =  qr(M)
_ , sv_r, Vt_r = svd(R)
mask = sv_r .> 1e-14
#println(mask)
V1t = Vt_r[ : , mask]
B * V1t


sv1 = svdvals((B * V1t), (B_int * V1t))
sv2 = svdvals(B, B_int)



GESVDVALS(B, B_int; eps=1e-14)
solution = svdvals(B, B_int)
size(B)
size(B_int)
solution.a
solution.b


M = reduce(vcat, [B, B_int]) #concatenate columns
solution = qr(M)
Q, R = qr(M) 
solution.Q, solution.R
solution = svd(R)
_ , sv_r, Vt_r = solution
println(sv_r)
mask = sv_r .> eps

println(mask)
V1t = Vt_r[ : , mask]


minimum(Diagonal(solution.D1))
d = @view solution.D1[:,:]
d = @view solution.D2[:,:]
d = diag(solution.D1) ./ diag(solution.D2)
d = diag(solution.D1).^2 .+ diag(solution.D2).^2

fig = plot_matrix(B_int)

R

#solver = ParticularSolutionsMethod(7.0,20.0,10.0)
#fb_basis

#=
function k_sweep(solver::ParticularSolutionsMethod, billiard::AbsBilliard, basis::AbsBasis, sampler, ks; eps=1e-14)
    L = real_length(billiard)
    N = round(Int, maximum(ks)*L*solver.dim_scaling_factor/(2*pi))
    println(N)
    #basis = rescale_basis(basis, N)
    res = similar(ks)
    for (i,k) in enumerate(ks)
        pts = evaluate_points(solver, billiard, sampler, k)
        B, B_int = construct_matrices(basis, pts, k)
        solution = GESVDVALS(B, B_int; eps=eps)

        res[i] = minimum(solution)
    end
    return res
end
=#
#=
function k_minimize(solver::ParticularSolutionsMethod, pts::PSMPoints, basis::AbsBasis, k_min, k_max)
    function f(k)
        B, B_int = construct_matrices(solver, basis, pts, k)
        solution = GESVDVALS(B, B_int; eps=solver.eps)
        return minimum(solution)
    end
    return optimize(f, k_min, k_max)
end
=#
#=
function k_sweep(solver::ParticularSolutionsMethod, pts::PSMPoints, basis::AbsBasis, ks)
    res = similar(ks)
    for (i,k) in enumerate(ks)
        B, B_int = construct_matrices(solver, basis, pts, k)
        solution = GESVDVALS(B, B_int; eps=solver.eps)
        res[i] = minimum(solution)
    end
    return res
end
=#