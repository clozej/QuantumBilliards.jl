using MKL
using GLMakie
using Revise, BenchmarkTools
using LinearAlgebra, StaticArrays, Optim
using Statistics
include("../../src/billiards/triangle.jl")
include("../../src/billiards/billiard.jl")
include("../../src/basis/fourierbessel.jl")
include("../../src/plotting/plottingmakie.jl")
include("../../src/plotting/matrixplotting.jl")
include("../../src/solvers/scalingmethod.jl")
include("../../src/solvers/samplers.jl")
include("../../src/solvers/decompositions.jl")
include("../../src/utils/angleutils.jl")
include("../../src/utils/benchmarkutils.jl")
include("../../src/spectra/spectralutils.jl")

function benchmark_triangle(gamma,chi,k,dk; d = 10.0, b = 10.0, edge_i=1, sampler = gauss_legendre_nodes,
    btimes = 1, print_info=true, plot_info=true, fig_res=(1000,1000), log_mat=false)
    cor = Triangle(gamma,chi).corners
    x0,y0 = cor[mod1(edge_i+2,3)]
    #println("$x0, $y0")

    solver = ScalingMethod(d,b)
    re = [:Virtua, :Virtual, :Virtual]
    re[edge_i] = :Real 
    tr = Triangle(gamma,chi; curve_types = re, x0 = x0, y0 =y0)
    #dim = round(Int, tr.boundary[i].length*k*solver.dim_scaling_factor/(2*pi))
    basis = CornerAdaptedFourierBessel(1, adapt_basis(tr,edge_i+2,3)...) 
    #pts = evaluate_points(solver, tr, sampler, k)
    
    pad = (-0.2,0.2)
    limits = [extrema(row).+pad for row in eachrow(reduce(hcat, tr.corners))]
    
    ks, ten, decomp_time, mat_time=benchmark_solver(solver, basis, tr, sampler, k, dk; 
    btimes = btimes, print_info=print_info, plot_info=plot_info, fig_res=fig_res, log_mat=log_mat, limits=limits)

    return ks, ten, decomp_time, mat_time
end

function compute_triangle_spectrum(gamma,chi,k1,k2,dk; d = 3.0, b = 3.0, i=1, sampler = gauss_legendre_nodes,tol=1e-5)
    k = k2
    cor = Triangle(gamma,chi).corners
    x0,y0 = cor[mod1(i+2,3)]
    #println("$x0, $y0")
    solver = ScalingMethod(d,b)
    re = [:Virtua, :Virtual, :Virtual]
    re[i] = :Real 
    tr = Triangle(gamma,chi; curve_types = re, x0 = x0, y0 =y0)
    dim = round(Int, tr.boundary[i].length*k*solver.dim_scaling_factor/(2*pi))
    basis = CornerAdaptedFourierBessel(dim, adapt_basis(tr,i+2,3)...) 
    pts = evaluate_points(solver, tr, sampler, k)

    k_res, ten_res, check_res = compute_spectrum_test(solver, basis, pts, k1,k2,dk;tol=tol)
    return k_res, ten_res, check_res
end

#Strong-Mixing
gamma = pi/2*(1+sqrt(5))/2
chi  = (1+sqrt(2))/2

#Weak-Mixing
gamma = sqrt(2)/2 * pi
chi  = 2.0

#Veech
kind =7
n = 4 + kind
gamma = pi/2
chi = (n-2)/2

k = 4000.00001
dk = 0.005
ks, ten, decomp_time, mat_time = benchmark_triangle(gamma,chi,k,dk; d = 3.0, b = 10.0, edge_i=3, sampler = gauss_legendre_nodes)

k1 = k
k2 = k + 3*dk
ks, ten, control = compute_triangle_spectrum(gamma,chi,k1,k2,dk; d = 3.0, b = 4.0, i=3, sampler = gauss_legendre_nodes,tol=1e-5)
control


#=




function run_tests(gamma,chi,k,dk; d = 40.0, b=40.0, sampler=gauss_legendre_nodes)
    solver = ScalingMethod(d,b)
    #f = Figure(resolution = (1000,500));
    #axis = [Axis(f[1,i]) for i in 1:3]
    #ax = Axis(f[2:3,:],xlabel=L"k", ylabel=L"sv(k)")
    Fs = []
    Fks = []
    for i in 1:3
        println("")
        println("Computing single edge $i")
        re = [:Virtua, :Virtual, :Virtual]
        re[i] = :Real 
        tr = Triangle(gamma,chi; curve_types = re)
        dim = round(Int, tr.boundary[i].length*k*solver.dim_scaling_factor/(2*pi))
        al,phi0,x0,y0 = adapt_basis(tr,i+2,3)
        #println("angle=$alpha,  phi=$phi0")
        basis = CornerAdaptedFourierBessel(dim, al,phi0,x0,y0) 
        #push!(basis_vec,basis)
        #println(basis)
        #plot_basis_function!(axis[i],basis, 1, k; xlim=xlim.+ pad, ylim=ylim .+ pad)
        #plot_angle!(axis[i], al,phi0,x0,y0)
        #plot_billiard!(axis[i], tr)
        F, Fk, decomp_time, mat_time = benchmark_solver(solver, basis, tr, sampler, k, dk)
        push!(Fs,F)
        push!(Fks,Fk)
    end
    return Fs, Fks
end
=#





#Fs, Fks = run_tests(gamma,chi,k,dk;d = 40.0,b=40.0, sampler=gauss_legendre_nodes)

function compute_spectrum(gamma,chi,k1,k2,dk;d=4.0,b=5.0,id = 1)
    #Vebles method of merging
    #tol=0.2*dk
    x(k,k0,dk) = (k - k0)/dk
    rho(x) = (x>=0.0 && x<=1.0) ? 4.0*x*(1.0-x) : 0.0
    sig(ks) = sum(rho.(x.(ks,k,dk)))

    f = Figure(resolution = (1000,800));
    ax = Axis(f[1,1])
    ks1, ten1  = test(gamma,chi,k1,dk; d = d, b = b, i=id, sampler = gauss_legendre_nodes)
    ks2, ten2  = test(gamma,chi,k1+dk,dk; d = d, b = b, i=id, sampler = gauss_legendre_nodes)
    ks3, ten3  = test(gamma,chi,k1+0.5*dk,dk; d = d, b = b, i=id, sampler = gauss_legendre_nodes)

    scatter!(ax, ks1, log10.(ten1), color=(:blue, 0.5))
    errorbars!(ax, ks1, log10.(ten1), (0.5*ten1), (0.5*ten1); color=(:blue, 0.5), direction= :x)
    scatter!(ax, ks2, log10.(ten2), color=(:green, 0.5))
    errorbars!(ax, ks2, log10.(ten2), (0.5*ten2), (0.5*ten2); color=(:green, 0.5), direction= :x)
    scatter!(ax, ks3, log10.(ten3), color=(:red, 0.5))
    errorbars!(ax, ks3, log10.(ten3), (0.5*ten3), (0.5*ten3); color=(:red, 0.5), direction= :x)

    
    scatter!(ax, ks1,zeros(length(ks1)), color=(:blue, 0.5))
    errorbars!(ax, ks1, zeros(length(ks1)), (0.5*ten1), (0.5*ten1); color=(:blue, 0.5), direction= :x)
    scatter!(ax, ks2, zeros(length(ks2)), color=(:green, 0.5))
    errorbars!(ax, ks2, zeros(length(ks2)), (0.5*ten2), (0.5*ten2); color=(:green, 0.5), direction= :x)
    scatter!(ax, ks3, zeros(length(ks3)), color=(:red, 0.5))
    errorbars!(ax, ks3, zeros(length(ks3)), (0.5*ten3), (0.5*ten3); color=(:red, 0.5), direction= :x)


    #ax2 = Axis(f[1,2])
    #println((0.5*(sig(ks1)+sig(ks2))))    
    display(f)
end

compute_spectrum(gamma,chi,k,0.0,dk;d=4.0,b=5.0,id = 3)


#F1, F1k, sol1 = test(gamma,chi,k,dk; d = 4.0, b = 4.0, i=1, sampler = linear_nodes)
ks1, ten1  = test(gamma,chi,k,dk; d = 4.0, b = 5.0, i=1, sampler = gauss_legendre_nodes)
ks2, ten2  = test(gamma,chi,k,dk; d = 3.0, b = 5.0, i=2, sampler = gauss_legendre_nodes)
ks3, ten3  = test(gamma,chi,k,dk; d = 3.0, b = 5.0, i=3, sampler = gauss_legendre_nodes)
#minimum(ten)
x = collect(LinRange(1.0, 3.0, 10))
deleteat!(x,3:length(x))


using IntervalArithmetic
function is_equal(x,dx,y,dy)
    X = x ± dx
    Y = y ± dy 
    Z = X ∩ Y
    return  ~(Z == ∅)
end

function match_wavenumbers(ks_l,ts_l,ks_r,ts_r)
    #vectors ks_l and_ks_r must be sorted
    i = j = 1 #counting index
    control = Vector{Bool}()#control bits
    ks = Vector{eltype(ks_l)}()#final wavenumbers
    ts = Vector{eltype(ts_l)}()#final tensions
    while i <= length(ks_l) && j <= length(ks_r)
        x, dx = ks_l[i], ts_l[i]
        y, dy = ks_r[j], ts_r[j]
        if  is_equal(x,dx,y,dy) #check equality with errorbars
            i += 1 
            j += 1
            if dx < dy
                push!(ks, x)
                push!(ts, dx)
                push!(control, true)
            else
                push!(ks, y)
                push!(ts, dy)
                push!(control, true)
            end
        elseif x < y
            i += 1
            push!(ks, x)
            push!(ts, dx)
            push!(control, false)
        else 
            j += 1
            push!(ks, y)
            push!(ts, dy)
            push!(control, false)
        end
    end
    return ks, ts, control 
end

function overlap_and_merge!(k_left, ten_left, k_right, ten_right, control_left, kl, kr; tol=1e-3)
    #find overlaps in interval [k1,k2]
    idx_l = k_left .> (kl-tol) .&& k_left .< (kr+tol)
    idx_r = k_right .> (kl-tol) .&& k_right .< (kr+tol)
    
    ks_l,ts_l,ks_r,ts_r = k_left[idx_l], ten_left[idx_l], k_right[idx_r], ten_right[idx_r]
    #check if wavnumbers match in overlap interval
    ks, ts, control = match_wavenumbers(ks_l,ts_l,ks_r,ts_r)
    #println("left: $ks_l")
    #println("right: $ks_r")
    #println("overlaping: $ks")
    #i_l = idx_l[1]
    #i_r = idx_r[end]+1
    deleteat!(k_left, idx_l)
    append!(k_left, ks)
    deleteat!(ten_left, idx_l)
    append!(ten_left, ts)
    deleteat!(control_left, idx_l)
    append!(control_left, control)

    idx_last = findlast(idx_r) + 1
    append!(k_left, k_right[idx_last:end])
    append!(ten_left, ten_right[idx_last:end])
    append!(control_left, [false for i in idx_last:length(k_right)])

end

function compute_spectrum(solver::ScalingMethod, basis::AbsBasis, pts::SMPoints,k1,k2,dk;tol=1e-4)
    k0 = k1
    #initial computation
    k_res, ten_res = solve(solver, basis, pts, k0, dk+tol)
    control = [false for i in 1:length(k_res)]
    cycle = Makie.wong_colors()[1:6]
    f = Figure(resolution = (1000,1000));
    ax = Axis(f[1,1])
    scatter!(ax, k_res, log10.(ten_res),color=(cycle[1], 0.5))
    scatter!(ax, k_res, zeros(length(k_res)),color=(cycle[1], 0.5))
    #println("iteration 0")
    #println("merged: $k_res")
    
    #println("overlaping: $ks")
    i=1
    while k0 < k2
        #println("iteration $i")
        k0 += dk
        k_new, ten_new = solve(solver, basis, pts, k0, dk+tol)
        scatter!(ax, k_new, log10.(ten_new),color=(cycle[mod1(i+1,6)], 0.5))
        scatter!(ax, k_new, zeros(length(k_new)),color=(cycle[mod1(i+1,6)], 0.5))
        i+=1
        #println("new: $k_new")
        #println("overlap: $(k0-dk), $k0")
        overlap_and_merge!(k_res, ten_res, k_new, ten_new, control, k0-dk, k0; tol=tol)
        #println("merged: $k_res")
        #println("control: $control")
    end
    
    scatter!(ax, k_res, log10.(ten_res), color=(:black, 1.0), marker=:x,  ms = 100)
    scatter!(ax, k_res, zeros(length(k_res)), color=(:black, 1.0), marker=:x, ms = 100)
    display(f)
    return k_res, ten_res, control
end

function compute_triangle_spectrum(gamma,chi,k1,k2,dk; d = 3.0, b = 3.0, i=1, sampler = gauss_legendre_nodes,tol=1e-4)
    k = k2
    cor = Triangle(gamma,chi).corners
    x0,y0 = cor[mod1(i+2,3)]
    #println("$x0, $y0")
    solver = ScalingMethod(d,b)
    re = [:Virtua, :Virtual, :Virtual]
    re[i] = :Real 
    tr = Triangle(gamma,chi; curve_types = re, x0 = x0, y0 =y0)
    dim = round(Int, tr.boundary[i].length*k*solver.dim_scaling_factor/(2*pi))
    basis = CornerAdaptedFourierBessel(dim, adapt_basis(tr,i+2,3)...) 
    pts = evaluate_points(solver, tr, sampler, k)

    k_res, ten_res, check_res = compute_spectrum(solver, basis, pts, k1,k2,dk;tol=tol)
    return k_res, ten_res, check_res
end

k = 100.001
dk = 0.10014
ks, ts, check =compute_triangle_spectrum(gamma,chi,k,k+10*dk,dk; d = 5.0, b = 10.0, i=3, sampler = gauss_legendre_nodes)

rho(x) = (x>=0.0 && x<=1.0) ? 4.0*x*(1.0-x) : 0.0
f = Figure(resolution = (1000,1000));
ax = Axis(f[1,1])
xs = LinRange(-1.0,2.0,100)
lines!(ax,xs,rho.(xs))
display(f)

scatter!(ax, ks1, log10.(ten1), color=(:blue, 0.5))
scatter!(ax, ks2, log10.(ten2), color=(:green, 0.5))
scatter!(ax, ks3, log10.(ten3), color=(:red, 0.5))

x(k,k0,dk) = (k .- k0)./dk
cycle = Makie.wong_colors()[1:6]
x(ks,k,dk)
check

is_equal(0.1,0.01,0.1,0.02)

X = 1.0 ± 0.5

idx = x .> 1.3 .&& x .< 2.5  
x[idx]

