#using  CoordinateTransformations, Rotations, StaticArrays
using GLMakie
using FFTW
using Revise, BenchmarkTools
#using Revise
#include("../src/billiards/coordinatesystems.jl")
#include("../src/billiards/geometry.jl")
include("../src/billiards/triangle.jl")
include("../src/billiards/billiard.jl")

include("../src/solvers/acceleratedmethods/acceleratedmethods.jl")
include("../src/solvers/acceleratedmethods/scalingmethod.jl")
include("../src/solvers/sweepmethods/sweepmethods.jl")
include("../src/solvers/sweepmethods/decompositionmethod.jl")
include("../src/states/eigenstates.jl")
include("../src/states/basisstates.jl")
include("../src/states/wavefunctions.jl")
include("../src/states/boundarycoherentstates.jl")
include("../src/states/gradients.jl")
include("../src/states/boundaryfunctions.jl")
include("../src/plotting/plottingmakie.jl")
include("../src/utils/benchmarkutils.jl")

s_grid = collect(range(0.0,1.0,512))
L = 1.0
q = 0.5
p = 1.0
k = 100.0

cs = collect(coherent(q,p,k,s,L,0) for s in s_grid)

f = Figure(resolution = (1000,500));
ax = Axis(f[1,1],xlabel=L"s", ylabel=L"c")
lines!(ax, s_grid, real.(cs))
lines!(ax, s_grid, imag.(cs))
display(f)

f = Figure(resolution = (1000,500));
ax = Axis(f[1,1],xlabel=L"s", ylabel=L"c")
lines!(ax, s_grid, abs2.(cs))
display(f)



gamma = 2/3*pi
chi  = 2.0

gamma = pi/2*(1+sqrt(5))/2
chi  = (1+sqrt(2))/2

billiard, basis = make_triangle_and_basis(gamma, chi)
basis32 = Float32(basis)
d = 5.0
b = 5.0
acc_solver = ScalingMethod(d,b)
sw_solver = DecompositionMethod(d,b)


k0 = 5.5
dk = 0.5
k, ten = solve_wavenumber(sw_solver,basis,billiard,k0,dk)
state = compute_eigenstate(sw_solver, basis, billiard, k)
state32 = Eigenstate(Float32(state.k), Float32.(state.vec))


k0 = 3005.8
dk = 0.05
k, ten = solve_wavenumber(acc_solver,basis,billiard,k0,dk)
state = compute_eigenstate(sw_solver, basis, billiard, k)
state32 = Eigenstate(Float32(state.k), Float32.(state.vec))
k

f = Figure(resolution = (1000,1000));
plot_probability!(f,state32, basis32, billiard; b=10.0, log =(true,-5), plot_normal=false, inside_only=true) 
display(f)

f = Figure(resolution = (1000,500));
ax = Axis(f[1,1],xlabel=L"s", ylabel=L"u")
plot_boundary_function!(ax,state, basis, billiard; b=50.0) 
display(f)

memory_size(state.vec)

u, s = boundary_function(state, basis, billiard; b = 10.0, sampler=fourier_nodes)
#regularize!(u)
fu = rfft(u)
sr = 1.0/diff(s)[1]
freq = rfftfreq(length(s),sr)

f = Figure(resolution = (1000,500));
ax = Axis(f[1,1],xlabel=L"s", ylabel=L"u")
lines!(ax, s, u)
vlines!(ax, curve_edge_lengths(billiard); color=:black, linewidth=0.5)
ax2 = Axis(f[2,1],xlabel=L"\frac{k}{2\pi}", ylabel=L"ft")
lines!(ax2,freq, log10.(abs.(fu)))
vlines!(ax2, [k/(2*pi)]; color=:black, linewidth=0.5)
display(f)

ft_bas = rfft(state.vec)
freq_bas = rfftfreq(length(state.vec))

f = Figure(resolution = (1000,500));
ax = Axis(f[1,1],xlabel=L"n", ylabel=L"c_n")
lines!(ax, state.vec)
ax2 = Axis(f[2,1],xlabel=L"k", ylabel=L"c_k")
lines!(ax2,freq_bas, abs.(ft_bas))
display(f)

findall(iszero,state.vec)


function antisym_vec(x)
    v = reverse(-x[2:end])
    return append!(v,x)
end

function husimi(k,u,s; c = 10.0, w = 7.0)
    #compute coherrent state weights
    N = length(s)
    sig = one(k)/sqrt(k) #width of the gaussian
    x = s[s.<=w*sig]
    idx = length(x) #do not change order here
    x = antisym_vec(x)
    a = one(k)/(2*pi*sqrt(pi*k)) #normalization factor in this version Hsimi is not noramlized to 1
    ds = (x[end]-x[1])/length(x)
    uc = CircularVector(u) #allows circular indexing
    gauss = @. exp(-k/2*x^2)*ds
    #construct evaluation points in p coordinate
    ps = collect(range(0.0,1.0,step = sig/c))
    #construct evaluation points in q coordinate
    q_stride = length(s[s.<=sig/c])
    q_idx = collect(1:q_stride:N)
    push!(q_idx,N) #add last point
    #println(length(q_idx))
    #println(q_idx)
    #println(N)
    qs = s[q_idx]
    #println(length(qs))
    H = zeros(typeof(k),length(qs),length(ps))
    for i in eachindex(ps)   
        cr = @. cos(ps[i]*k*x)*gauss #real part of coherent state
        ci = @. sin(ps[i]*k*x)*gauss #imag part of coherent state
        for j in eachindex(q_idx)
            u_w = uc[q_idx[j]-idx+1:q_idx[j]+idx-1] #window with relevant values of u
            hr = sum(cr.*u_w)
            hi = sum(ci.*u_w)
            H[j,i] = a*(hr*hr + hi*hi)
        end
    end
    return H, qs, ps    
end

range(1,100,step=10)

H, s, ps = husimi(k,u,s; w = 5.0)    
@code_warntype husimi(k,u,s; w = 7.0)
@btime husimi(k,u,s; w = 7.0)
sum(abs2.(u))
f = Figure(resolution = (1000,500));
hmap, ax = plot_heatmap!(f,s,ps,H)
vlines!(ax, curve_edge_lengths(billiard); color=:black, linewidth=0.5)
display(f)
H