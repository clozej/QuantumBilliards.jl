using GLMakie
using FFTW
using DSP
using Revise, BenchmarkTools
include("../src/solvers/samplers.jl")


lam = 0.5
k = 2*pi/lam
b = 100.0
N = round(Int,b/lam)
t, dt = fourier_nodes(N)
s = t #samples
sr = 1.0/dt[1]

u = 1.0*sin.(k.*s) .+ sin.(10.0*k.*s)
@btime rfft(u)
fu = rfft(u)
freq = rfftfreq(length(s),sr)
freq

length(s)

f = Figure(resolution = (1000,1000));
ax = Axis(f[1,1],xlabel=L"s", ylabel=L"u")
lines!(ax, s, u)
scatter!(ax, s, u)
ax2 = Axis(f[2,1],xlabel=L"\k", ylabel=L"ft")
lines!(ax2,freq, abs.(fu))
vlines!(ax2, [1/lam, 10.0/lam]; color=:black, linewidth=0.5)
display(f)


pg = periodogram(u, fs=sr)
pg.freq
f = Figure(resolution = (1000,1000));
ax = Axis(f[1,1],xlabel=L"s", ylabel=L"u")
lines!(ax, s, u)
scatter!(ax, s, u)
ax2 = Axis(f[2,1],xlabel=L"\lambda", ylabel=L"ft")
lines!(ax2,pg.freq, pg.power)
vlines!(ax2, [lam]; color=:black, linewidth=0.5)
display(f)


f = Figure(resolution = (1000,500));

display(f)

f = Figure(resolution = (1000,500));
ax = Axis(f[1,1],xlabel=L"k", ylabel=L"ft")
lines!(ax,ks, real(fu))
lines!(ax,ks, imag(fu))

display(f)

ks