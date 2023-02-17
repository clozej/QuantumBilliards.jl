include("../abstracttypes.jl")
include("../billiards/billiard.jl")
include("../utils/gridutils.jl")
using StatsBase
using FastGaussQuadrature
using StaticArrays

function linear_nodes(N::Int)
    t = midpoints(range(0,1.0,length = (N+1)))
    dt = diff(range(0,1.0,length =(N+1)))
    return t, dt
end

function fourier_nodes(N::Int; primes=(2,3,5)) #starts at 0 ends at 
    M = nextprod(primes,N)
    t = collect(i/M for i in 0:(M-1))
    dt = diff(t)
    dt = push!(dt,dt[1])
    return t, dt
end

function gauss_legendre_nodes(N::Int)
    x, w = gausslegendre(N)
    t = 0.5 .* x  .+ 0.5
    dt = w .* 0.5 
    return t, dt
end

function chebyshev_nodes(N::Int)
    x = [cos((2*i-1)/(2*N)*pi) for i in 1:N]
    t = 0.5 .* x  .+ 0.5
    dt = ones(N)  #wrong
    return t, dt
end

function random_interior_points(billiard::AbsBilliard, N::Int; grd::Int = 1000)
    xlim,ylim = boundary_limits(billiard.boundary; grd=grd)
    dx =  xlim[2] - xlim[1]
    dy =  ylim[2] - ylim[1]
    pts = []
 
    #println(length(pts))
    while length(pts)<=N
        x = (dx .* rand() .+ xlim[1]) 
        y = (dy .* rand() .+ ylim[1])
        pt = SVector(x,y)
        if is_inside(billiard, pt)
            push!(pts,pt)
        end
    end
    return pts
end

