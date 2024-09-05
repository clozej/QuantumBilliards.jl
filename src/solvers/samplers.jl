
using FastGaussQuadrature
using StaticArrays
using StatsBase

"""
This module provides a set of sampling strategies for numerical integration and other computational methods where specific distributions of points are required. It defines several types of samplers, each corresponding to a different method of generating sample points. The main samplers include:

- `LinearNodes`: Uniform sampling over an interval.
- `GaussLegendreNodes`: Nodes used in Gauss-Legendre quadrature for optimal polynomial integration.
- `FourierNodes`: Complex sampling method that can handle periodic domains and considers prime factors and custom lengths.

Each sampler type has an associated method for generating sample points, and additional utilities for working with random points inside defined boundaries.
"""

struct LinearNodes <: AbsSampler 
end 

"""
Generate sample points and corresponding intervals using uniform linear spacing.

# Arguments
- `sampler::LinearNodes`: The linear sampler type.
- `N::Int`: The number of intervals.

# Returns
- `t`: Vector of midpoints of the intervals between 0 and 1.
- `dt`: Vector of interval lengths between consecutive points. They are calculated via the diff method
"""
function sample_points(sampler::LinearNodes, N::Int)
    t = midpoints(range(0,1.0,length = (N+1)))
    dt = diff(range(0,1.0,length =(N+1)))
    return t, dt
end

struct GaussLegendreNodes <: AbsSampler 
end 

"""
Generate sample points and corresponding intervals using Gauss-Legendre quadrature.

# Arguments
- `sampler::GaussLegendreNodes`: The Gauss-Legendre sampler type.
- `N::Int`: The number of quadrature points.

# Returns
- `t`: Scaled Gauss-Legendre nodes mapped to the interval [0, 1].
- `dt`: Scaled weights corresponding to each node, used for integration.
"""
function sample_points(sampler::GaussLegendreNodes, N::Int)
    x, w = gausslegendre(N)
    t = 0.5 .* x  .+ 0.5
    dt = w .* 0.5 
    return t, dt
end


# TODO: Generate sample points using Chebyshev nodes
#=
function chebyshev_nodes(N::Int)
    x = [cos((2*i-1)/(2*N)*pi) for i in 1:N]
    t = 0.5 .* x  .+ 0.5
    dt = ones(N)  #wrong
    return t, dt
end
=#

struct FourierNodes <: AbsSampler where T<:Real
    primes::Union{Vector{Int64},Nothing}
    lengths::Union{Vector{Float64},Nothing} 
end 

FourierNodes() = FourierNodes(nothing,nothing)
FourierNodes(lengths::Vector{Float64}) = FourierNodes(nothing,lengths)

"""
Generate sample points and corresponding intervals using Fourier quadrature.

# Arguments
- `sampler::FourierNodes`: The Fourier sampler type.
- `N::Int`: The number of quadrature points.

# Returns
- `t`: Scaled Fourier nodes mapped to the interval [0, 1].
- `dt`: Scaled weights corresponding to each node, used for integration.
"""
function sample_points(sampler::FourierNodes, N::Int)
    if isnothing(sampler.primes) 
        M = N
    else
        M = nextprod(sampler.primes,N)
    end

    ts = Vector{Vector{Float64}}(undef,0)
    dts = Vector{Vector{Float64}}(undef,0)
    t::Vector{Float64} = Vector{Float64}(undef,0)
    dt::Vector{Float64} = Vector{Float64}(undef,0)
    if isnothing(sampler.lengths)
        t = collect(i/M for i in 0:(M-1))
        dt = diff(t)
        dt = push!(dt,dt[1])
        push!(ts,t)
        push!(dts,dt)
    else
        crv_lengths::Vector{Float64} = sampler.lengths
        L::Float64 = sum(crv_lengths)

        start::Float64 = 0.0
        dt_end::Float64 = 0.0
        ds::Float64 = 0.0
        for l in crv_lengths
            ds = L/(l*M) 
            println(start*ds)
            t = collect(range(start*ds,1.0,step=ds))
            #println(t)
            dt_end = 1.0 - t[end]
            start = (ds - dt_end)/ds
            push!(ts,t)
            dt = diff(t)
            push!(dt,dt_end)
            push!(dts,dt)
        end
    end
    
    return ts,dts
end

"""
Generate random points within a 2D boundary.

# Arguments
- `billiard::AbsBilliard`: The object representing the 2D boundary.
- `N::Int`: The number of random points to generate.
- `grd::Int`: Grid resolution used to determine the boundary limits. Defaults to 1000.

# Returns
- `pts`: A vector of randomly generated points within the boundary.
"""
function random_interior_points(billiard::AbsBilliard, N::Int; grd::Int = 1000)
    xlim,ylim = boundary_limits(billiard.fundamental_boundary; grd=grd)
    dx =  xlim[2] - xlim[1]
    dy =  ylim[2] - ylim[1]
    pts = []
 
    #println(length(pts))
    while length(pts)<N
        x = (dx .* rand() .+ xlim[1]) 
        y = (dy .* rand() .+ ylim[1])
        pt = SVector(x,y)
        if is_inside(billiard, [pt])[1] #rework this
            push!(pts,pt)
        end
    end
    return pts
end

# TODO: 
"""
Generate sample points based on Fourier nodes with prime factor considerations.

# Arguments
- `N::Int`: The base number of sample points.
- `primes`: Tuple of prime factors used to adjust the number of sample points.

# Returns
- `t`: Sample points distributed according to Fourier analysis.
- `dt`: Corresponding intervals between the points.
"""
#=
#needs some work
function fourier_nodes(N::Int; primes=(2,3,5)) #starts at 0 ends at 
    if primes == false
        M = N
    else
        M = nextprod(primes,N)
    end
    t = collect(i/M for i in 0:(M-1))
    dt = diff(t)
    dt = push!(dt,dt[1])
    return t, dt
end

function fourier_nodes(N::Int, crv_lengths; primes=(2,3,5)) #starts at 0 ends at 
    if primes == false
        M = N
    else
        M = nextprod(primes,N)
    end
    L = sum(crv_lengths)
    ts =Vector{Vector{typeof(L)}}(undef,0)
    dts =Vector{Vector{typeof(L)}}(undef,0)
    start = 0.0
    for l in crv_lengths
        ds = L/(l*M) 
        println(start*ds)
        t = collect(range(start*ds,1.0,step=ds))
        #println(t)
        dt_end = 1.0 - t[end]
        start = (ds - dt_end)/ds
        push!(ts,t)
        dt = diff(t)
        push!(dt,dt_end)
        push!(dts,dt)
    end
    return ts,dts
end
=#