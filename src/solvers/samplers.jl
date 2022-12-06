include("../abstracttypes.jl")
include("../utils/domainutils.jl")
using StatsBase
using FastGaussQuadrature

function linear_nodes(N)
    t = midpoints(range(0,1.0,length = (N+1)))
    dt = diff(range(0,1.0,length =(N+1)))
    return t, dt
end

function gauss_legendre_nodes(N)
    x, w = gausslegendre(N)
    t = 0.5 .* x  .+ 0.5
    dt = w .* 0.5 
    return t, dt
end

function chebyshev_nodes(N)
    x = [cos((2*i-1)/(2*N)*pi) for i in 1:N]
    t = 0.5 .* x  .+ 0.5
    dt = ones(N)  
    return t, dt
end

using PolygonOps
using StaticArrays
function random_interior_points(billiard::AbsBilliard, N; grd = 1000)
    A =  billiard.area
    #println(A)
    #println(N_int)
    xlim,ylim,dx,dy = boundary_limits(billiard.boundary; grd=grd)
    N_int = round(Int, N*dx*dy/A)
    pts_x = (dx .* rand(N_int) .+ xlim[1]) 
    pts_y = (dy .* rand(N_int) .+ ylim[1]) 
    #println(length(pts))
        
    function adjust_number_of_pts(pts_x, pts_y)

        mask = billiard.domain(pts_x, pts_y) 
        int_pts_x = pts_x[mask]
        int_pts_y = pts_y[mask]
        N_i = length(pts_x)
        if N_i < N
            #new_pts = [ ([dx, dy] .* v) .+ [xlim[1],ylim[1]] for v in rand(SVector{2, Float64}, )]
            new_pts_x = (dx .* rand(5*(N-N_i)) .+ xlim[1]) 
            new_pts_y = (dy .* rand(5*(N-N_i)) .+ ylim[1]) 
            
            mask = billiard.domain(new_pts_x, new_pts_y) 
            append!(int_pts_x, new_pts_x[mask])
            append!(int_pts_y, new_pts_y[mask])
        end
               
        if length(int_pts_x) >= N
            return int_pts_x[1:N], int_pts_y[1:N]
        else
            return adjust_number_of_pts(int_pts_x, int_pts_y)
        end
    end
    #println(mask)
    int_pts_x, int_pts_y = adjust_number_of_pts(pts_x, pts_y)
    
    return int_pts_x, int_pts_y
end

#=
function chebyshev_nodes(N)
    t = midpoints(range(0,1.0,N+1))
    dt = diff(range(0,1.0,N+1))
    return t, dt
end

function chebyshev_alt(N)
    t = midpoints(range(0,1.0,N+1))
    dt = diff(range(0,1.0,N+1))
    return t, dt
end

function half_chebyshev_nodes(N)
    t = midpoints(range(0,1.0,N+1))
    dt = diff(range(0,1.0,N+1))
    return t, dt
end
=#
#t, dt = gauss_legendre_nodes(10)
#t, dt = linear_trapez(10)

#dt
#t