include("../abstracttypes.jl")
#include("curves.jl")

using PolygonOps
using StaticArrays

function approx_polygon(curves; grd=1000)
    x_bnd = Float64[]
    y_bnd = Float64[]
    for curve in curves
        L = curve.length
        N_bnd = round(Int, grd/L)
        t = LinRange(0.0,1.0, N_bnd)[1:end-1]
        x, y = curve.r(t)
        append!(x_bnd, x)
        append!(y_bnd, y)
    end
    x_bnd[end] = x_bnd[1]
    y_bnd[end] = y_bnd[1]
    polygon = [SVector(x,y) for (x,y) in zip(x_bnd ,y_bnd)]
    A =  PolygonOps.area(polygon)
    return polygon, A
end

#mask = [inpolygon(p, polygon; in=true, on=true, out=false) for p in new_pts]


function boundary_limits(curves; grd=1000)
    x_bnd = Float64[]
    y_bnd = Float64[]
    for curve in curves
        L = curve.length
        N_bnd = round(Int, grd/L)
        t = LinRange(0.0,1.0, N_bnd)[1:end-1]
        x, y = curve.r(t)
        append!(x_bnd, x)
        append!(y_bnd, y)
    end
    x_bnd[end] = x_bnd[1]
    y_bnd[end] = y_bnd[1]
    xlim = extrema(x_bnd)
    dx =  xlim[2] - xlim[1]
    ylim = extrema(y_bnd)
    dy =  ylim[2] - ylim[1]
    return xlim,ylim,dx,dy
end