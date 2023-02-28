#using LazyGrids
#used for plotting etc.
#include("../abstracttypes.jl")
#include("billiardutils.jl")
using StaticArrays



struct GridPoint{T} <:AbsGrid where T <: Number
    xy::SVector{2,T} #point coords
    inside::Bool #true means point is inside billiard
    idx::CartesianIndex{2}
    #idx::Int #consecutive index
end

function make_grid_point(idx,x_grid,y_grid,billiard) 
    x, y = x_grid[idx[1]],y_grid[idx[2]] #for capturing variables
    pt = SVector(x,y)
    inside = is_inside(billiard,pt) 
    return GridPoint(pt, inside ,idx)      
end

function interior_grid(billiard::AbsBilliard,size::Tuple{Int,Int},xlim,ylim)
    type = eltype(xlim)
    x_grid = collect(type,range(xlim... , size[1]))
    y_grid = collect(type,range(ylim... , size[2]))

    grid_indices = CartesianIndices(size)
    generator = (make_grid_point(idx,x_grid,y_grid,billiard) for idx in grid_indices)
    return x_grid, y_grid, generator
end

function interior_grid(billiard::AbsBilliard,size::Tuple{T,T},xlim,ylim) where T <:AbstractFloat
    dx, dy = size
    type = eltype(xlim)
    x_grid = collect(type,range(xlim... ;step=dx))
    y_grid = collect(type,range(ylim... ;step=dy))
    grid_indices = CartesianIndices((length(x_grid),length(y_grid)))
    generator = (make_grid_point(idx,x_grid,y_grid,billiard) for idx in grid_indices)
    return x_grid, y_grid, generator
end

function interior_grid(billiard::AbsBilliard,size)
    xlim,ylim = boundary_limits(billiard.boundary; grd=1000)
    println(xlim,ylim)
    return interior_grid(billiard,size,xlim,ylim)
end

#lazy grid generators


function make_grid(size::Tuple{Int,Int},xlim,ylim)
    type = eltype(xlim)
    x_grid = collect(type,range(xlim... , size[1]))
    y_grid = collect(type,range(ylim... , size[2]))
    grid_indices = CartesianIndices((length(x_grid),length(y_grid)))
    generator = (SVector( x_grid[idx[1]],y_grid[idx[2]]) for idx in grid_indices)
    return x_grid, y_grid, generator
end

function make_grid(size::Tuple{T,T},xlim,ylim) where T <:AbstractFloat
    dx, dy = size
    type = eltype(xlim)
    x_grid = collect(type,range(xlim... ;step=dx))
    y_grid = collect(type,range(ylim... ;step=dy))
    grid_indices = CartesianIndices((length(x_grid),length(y_grid)))
    generator = (SVector( x_grid[idx[1]],y_grid[idx[2]]) for idx in grid_indices)
    return x_grid, y_grid, generator
end

#=
struct InteriorGrid{T,B} <: AbsGrid where {T <: Number, B <: AbsBilliard}
    billiard::Bi
    size::Tuple{Int,Int}
    x_grid::AbstractVector{T}
    y_grid::AbstractVector{T}
    #count::Int
end

function InteriorGrid(billiard,size,xlim,ylim)
    x_grid = range(xlim... , size[1])
    y_grid = range(ylim... , size[2])
    return InteriorGrid(billiard,size,x_grid,y_grid)
end

function InteriorGrid(billiard,size)
    x_grid = range(xlim... , size[1])
    y_grid = range(ylim... , size[2])
    return InteriorGrid(billiard,size,x_grid,y_grid)
end

=#

    #=
function Base.iterate(g::InteriorGrid)
    xy = SVector(g.x_grid[1],g.y_grid[1])
    inside = is_inside(g.billiard, xy)
    return GridPoint(xy,inside,1,1,1), 1
end

function Base.iterate(g::InteriorGrid, state=1 )
    if idx
    nx,ny = g.size
end
=#
#=
function lazy_grid(xlim, ylim, grd::Tuple)
    x = collect(range(xlim...,grd[1]))
    y = collect(range(ylim...,grd[2]))
    (xg, yg) = ndgrid(x, y)
    return x, y, xg, yg
end
=#
#x, y, xg, yg = lazy_grid((0,1), (0,1), (100,100))

