#using LazyGrids
#used for plotting etc.
#include("../abstracttypes.jl")
#include("billiardutils.jl")
using StaticArrays

"""
A struct representing a point on a grid within a 2D billiard domain.

# Description
`GridPoint` represents a point on a 2D grid, including its coordinates, whether it is inside a billiard domain, and its index on the grid.

# Fields
- `xy::SVector{2,T}`: The 2D coordinates of the grid point.
- `inside::Bool`: A boolean indicating whether the point lies inside the billiard domain (`true`) or outside (`false`).
- `idx::CartesianIndex{2}`: The index of the point on the grid in Cartesian coordinates.
"""
struct GridPoint{T} <:AbsGrid where T <: Number
    xy::SVector{2,T} #point coords
    inside::Bool #true means point is inside billiard
    idx::CartesianIndex{2}
    #idx::Int #consecutive index
end

"""
Create a `GridPoint` at a specified grid index within a billiard domain.

# Description
This function creates a `GridPoint` object by determining the coordinates at the given index on the grid, checking whether the point is inside the billiard domain, and then packaging this information into a `GridPoint` struct.

# Logic
- The x and y coordinates are extracted from `x_grid` and `y_grid` using the provided index.
- The point is then checked to see if it lies inside the billiard domain.
- A `GridPoint` object is returned with the point's coordinates, inside status, and index.

# Arguments
- `idx`: The Cartesian index on the grid.
- `x_grid`: The grid of x-coordinates.
- `y_grid`: The grid of y-coordinates.
- `billiard`: The billiard domain used to check if the point is inside.

# Returns
- A `GridPoint` object representing the point at the specified index.
"""
function make_grid_point(idx,x_grid,y_grid,billiard) 
    x, y = x_grid[idx[1]],y_grid[idx[2]] #for capturing variables
    pt = SVector(x,y)
    inside = is_inside(billiard,pt) 
    return GridPoint(pt, inside ,idx)      
end

"""
Generate a grid of points within a billiard domain with specified limits.

# Description
This function generates a 2D grid of points within the given billiard domain, constrained by specified x and y limits. It returns the x and y grids along with a generator for creating `GridPoint` objects at each grid index.

# Logic
- The x and y grids are created using `range` within the provided limits and grid size.
- A generator is created to produce `GridPoint` objects for each Cartesian index in the grid.

# Arguments
- `billiard`: The billiard domain used to check if points are inside.
- `size::Tuple{Int,Int}`: The size of the grid (number of points in x and y directions).
- `xlim`: The limits of the x-coordinates.
- `ylim`: The limits of the y-coordinates.

# Returns
- A tuple `(x_grid, y_grid, generator)` where:
  - `x_grid`: The grid of x-coordinates.
  - `y_grid`: The grid of y-coordinates.
  - `generator`: A generator that produces `GridPoint` objects for each point in the grid.
"""
function interior_grid(billiard::AbsBilliard,size::Tuple{Int,Int},xlim,ylim)
    type = eltype(xlim)
    x_grid = collect(type,range(xlim... , size[1]))
    y_grid = collect(type,range(ylim... , size[2]))

    grid_indices = CartesianIndices(size)
    generator = (make_grid_point(idx,x_grid,y_grid,billiard) for idx in grid_indices)
    return x_grid, y_grid, generator
end

"""
Generate a grid of points within a billiard domain using step sizes.

# Description
This function generates a 2D grid of points within the given billiard domain, using the specified step sizes for x and y coordinates. It returns the x and y grids along with a generator for creating `GridPoint` objects at each grid index.

# Logic
- The x and y grids are created using `range` with the provided step sizes.
- A generator is created to produce `GridPoint` objects for each Cartesian index in the grid.

# Arguments
- `billiard`: The billiard domain used to check if points are inside.
- `size::Tuple{T,T}`: The step sizes for the x and y directions.
- `xlim`: The limits of the x-coordinates.
- `ylim`: The limits of the y-coordinates.

# Returns
- A tuple `(x_grid, y_grid, generator)` where:
  - `x_grid`: The grid of x-coordinates.
  - `y_grid`: The grid of y-coordinates.
  - `generator`: A generator that produces `GridPoint` objects for each point in the grid.
"""
function interior_grid(billiard::AbsBilliard,size::Tuple{T,T},xlim,ylim) where T <:AbstractFloat
    dx, dy = size
    type = eltype(xlim)
    x_grid = collect(type,range(xlim... ;step=dx))
    y_grid = collect(type,range(ylim... ;step=dy))
    grid_indices = CartesianIndices((length(x_grid),length(y_grid)))
    generator = (make_grid_point(idx,x_grid,y_grid,billiard) for idx in grid_indices)
    return x_grid, y_grid, generator
end

"""
Generate a grid of points within a billiard domain with automatically determined limits using boundaryutils.jl.

# Description
This function generates a 2D grid of points within the given billiard domain, automatically determining the limits based on the billiard's boundary. It returns the x and y grids along with a generator for creating `GridPoint` objects at each grid index.

# Logic
- The x and y limits are automatically determined from the billiard's boundary.
- The `interior_grid` function is then called with these limits and the specified grid size. This function is a lower level helper.

# Arguments
- `billiard`: The billiard domain used to check if points are inside.
- `size`: The size of the grid.

# Returns
- A tuple `(x_grid, y_grid, generator)` where:
  - `x_grid`: The grid of x-coordinates.
  - `y_grid`: The grid of y-coordinates.
  - `generator`: A generator that produces `GridPoint` objects for each point in the grid.
"""
function interior_grid(billiard::AbsBilliard,size)
    xlim,ylim = boundary_limits(billiard.fundamental_boundary; grd=1000)
    #println(xlim,ylim)
    return interior_grid(billiard,size,xlim,ylim)
end

#lazy grid generators

"""
Generate a grid of 2D points within specified limits.

# Description
This function generates a 2D grid of points within the given x and y limits. It returns the x and y grids along with a generator for creating 2D points (`SVector`) at each grid index.

# Logic
- The x and y grids are created using `range` within the provided limits and wanted grid size.
- A generator is created to produce `SVector{2,T}` points for each Cartesian index in the grid.

# Arguments
- `size::Tuple{Int,Int}`: The size of the grid (number of points in x and y directions).
- `xlim`: The limits of the x-coordinates.
- `ylim`: The limits of the y-coordinates.

# Returns
- A tuple `(x_grid, y_grid, generator)` where:
  - `x_grid`: The grid of x-coordinates.
  - `y_grid`: The grid of y-coordinates.
  - `generator`: A generator that produces `SVector{2,T}` points for each point in the grid.
"""
function make_grid(size::Tuple{Int,Int},xlim,ylim)
    type = eltype(xlim)
    x_grid = collect(type,range(xlim... , size[1]))
    y_grid = collect(type,range(ylim... , size[2]))
    grid_indices = CartesianIndices((length(x_grid),length(y_grid)))
    generator = (SVector( x_grid[idx[1]],y_grid[idx[2]]) for idx in grid_indices)
    return x_grid, y_grid, generator
end

"""
Generate a grid of 2D points within specified limits using step sizes. This is the abstracted version of the `Int` version

# Description
This function generates a 2D grid of points within the given x and y limits, using the specified step sizes. It returns the x and y grids along with a generator for creating 2D points (`SVector`) at each grid index.

# Logic
- The x and y grids are created using `range` with the provided step sizes.
- A generator is created to produce `SVector{2,T}` points for each Cartesian index in the grid.

# Arguments
- `size::Tuple{T,T}`: The step sizes for the x and y directions.
- `xlim`: The limits of the x-coordinates.
- `ylim`: The limits of the y-coordinates.

# Returns
- A tuple `(x_grid, y_grid, generator)` where:
  - `x_grid`: The grid of x-coordinates.
  - `y_grid`: The grid of y-coordinates.
  - `generator`: A generator that produces `SVector{2,T}` points for each point in the grid.
"""
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

