
using CoordinateTransformations
using StaticArrays

"""
Reflection across the y-axis.

# Description
`reflect_x` is a `LinearMap` object that represents a reflection across the y-axis in a 2D coordinate system. This reflection inverts the x-coordinate of any point, effectively mirroring it across the y-axis.
- Applying this transformation to a point `(x, y)` results in the point `(-x, y)`.
"""
reflect_x = LinearMap(SMatrix{2,2}([-1.0 0.0;0.0 1.0]))

"""
Reflection across the x-axis.

# Description
`reflect_y` is a `LinearMap` object that represents a reflection across the x-axis in a 2D coordinate system. This reflection inverts the y-coordinate of any point, effectively mirroring it across the x-axis.
- Applying this transformation to a point `(x, y)` results in the point `(x, -y)`.
"""
reflect_y = LinearMap(SMatrix{2,2}([1.0 0.0;0.0 -1.0]))

"""
A struct representing a reflection symmetry in a 2D coordinate system.

# Description
The `Reflection` struct defines a reflection symmetry transformation in a 2D space. It includes the reflection mapping (`sym_map`), the parity of the transformation, and the axis along which the reflection occurs.

# Fields
- `sym_map::LinearMap{SMatrix{2, 2, Float64, 4}}`: The linear map representing the reflection transformation.
- `parity::Int64`: The parity of the transformation, which affects the sign of the wavefunction during reflection.
- `axis::Symbol`: The axis along which the reflection occurs. Possible values are `:y_axis`, `:x_axis`, or `:origin`.

# Arguments
- `sym_map`: A linear map representing the reflection.
- `parity`: An integer value representing the parity of the reflection.
- `axis`: A symbol indicating the axis of reflection.
"""
struct Reflection <: AbsSymmetry
    sym_map::LinearMap{SMatrix{2, 2, Float64, 4}}
    parity::Int64
    axis::Symbol
end

"""
Create a reflection symmetry across the y-axis.

# Description
This function returns a `Reflection` object that represents a reflection symmetry across the y-axis. The reflection is defined by the `reflect_x` linear map, which inverts the x-coordinates of the points.

# Arguments
- `parity`: The parity value that affects the sign of the wavefunction during the reflection.

# Returns
- A `Reflection` object representing the y-axis reflection.
"""
function XReflection(parity)
    return Reflection(reflect_x, parity, :y_axis)
end

"""
Create a reflection symmetry across the x-axis.

# Description
This function returns a `Reflection` object that represents a reflection symmetry across the x-axis. The reflection is defined by the `reflect_y` linear map, which inverts the y-coordinates of the points.

# Arguments
- `parity`: The parity value that affects the sign of the wavefunction during the reflection.

# Returns
- A `Reflection` object representing the x-axis reflection.
"""
function YReflection(parity)
    return Reflection(reflect_y, parity, :x_axis)
end

"""
Create a reflection symmetry across both the x-axis and y-axis (origin reflection).

# Description
This function returns a `Reflection` object that represents a reflection symmetry across both the x-axis and y-axis, effectively reflecting across the origin. The reflection is defined by the composition of `reflect_x` and `reflect_y`, and the parity is determined by the product of `parity_x` and `parity_y`.

# Arguments
- `parity_x`: The parity value for the x-axis reflection.
- `parity_y`: The parity value for the y-axis reflection.

# Returns
- A `Reflection` object representing the origin reflection.
"""
function XYReflection(parity_x, parity_y)
    return Reflection(reflect_x ∘ reflect_y, parity_x*parity_y, :origin)
end

"""
Reflect a wavefunction across specified symmetries.

# Description
This function reflects a wavefunction `Psi` across specified symmetries (e.g., x-axis, y-axis, origin) and adjusts the corresponding grids (`x_grid`, `y_grid`) accordingly. The function iterates over the provided symmetries and applies the appropriate reflection to both the wavefunction and the grids.

# Logic
- For each symmetry in `symmetries`:
  - If the symmetry is across the y-axis, the function reverses the `x_grid` and reflects the wavefunction along the x-axis.
  - If the symmetry is across the x-axis, the function reverses the `y_grid` and reflects the wavefunction along the y-axis.
- The wavefunction and grids are updated to include the reflected components.

# Arguments
- `Psi`: The wavefunction to be reflected, represented as a 2D array.
- `x_grid`: The grid of x-coordinates corresponding to the wavefunction.
- `y_grid`: The grid of y-coordinates corresponding to the wavefunction.
- `symmetries`: A list of `Reflection` objects representing the symmetries across which to reflect the wavefunction.

# Returns
- A tuple `(Psi, x_grid, y_grid)` where `Psi` is the reflected wavefunction, and `x_grid` and `y_grid` are the updated grids after reflection.
"""
function reflect_wavefunction(Psi,x_grid,y_grid,symmetries)
    for sym in symmetries
        if sym.axis == :y_axis
            x = -reverse(x_grid)
            Psi_ref = reverse(sym.parity.*Psi; dims=1)

            Psi = vcat(Psi_ref,Psi)
            x_grid = append!(x,x_grid)
        end
        if sym.axis == :x_axis
            y = -reverse(y_grid)
            Psi_ref = reverse(sym.parity.*Psi; dims=2)

            Psi = hcat(Psi_ref,Psi)
            y_grid = append!(y,y_grid)
        end
    end
    return Psi, x_grid, y_grid
end


#=
function reflect_wavefunction(Psi,x_grid,y_grid,symmetries)
    for sym in symmetries
        if sym.axis == :y_axis
            if x_grid[1] == zero(eltype(x_grid))
                x = -reverse(x_grid[2:end])
                Psi_ref = reverse(sym.parity.*Psi[2:end,:]; dims=1)
            else
                x = -reverse(x_grid)
                Psi_ref = reverse(sym.parity.*Psi; dims=1)
            end
            Psi = vcat(Psi_ref,Psi)
            x_grid = append!(x,x_grid)
        end
        if sym.axis == :x_axis
            if y_grid[1] == zero(eltype(y_grid))
                y = -reverse(y_grid[2:end])
                Psi_ref = reverse(sym.parity.*Psi[:,2:end]; dims=2)
            else
                y = -reverse(y_grid)
                Psi_ref = reverse(sym.parity.*Psi; dims=2)
            end
            Psi = hcat(Psi_ref,Psi)
            y_grid = append!(y,y_grid)
        end
    end
    return Psi, x_grid, y_grid
end
=#