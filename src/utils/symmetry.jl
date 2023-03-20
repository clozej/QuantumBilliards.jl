
using CoordinateTransformations
using StaticArrays
#= #unfinished
function construct_reflection(line::SVector{2,T}) where T<:Real

    return LinearMap()
end

struct NoSym <: AbsSymmetry
end

=#

reflect_x = LinearMap(SMatrix{2,2}([-1.0 0.0;0.0 1.0]))
reflect_y = LinearMap(SMatrix{2,2}([1.0 0.0;0.0 -1.0]))

struct Reflection<: AbsSymmetry
    sym_map::LinearMap{SMatrix{2, 2, Float64, 4}}
    parity::Int64
end

function XReflection(parity)
    return Reflection(reflect_x, parity)
end

function YReflection(parity)
    return Reflection(reflect_y, parity)
end

function XYReflection(parity_x, parity_y)
    return Reflection(reflect_x ∘ reflect_y, parity_x*parity_y)
end

#=

using CairoMakie
x_grid = collect(range(0.0,2.0,400))
y_grid = collect(range(0.0,2.0,400))

pts = [SVector(x,y) for y in y_grid for x in x_grid]

pts1 = reflect_x.(pts)

f = Figure(resolution = (1500,900))
ax = Axis(f[1,1])
#ax, hmap = plot_heatmap_balaced!(f[2,3],x_grid,y_grid,dk_alt)
scatter!(ax, pts1)
display(f)
=#