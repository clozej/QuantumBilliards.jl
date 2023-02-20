#=
struct Billiard  <: AbsBilliard
    boundary :: Vector{<:AbsCurve}
    area :: Float64
    length:: Float64
    domain:: Function
end
=#

include("../abstracttypes.jl")

function real_length(billiard::Bi) where Bi<:AbsBilliard
    L = 0.0
    for curve in billiard.boundary
        if typeof(curve) <: AbsRealCurve
            L += curve.length
        end
    end
    return L 
end

function virtual_length(billiard::Bi) where Bi<:AbsBilliard
    L = 0.0
    for curve in billiard.boundary
        if typeof(curve) <: AbsVirtualCurve
            L += curve.length
        end
    end
    return L 
end

function curve_edge_lengths(billiard::Bi) where Bi<:AbsBilliard
    L = 0.0
    res = [L]
    for curve in billiard.boundary
        L += curve.length
        push!(res,L) 
    end
    return res
end


function is_inside(billiard::Bi, pt) where Bi<:AbsBilliard
    return all(is_inside(crv, pt) for crv in billiard.boundary) 
end


function is_inside(billiard::Bi, pts::AbstractArray) where Bi<:AbsBilliard
    let curves = billiard.boundary
        inside = is_inside(curves[1], pts)
        for i in 2:length(curves)
            inside = inside .& is_inside(curves[i], pts)
        end
        return inside
    end 
end

function is_inside(billiard::Bi, x_grid::AbstractArray, y_grid::AbstractArray) where Bi<:AbsBilliard
    let curves = billiard.boundary
        inside = is_inside(curves[1], x_grid, y_grid)
        for i in 2:length(curves)
            inside = inside .& is_inside(curves[i], x_grid, y_grid)
        end
        return inside
    end 
end
