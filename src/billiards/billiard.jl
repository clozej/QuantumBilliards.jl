
struct Billiard  <: AbsBilliard
    boundary :: Vector{<:AbsCurve}
    area :: Float64
    length:: Float64
    domain:: Function
end

function real_length(billiard::AbsBilliard)
    L = 0
    for curve in billiard.boundary
        if typeof(curve) <: AbsRealCurve
            L += curve.length
        end
    end
    return L 
end

function virtual_length(billiard::AbsBilliard)
    L = 0
    for curve in billiard.boundary
        if typeof(curve) <: AbsVirtualCurve
            L += curve.length
        end
    end
    return L 
end

function is_inside(x,y,billiard::AbsBilliard)
    return billiard.domain(x,y) .> 0.0 
end
