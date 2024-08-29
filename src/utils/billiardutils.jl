"""
Calculates the total arclength of the desymmetrized boundary of the given billiard. This is the actual ("physical") boundary

# Logic:
- For each curve in the desymmetrized composite boundary we check if it a real boundary element (instead of a virtual line segment that forms due to symmetries (like a Neumann/Dirichlet ficticious boundary)). It that is the case we add this length to the total lenght of the real/physical boundary

# Returns
The real boundary length of the fundamental boundary
"""
function real_length(billiard::Bi) where Bi<:AbsBilliard
    L = 0.0
    for curve in billiard.fundamental_boundary
        if typeof(curve) <: AbsRealCurve
            L += curve.length
        end
    end
    return L 
end

"""
Calculates the total arclength of the virtual desymmetrized boundary of the given billiard. This is the virtual ("non-physical") boundary that togeher with the real ("physical") boundary forms a closed boundary

# Logic:
- For each curve in the desymmetrized composite boundary we check if it a virtual boundary element (instead of a real line segment that is given by initial construction. It that is the case we add this length to the total lenght of the virtual boundary

# Returns
The total virtual boundary length of the fundamental boundary
"""
function virtual_length(billiard::Bi) where Bi<:AbsBilliard
    L = 0.0
    for curve in billiard.fundamental_boundary
        if typeof(curve) <: AbsVirtualCurve
            L += curve.length
        end
    end
    return L 
end

"""
Calculates all the arclengths of the desymmetrized boundary of the given billiard. This is a vector of real ("physical") boundares that together with the virtual boundary forms a closed boundary

# Logic:
- For each curve in the desymmetrized composite boundary we check if it a virtual boundary element (instead of a real line segment that is given by initial construction. It that is the case we add this length to the total lenght of the virtual boundary

# Returns
A vector of real ("physical") and virtual boundary lengths
"""
function curve_edge_lengths(billiard::Bi) where Bi<:AbsBilliard
    L = 0.0
    res = [L]
    for crv in billiard.full_boundary
        if (typeof(crv) <: AbsRealCurve)
            L += crv.length
            push!(res,L)
        end 
    end
    return res
end

"""
High level wrapper for checking if a given point is inside a boundary (either desymetrized/fundamental or the full one). Under the hood it calls a check if for each curve (either real or virtual/desymmetrized bases in the fundamental_boundary kwarg) the given point is inside. If true it must be inside a the given fundamental or full boundary 

# Returns
A boolean indicating whether the point is inside the boundary
"""
function is_inside(billiard::Bi, pt; fundamental_domain = true ) where Bi<:AbsBilliard
    if fundamental_domain 
        boundary = billiard.fundamental_boundary  
    else
        boundary = billiard.full_boundary
    end
    return all(is_inside(crv, pt) for crv in boundary) 
end

"""
High level wrapper for checking if a given vector of points is inside a fundamental boundary or of a full boundary. 
    
# Logic
- Under the hood it checks whether each point in the vector of points `pts` is inside the curves that form the full or desymmetrized/fundamental boundary. If true all the `pts` must be inside a the given boundary

# Returns
A boolean indicating whether all the points are inside the boundary
"""
function is_inside(billiard::Bi, pts::AbstractArray; fundamental_domain = true) where Bi<:AbsBilliard
    let 
        if fundamental_domain 
            curves = billiard.fundamental_boundary  
        else
            curves = billiard.full_boundary
        end
    
        inside = is_inside(curves[1], pts)
        for i in 2:length(curves)
            inside = inside .& is_inside(curves[i], pts)
        end
        return inside
    end 
end

"""
Calculates the extremal values (min/max) of the x and y coordinates of the curves forming the boundary. 

# Logic
- Each curve is homogenously discretrized with at least 512 points and the (x,y) coordianates are appended to the x and y vectors that contain index wise the x and y coordinates of the curves.
- We close the x and y vectors of points (if not already closed) since we must have a closed boundary
- Use the extrema function to calculate the minimum and maximum values for x and y

# Returns 
The (min,max) value of two vectors x and y as a 2 vector (x_min, x_max) , (y_min, y_max)
"""
function boundary_limits(curves; grd=1000) 
    x_bnd = Vector{Any}()
    y_bnd = Vector{Any}()
    for crv in curves #names of variables not very nice
        L = crv.length
        N_bnd = max(512,round(Int, grd/L))
        t = range(0.0,1.0, N_bnd)[1:end-1]
        pts = curve(crv,t)
        append!(x_bnd, getindex.(pts,1))
        append!(y_bnd, getindex.(pts,2))
    end
    x_bnd[end] = x_bnd[1]
    y_bnd[end] = y_bnd[1]
    xlim = extrema(x_bnd)
    #dx =  xlim[2] - xlim[1]
    ylim = extrema(y_bnd)
    #dy =  ylim[2] - ylim[1]
    return xlim, ylim #,dx,dy
end

