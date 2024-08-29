
#include("curves.jl")
"""
The r(t) representation of the limacon curve
"""
liamcon_radial(t, lam)  = (1.0 .+ lam .* cos.(t .* pi))
"""
The x projection of the r(t) for the limacon
"""
limacon_x(t, lam ; x0 = zero(t), kwargs...) =  cos.(t .* pi) .* liamcon_radial(t, lam) .- x0
"""
The y projection of the r(t) for the limacon
"""
limacon_y(t, lam ; y0 = zero(t), kwargs...) = sin.(t .* pi) .* liamcon_radial(t, lam) .- y0
"""
Returns the arclenght of the limacon curve via the elliptic integral
"""
limacon_arclength(t, lam) = 2*(1+lam) .* Elliptic.E.(t .* pi/2, 4*lam/(1+lam)^2)

"""
A structure representing a Limacon curve. The curve can be parameterized with a specific `lam` parameter, and optional translation and orientation.

# Fields
- `r::Function`: A function representing the curve in Cartesian coordinates.
- `n::Function`: A function representing the normal vector along the curve.
- `s::Function`: A function representing the arc length of the curve as a function of the parameter `t`.
- `ds::Function`: A function representing the derivative of the arc length (tangent length) along the curve.
- `orientation::Int`: The orientation of the curve, which affects the direction of the normal vector.
- `length::Real`: The total length of the curve.

# Constructor Arguments
- `lam`: Parameter controlling the shape of the Limacon curve.
- `x0`: Optional x-coordinate translation (default is `0.0`).
- `y0`: Optional y-coordinate translation (default is `0.0`).
- `orientation`: The orientation of the curve (default is `1`).

# Returns
- `LimaconCurve`: A `LimaconCurve` struct representing the specified Limacon curve.
"""
struct LimaconCurve  <: AbsRealCurve
    r::Function
    n::Function
    s::Function
    ds::Function
    orientation::Int
    length::Real
    #domain_function::Function
    function LimaconCurve(lam; x0 = 0.0, y0 = 0.0, orientation = 1)
        r = t -> curve(t, limacon_x, limacon_y, lam; x0=x0,y0=y0)
        n = t -> orientation .* normal_vec(t,limacon_x,limacon_y, lam)
        s = t -> limacon_arclength(t, lam)
        ds = t -> tangent_length(t, limacon_x, limacon_y, lam)
        L = s(1.0)
        return new(r,n,s,ds,orientation,L)
    end
end