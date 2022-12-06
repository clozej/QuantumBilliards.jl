
include("curves.jl")

liamcon_radial(t, lam)  = (1.0 .+ lam .* cos.(t .* pi))
limacon_x(t, lam ; x0 = zero(t), kwargs...) =  cos.(t .* pi) .* liamcon_radial(t, lam) .- x0
limacon_y(t, lam ; y0 = zero(t), kwargs...) = sin.(t .* pi) .* liamcon_radial(t, lam) .- y0
limacon_arclength(t, lam) = 2*(1+lam) .* Elliptic.E.(t .* pi/2, 4*lam/(1+lam)^2)

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