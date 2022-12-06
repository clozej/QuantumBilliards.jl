
include("../abstracttypes.jl")
#using BenchmarkTools
#using SpecialFunctions
using Bessels

#Jv_alt(nu, r) = SpecialFunctions.besselj.(nu,r)
Jv(nu, r) = Bessels.besselj(nu, r)
ca_fb(nu,k,r,phi) = Jv(nu, k*r)*sin(nu*phi)
Jvp(nu, r) = 0.5*(Jv(nu-1,r) - Jv(nu+1,r))
ca_fb_dk(nu,k,r,phi) = r * Jvp(nu, k*r) * sin(nu * phi)

#=
function _bessel_diff_formula(v, z, n, L, phase)
    # from AMS55.
    # L(v, z) = J(v, z), Y(v, z), H1(v, z), H2(v, z), phase = -1
    # L(v, z) = I(v, z) or exp(v*pi*i)K(v, z), phase = 1
    # For K, you can pull out the exp((v-k)*pi*i) into the caller
    p = 1.0
    s = L(v-n, z)
    for i in 1:n
        p = phase * (p * (n-i+1)) / i   # = choose(k, i)
        s .+= p*L(v-n + i*2, z)
    end
    return s ./ (2.0^n)
end
=#


struct CornerAdaptedFourierBessel <: AbsBasis
    dim::Int
    type::Symbol
    nu::Float64 #order constant, order=nu*i
    phi0::Float64 #roatation angle
    x0::Float64
    y0::Float64
    branch_angle::Float64
    #sym_x:: Symbol
    #sym_y:: Symbol
    function CornerAdaptedFourierBessel(dim, corner_angle, phi0, x0, y0;branch_angle = 0.0 )
        nu = pi/corner_angle
        return new(dim, :FourrierBessel, nu, phi0, x0, y0, branch_angle)
    end
end

function rescale_basis(basis::CornerAdaptedFourierBessel, dim::Int)
    if basis.dim == dim
        return basis
    else
        return CornerAdaptedFourierBessel(dim, pi/basis.nu, basis.phi0, basis.x0, basis.y0;branch_angle = basis.branch_angle)
    end
end

@inline function basis_fun(basis::CornerAdaptedFourierBessel, i::Int, k, x::T, y::T) where T<:Number
    #translation of coordiante origin
    X = x - basis.x0 
    Y = y - basis.y0
    rot_angle =  basis.phi0 #.+
    #X = X .* cos(rot_angle) + Y .* sin(rot_angle)
    #Y = -X .* sin(rot_angle) + Y .* cos(rot_angle)
    phi = rem2pi(atan(Y, X) - rot_angle, RoundNearest)
    r = hypot(X, Y)
    norm = 1.0/sqrt(basis.dim)  
    return ca_fb(basis.nu*i, k, r, phi)/norm
end

@inline function basis_fun(basis::CornerAdaptedFourierBessel, i::Int, k, x::Vector{T}, y::Vector{T}) where T<:Number
    #translation of coordiante origin
    X = x .- basis.x0 
    Y = y .- basis.y0
    rot_angle =  basis.phi0 #.+
    #X = X .* cos(rot_angle) + Y .* sin(rot_angle)
    #Y = -X .* sin(rot_angle) + Y .* cos(rot_angle)
    phi = rem2pi.(atan.(Y, X) .- rot_angle, RoundNearest)
    r = hypot.(X, Y)
    norm = 1.0/sqrt.(basis.dim)  
    return ca_fb.(basis.nu*i, k, r, phi)./norm
end



#=
function basis_fun!(out_vec, basis::CornerAdaptedFourierBessel, i::Int, k, x::Vector{T}, y::Vector{T}) where T<:Number
    #translation of coordiante origin
    X = x .- basis.x0 
    Y = y .- basis.y0
    rot_angle =  basis.phi0 #.+
    #X = X .* cos(rot_angle) + Y .* sin(rot_angle)
    #Y = -X .* sin(rot_angle) + Y .* cos(rot_angle)
    phi = rem2pi.(atan.(Y, X) .- rot_angle, RoundNearest)
    r = hypot.(X, Y)
    norm = 1.0/sqrt(basis.dim)  
    out_vec .= ca_fb(basis.nu*i, k, r, phi) ./ norm
end
=#
@inline function dk_fun(basis::CornerAdaptedFourierBessel, i::Int, k, x::T, y::T) where T<:Number
    #translation of coordiante origin
    X = x - basis.x0 
    Y = y - basis.y0
    rot_angle =  basis.phi0 #.+
    #X = X .* cos(rot_angle) + Y .* sin(rot_angle)
    #Y = -X .* sin(rot_angle) + Y .* cos(rot_angle)
    phi = rem2pi(atan(Y, X) - rot_angle, RoundNearest)
    r = hypot(X, Y)
    norm = 1.0/sqrt(basis.dim) 
    return ca_fb_dk(basis.nu*i, k, r, phi)/norm
end

@inline function dk_fun(basis::CornerAdaptedFourierBessel, i::Int, k, x::Vector{T}, y::Vector{T}) where T<:Number
    #translation of coordiante origin
    X = x .- basis.x0 
    Y = y .- basis.y0
    rot_angle =  basis.phi0 #.+
    #X = X .* cos(rot_angle) + Y .* sin(rot_angle)
    #Y = -X .* sin(rot_angle) + Y .* cos(rot_angle)
    phi = rem2pi.(atan.(Y, X) .- rot_angle, RoundNearest)
    r = hypot.(X, Y)
    norm = 1.0/sqrt.(basis.dim)  
    return ca_fb_dk.(basis.nu*i, k, r, phi)./norm
end

#=

function dk_fun!(out_vec, basis::CornerAdaptedFourierBessel, i::Int, k, x::Vector{T}, y::Vector{T}) where T<:Number
    #translation of coordiante origin
    X = x .- basis.x0 
    Y = y .- basis.y0
    rot_angle =  basis.phi0 #.+
    #X = X .* cos(rot_angle) + Y .* sin(rot_angle)
    #Y = -X .* sin(rot_angle) + Y .* cos(rot_angle)
    phi = rem2pi.(atan.(Y, X) .- rot_angle, RoundNearest)
    r = hypot.(X, Y) 
    norm = 1.0/sqrt(basis.dim)
    out_vec .=  ca_fb_dk(basis.nu*i, k, r, phi) ./ norm
end

=#

function grad_fun(basis::CornerAdaptedFourierBessel, i::Int, k, x::T, y::T) where T<:Number
    #translation of coordiante origin
    X = x - basis.x0 
    Y = y - basis.y0
    rot_angle =  basis.phi0 #.+
    #X = X .* cos(rot_angle) + Y .* sin(rot_angle)
    #Y = -X .* sin(rot_angle) + Y .* cos(rot_angle)
    phi = rem2pi(atan(Y, X) - rot_angle, RoundNearest)
    r = hypot(X, Y)
    norm = 1.0/sqrt(basis.dim)
    j = Jv(basis.nu*i, r)
    #println(size(j))
    dj = Jvp(basis.nu*i, r)
    #println(size(dj))
    s = sin(basis.nu*i.*phi)
    c = cos(basis.nu*i.*phi)
    #println(size(s))
    dx = (dj*k*(X/r)*s-basis.nu*i*j*c*Y/(X^2+Y^2))/norm
    dy = (dj*k*(Y/r)*s+basis.nu*i*j*c*X/(X^2+Y^2))/norm
    return dx, dy
end

function grad_fun(basis::CornerAdaptedFourierBessel, i::Int, k, x::Vector{T}, y::Vector{T}) where T<:Number
    #translation of coordiante origin
    X = x .- basis.x0 
    Y = y .- basis.y0
    rot_angle =  basis.phi0 #.+
    #X = X .* cos(rot_angle) + Y .* sin(rot_angle)
    #Y = -X .* sin(rot_angle) + Y .* cos(rot_angle)
    phi = rem2pi.(atan.(Y, X) .- rot_angle, RoundNearest)
    r = hypot.(X, Y)
    norm = 1.0/sqrt.(basis.dim)  
    j = Jv.(basis.nu*i, r)
    #println(size(j))
    dj = Jvp.(basis.nu*i, r)
    #println(size(dj))
    s = sin.(basis.nu*i.*phi)
    c = cos.(basis.nu*i.*phi)
    #println(size(s))
    dx = @. (dj*k*(X/r)*s-basis.nu*i*j*c*Y/(X^2+Y^2))/norm
    dy = @. (dj*k*(Y/r)*s+basis.nu*i*j*c*X/(X^2+Y^2))/norm
    return dx, dy
end
