
#include("../../abstracttypes.jl")
#include("../../utils/coordinatesystems.jl")
#include("../../utils/gridutils.jl")
#using BenchmarkTools
#using SpecialFunctions
using Bessels
using CoordinateTransformations, Rotations, StaticArrays

#Jv(nu, r) = SpecialFunctions.besselj(nu,r)
Jv(nu, r) = Bessels.besselj(nu, r)
ca_fb(nu,k,r,phi) = Jv(nu, k*r)*sin(nu*phi)
function Jvp(nu, r::T) where {T<:Real}
    let
    j_minus = Jv(nu-one(T),r)
    j_plus = Jv(nu+one(T),r)
    return 0.5*( j_minus - j_plus)
    end
end

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

struct CornerAdaptedFourierBessel{T,Sy} <: AbsBasis where  {T<:Real,Sy<:Union{AbsSymmetry,Nothing}}
    cs::PolarCS{T}
    dim::Int64 #using concrete type
    corner_angle::T
    nu::T #order constant, order=nu*i
    symmetries::Union{Vector{Sy},Nothing}
end

function CornerAdaptedFourierBessel(dim, corner_angle, origin, rot_angle)
    cs = PolarCS(origin, rot_angle)
    nu = pi/corner_angle
    return CornerAdaptedFourierBessel{Float64,Nothing}(cs, dim, corner_angle, nu, nothing)
end

function CornerAdaptedFourierBessel(dim, corner_angle, cs::CoordinateSystem)
    nu = pi/corner_angle
    return CornerAdaptedFourierBessel{Float64,Nothing}(cs, dim, corner_angle, nu, nothing)
end

toFloat32(basis::CornerAdaptedFourierBessel) = CornerAdaptedFourierBessel(basis.dim, Float32(basis.corner_angle), Float32.(basis.cs.origin), Float32(basis.cs.rot_angle))

function resize_basis(basis::CornerAdaptedFourierBessel, billiard::Bi, dim::Int, k) where {Bi<:AbsBilliard}
    if basis.dim == dim
        return basis
    else
        return CornerAdaptedFourierBessel(dim, basis.corner_angle, basis.cs)
    end
end

#evaluation functions
#=
function local_coords(basis::CornerAdaptedFourierBessel{T}, pt::SVector{2,T}) where {T<:Number}
     #new_pt::SVector{2,T} = cartesian_to_polar(f(pt))
    return cartesian_to_polar(basis.cs.local_map(pt))
end  
=#

@inline function basis_fun(basis::CornerAdaptedFourierBessel{T}, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map, nu=basis.nu, pts=pts
        pt_pol = (cartesian_to_polar(pm(pt)) for pt in pts)
        #norm::T = one(T)/sqrt(basis.dim)
        return collect(ca_fb(nu*i, k, pt[1], pt[2]) for pt in pt_pol)
    end
end

@inline function basis_fun(basis::CornerAdaptedFourierBessel{T}, indices::AbstractArray, k::T, pts::AbstractArray; parallel_matrix = true) where {T<:Real}
    let pm = basis.cs.local_map, nu=basis.nu, pts=pts
        pt_pol = (cartesian_to_polar(pm(pt)) for pt in pts)
        #norm::T = one(T)/sqrt(basis.dim)
        M =  length(pts)
        N = length(indices)
        B = zeros(T,M,N)
        if parallel_matrix
            Threads.@threads for i in eachindex(indices)
                B[:,i] .= (ca_fb(nu*i, k, pt[1], pt[2]) for pt in pt_pol)
            end
        else
            for i in eachindex(indices)
                B[:,i] .= (ca_fb(nu*i, k, pt[1], pt[2]) for pt in pt_pol)
            end
        end
        return B 
    end
end


@inline function dk_fun(basis::CornerAdaptedFourierBessel{T}, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    #translation of coordiante origin
    let pm = basis.cs.local_map, nu=basis.nu, pts=pts
        pt_pol = [cartesian_to_polar(pm(pt)) for pt in pts]
        #norm::T = one(T)/sqrt(basis.dim)
        r = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2) 
        dj = @. Jvp(nu*i, k*r)
        s = @. sin(nu*i*phi)
        dk = @. r*dj*s
        return dk
    end
end
    

@inline function dk_fun(basis::CornerAdaptedFourierBessel{T}, indices::AbstractArray, k::T, pts::AbstractArray; parallel_matrix = true) where {T<:Real}
    let pm = basis.cs.local_map, nu=basis.nu, pts=pts
        pt_pol = [cartesian_to_polar(pm(pt)) for pt in pts]
        #norm::T = one(T)/sqrt(basis.dim)
        r = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2)
        M =  length(pts)
        N = length(indices)
        dB_dk = zeros(T,M,N)
        if parallel_matrix
            Threads.@threads for i in eachindex(indices)
                dj = @. Jvp(nu*i, k*r)
                s = @. sin(nu*i*phi)
                dB_dk[:,i] .= @. r*dj*s
            end
        else
            for i in eachindex(indices)
                dj = @. Jvp(nu*i, k*r)
                s = @. sin(nu*i*phi)
                dB_dk[:,i] .= @. r*dj*s
            end
        end
        return dB_dk
    end
end

function gradient(basis::CornerAdaptedFourierBessel, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map, nu=basis.nu, pts=pts
        pt_xy = collect(pm(pt) for pt in pts)
        pt_pol = collect(cartesian_to_polar(pt) for pt in pt_xy) #local cartesian coords
        #norm::T = one(T)/sqrt(basis.dim)
        r = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2)
        X = getindex.(pt_xy,1)
        Y = getindex.(pt_xy,2)
        j = Jv.(nu*i, k*r)
        #println(size(j))
        dj = Jvp.(nu*i, k*r) 
        #println(size(dj))
        s = @. sin(nu*i*phi) 
        c = @. cos(nu*i*phi) 
        #println(size(s))
        dx = @. (dj*k*(X/r)*s-nu*i*j*c*Y/(r^2))
        dy = @. (dj*k*(Y/r)*s+nu*i*j*c*X/(r^2))
    return dx, dy
    end
end

function gradient(basis::CornerAdaptedFourierBessel, indices::AbstractArray, k::T, pts::AbstractArray; parallel_matrix=true) where {T<:Real}
    let pm = basis.cs.local_map, nu=basis.nu, pts=pts
        #local cartesian coords
        pt_xy = collect(pm(pt) for pt in pts)
        pt_pol = collect(cartesian_to_polar(pt) for pt in pt_xy)
        #norm::T = one(T)/sqrt(basis.dim)
        r = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2)
        X = getindex.(pt_xy,1)
        Y = getindex.(pt_xy,2)
        M = length(pts)
        N = length(indices)
        dB_dx = zeros(T,M,N)
        dB_dy = zeros(T,M,N)
        if parallel_matrix
            Threads.@threads for i in eachindex(indices)
                j = Jv.(nu*i, k*r)
                #println(size(j))
                dj = Jvp.(nu*i, k*r) 
                #println(size(dj))
                s = @. sin(nu*i*phi) 
                c = @. cos(nu*i*phi) 
                #println(size(s))
                dB_dx[:,i] .= @. (dj*k*(X/r)*s-nu*i*j*c*Y/(r^2))
                dB_dy[:,i] .= @. (dj*k*(Y/r)*s+nu*i*j*c*X/(r^2))
            end
        else
            for i in eachindex(indices)
                j = Jv.(nu*i, k*r)
                #println(size(j))
                dj = Jvp.(nu*i, k*r) 
                #println(size(dj))
                s = @. sin(nu*i*phi) 
                c = @. cos(nu*i*phi) 
                #println(size(s))
                dB_dx[:,i] .= @. (dj*k*(X/r)*s-nu*i*j*c*Y/(r^2))
                dB_dy[:,i] .= @. (dj*k*(Y/r)*s+nu*i*j*c*X/(r^2))
            end
        end #println(size(s))
    return dB_dx, dB_dy
    end
end


function basis_and_gradient(basis::CornerAdaptedFourierBessel, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map, nu=basis.nu, pts=pts
        pt_xy = collect(pm(pt) for pt in pts)
        pt_pol = collect(cartesian_to_polar(pt) for pt in pt_xy) #local cartesian coords
        #norm::T = one(T)/sqrt(basis.dim)
        r = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2)
        X = getindex.(pt_xy,1)
        Y = getindex.(pt_xy,2)
        j = Jv.(nu*i, k*r)
        #println(size(j))
        dj = Jvp.(nu*i, k*r) 
        #println(size(dj))
        s = @. sin(nu*i*phi) 
        c = @. cos(nu*i*phi) 

        #println(size(s))
        bf = @. j*s
        dx = @. (dj*k*(X/r)*s-nu*i*j*c*Y/(r^2))
        dy = @. (dj*k*(Y/r)*s+nu*i*j*c*X/(r^2))
    return bf, dx, dy
    end
end


function basis_and_gradient(basis::CornerAdaptedFourierBessel, indices::AbstractArray, k::T, pts::AbstractArray; parallel_matrix=true) where {T<:Real}
    let pm = basis.cs.local_map, nu=basis.nu, pts=pts
        #local cartesian coords
        pt_xy = collect(pm(pt) for pt in pts)
        pt_pol = collect(cartesian_to_polar(pt) for pt in pt_xy)
         
        #norm::T = one(T)/sqrt(basis.dim)
        r = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2)
        X = getindex.(pt_xy,1)
        Y = getindex.(pt_xy,2)
        M =  length(pts)
        N = length(indices)
        B = zeros(T,M,N)
        dB_dx = zeros(T,M,N)
        dB_dy = zeros(T,M,N)
        if parallel_matrix
            Threads.@threads for i in eachindex(indices)
                j = Jv.(nu*i, k*r)
                #println(size(j))
                dj = Jvp.(nu*i, k*r) 
                #println(size(dj))
                s = @. sin(nu*i*phi) 
                c = @. cos(nu*i*phi) 
                #println(size(s))
                B[:,i] .= @. j*s
                dB_dx[:,i] .= @. (dj*k*(X/r)*s-nu*i*j*c*Y/(r^2))
                dB_dy[:,i] .= @. (dj*k*(Y/r)*s+nu*i*j*c*X/(r^2))
            end
        else
            for i in eachindex(indices)
                j = Jv.(nu*i, k*r)
                #println(size(j))
                dj = Jvp.(nu*i, k*r) 
                #println(size(dj))
                s = @. sin(nu*i*phi) 
                c = @. cos(nu*i*phi) 
                #println(size(s))
                B[:,i] .= @. j*s
                dB_dx[:,i] .= @. (dj*k*(X/r)*s-nu*i*j*c*Y/(r^2))
                dB_dy[:,i] .= @. (dj*k*(Y/r)*s+nu*i*j*c*X/(r^2))
            end
        end
        #println(size(s))
    return B, dB_dx, dB_dy
    end
end

