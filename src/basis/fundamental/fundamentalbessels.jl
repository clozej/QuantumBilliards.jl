

using Bessels
using CoordinateTransformations, Rotations, StaticArrays

#aliases so one is able to change bessel library easily
Y0(x) = Bessels.bessely0(x)
dY0(x) = -Bessels.bessely1(x)

struct FundamentalBessel{T} <: AbsBasis where  {T<:Real}
    #cs::PolarCS{T} #not fully implemented
    dim::Int64 #using concrete type
    dilation::T
    positions::Vector{SVector{2,T}}
end

function FundamentalBessel()
    return FundamentalBessel(0, 0.0, [SVector(0.0,0.0)])
end

function resize_basis(basis::Ba, billiard::Bi, dim::Int, k) where {Ba<:FundamentalBessel,Bi<:AbsBilliard}
    k_min = 25.0
    k_basis = max(k_min,k)
    pos = dilated_boundary_points(billiard, dim, k_basis; sampler=fourier_nodes, include_virtual=false)
    return FundamentalBessel(length(pos),2*pi/k_basis,pos)
end


@inline function basis_fun(basis::FundamentalBessel{T}, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let pts_loc = [pt - basis.positions[i] for pt in pts]
        arg = [k*hypot(pt[1], pt[2]) for pt in pts_loc]
        #norm::T = one(T)/sqrt(basis.dim)
        return Y0.(arg)
    end
end

#try using strided for this
@inline function basis_fun(basis::FundamentalBessel{T}, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let M =  length(pts)
        N = length(indices)
        B = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            pts_loc = [pt - basis.positions[i] for pt in pts]
            arg = [k*hypot(pt[1], pt[2]) for pt in pts_loc]
            B[:,i] .= Y0.(arg)
        end
        return B 
    end
end

function basis_and_gradient(basis::FundamentalBessel{T}, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let pts_loc = [pt - basis.positions[i] for pt in pts]
        r = [hypot(pt[1], pt[2]) for pt in pts_loc]
        x = getindex.(pts_loc, 1)
        y = getindex.(pts_loc, 2) 
        arg = k.*r
        #norm::T = one(T)/sqrt(basis.dim)
        bf = Y0.(arg)
        dx = @. k*x/r*dY0(arg)
        dy = @. k*y/r*dY0(arg)
    return bf, dx, dy
    end
end


function basis_and_gradient(basis::FundamentalBessel{T}, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let M =  length(pts)
        N = length(indices)
        B = zeros(T,M,N)
        dB_dx = zeros(T,M,N)
        dB_dy = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            pts_loc = [pt - basis.positions[i] for pt in pts]
            r = [hypot(pt[1], pt[2]) for pt in pts_loc]
            x = getindex.(pts_loc, 1)
            y = getindex.(pts_loc, 2) 
            arg = k.*r
            #println(size(s))
            B[:,i] .= Y0.(arg)
            dB_dx[:,i] .= @. k*x/r*dY0(arg)
            dB_dy[:,i] .= @. k*y/r*dY0(arg)
        end
        #println(size(s))
    return B, dB_dx, dB_dy
    end
end

@inline function dk_fun(basis::FundamentalBessel{T}, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let pts_loc = [pt - basis.positions[i] for pt in pts]
        r = [hypot(pt[1], pt[2]) for pt in pts_loc]
        #norm::T = one(T)/sqrt(basis.dim)
        return @. r*dY0(k*r)
    end
end
    

@inline function dk_fun(basis::FundamentalBessel{T}, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let M =  length(pts)
        N = length(indices)
        dB_dk = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            pts_loc = [pt - basis.positions[i] for pt in pts]
            r = [hypot(pt[1], pt[2]) for pt in pts_loc]
            arg = k.*r
            dB_dk[:,i] .= @. r*dY0(arg)
        end
        #println(size(s))
    return dB_dk
    end
end