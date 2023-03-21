

using Bessels
using CoordinateTransformations, Rotations, StaticArrays

#aliases so one is able to change bessel library easily
Y0(x) = Bessels.bessely0(x)
dY0(x) = -Bessels.bessely1(x)

struct FundamentalBessel{T,Sy} <: AbsBasis where  {T<:Real, Sy<:Union{AbsSymmetry,Nothing}}
    #cs::PolarCS{T} #not fully implemented
    dim::Int64 #using concrete type
    dilation::T
    positions::Vector{SVector{2,T}}
    symmetries::Union{Vector{Sy},Nothing}
end

function FundamentalBessel()
    return FundamentalBessel{Float64,Nothing}(0, 0.0, [SVector(0.0,0.0)], nothing)
end

function FundamentalBessel(symmetries::Vector{Sy}) where  {Sy<:AbsSymmetry}
    return FundamentalBessel(0, 0.0, [SVector(0.0,0.0)], symmetries)
end

function resize_basis(basis::Ba, billiard::Bi, dim::Int, k) where {Ba<:FundamentalBessel,Bi<:AbsBilliard}
    k_min = 25.0
    k_basis = max(k_min,k)
    if ~isnothing(basis.symmetries)
        pos = dilated_boundary_points(billiard, dim, k_basis; sampler=linear_nodes, include_virtual=false)
    else
        pos = dilated_boundary_points(billiard, dim, k_basis; sampler=fourier_nodes, include_virtual=false)
    end
    return FundamentalBessel{typeof(k),eltype(basis.symmetries)}(length(pos),2*pi/k_basis,pos,basis.symmetries)
end


@inline function basis_fun(basis::FundamentalBessel{T}, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let symmetries = basis.symmetries
        pts_loc = [pt - basis.positions[i] for pt in pts]
        arg = [k*hypot(pt[1], pt[2]) for pt in pts_loc]
        b = Y0.(arg)
        if ~isnothing(symmetries)
            for sym in symmetries
                pts_sym = sym.sym_map.(pts)
                pts_loc = [pt - basis.positions[i] for pt in pts_sym]
                arg = [k*hypot(pt[1], pt[2]) for pt in pts_loc]
                b .+= sym.parity.*Y0.(arg)
            end
            b = b./(length(symmetries)+1.0)
        end #norm::T = one(T)/sqrt(basis.dim)
        return b
    end
end

#try using strided for this
@inline function basis_fun(basis::FundamentalBessel{T}, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let symmetries = basis.symmetries
        M =  length(pts)
        N = length(indices)
        B = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            pts_loc = [pt - basis.positions[i] for pt in pts]
            arg = [k*hypot(pt[1], pt[2]) for pt in pts_loc]
            B[:,i] .= Y0.(arg)
            if ~isnothing(symmetries)
                for sym in symmetries
                    pts_sym = sym.sym_map.(pts)
                    pts_loc = [pt - basis.positions[i] for pt in pts_sym]
                    arg = [k*hypot(pt[1], pt[2]) for pt in pts_loc]
                    B[:,i] .+= sym.parity.*Y0.(arg)
                end
            end
        end
        if ~isnothing(symmetries)
            return B./(length(symmetries)+1.0)
        else
            return B
        end 
    end
end

function basis_and_gradient(basis::FundamentalBessel{T}, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let symmetries = basis.symmetries
        pts_loc = [pt - basis.positions[i] for pt in pts]
        r = [hypot(pt[1], pt[2]) for pt in pts_loc]
        x = getindex.(pts_loc, 1)
        y = getindex.(pts_loc, 2) 
        arg = k.*r
        #norm::T = one(T)/sqrt(basis.dim)
        bf = Y0.(arg)
        dx = @. k*x/r*dY0(arg)
        dy = @. k*y/r*dY0(arg)
        if ~isnothing(symmetries)
            n = (length(symmetries)+1.0)
            for sym in symmetries
                pts_sym = sym.sym_map.(pts)
                pts_loc = [pt - basis.positions[i] for pt in pts_sym]
                r = [hypot(pt[1], pt[2]) for pt in pts_loc]
                x = getindex.(pts_loc, 1)
                y = getindex.(pts_loc, 2) 
                arg = k.*r
                bf .+= sym.parity.*Y0.(arg)
                dx .+= @. sym.parity*k*x/r*dY0(arg)
                dy .+= @. sym.parity*k*y/r*dY0(arg)
            end
            bf = bf./n
            dx = dx./n
            dy = dy./n
        end
        return bf, dx, dy
    end
end


function basis_and_gradient(basis::FundamentalBessel{T}, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let symmetries = basis.symmetries
        M =  length(pts)
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
            if ~isnothing(symmetries)
                for sym in symmetries
                    pts_sym = sym.sym_map.(pts)
                    pts_loc = [pt - basis.positions[i] for pt in pts_sym]
                    r = [hypot(pt[1], pt[2]) for pt in pts_loc]
                    x = getindex.(pts_loc, 1)
                    y = getindex.(pts_loc, 2) 
                    arg = k.*r
                    B[:,i] .+= sym.parity.*Y0.(arg)
                    dB_dx[:,i] .+= @. sym.parity*k*x/r*dY0(arg)
                    dB_dy[:,i] .+= @. sym.parity*k*y/r*dY0(arg)
                end
            end
        end
        if ~isnothing(symmetries)
            n = (length(symmetries)+1.0)
            return B./n, dB_dx./n, dB_dy./n
        else
        #println(size(s))
            return B, dB_dx, dB_dy
        end
    end
end

@inline function dk_fun(basis::FundamentalBessel{T}, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let symmetries = basis.symmetries
        pts_loc = [pt - basis.positions[i] for pt in pts]
        r = [hypot(pt[1], pt[2]) for pt in pts_loc]
        dk = @. r*dY0(k*r)
        if ~isnothing(symmetries)
            n = (length(symmetries)+1.0)
            for sym in symmetries
                pts_sym = sym.sym_map.(pts)
                pts_loc = [pt - basis.positions[i] for pt in pts_sym]
                r = [hypot(pt[1], pt[2]) for pt in pts_loc]
                dk .+= @. sym.parity*r*dY0(k*r)
            end
            dk = dk./n
        end
        #norm::T = one(T)/sqrt(basis.dim)
        return dk
    end
end
    

@inline function dk_fun(basis::FundamentalBessel{T}, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let symmetries = basis.symmetries
        M =  length(pts)
        N = length(indices)
        dB_dk = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            pts_loc = [pt - basis.positions[i] for pt in pts]
            r = [hypot(pt[1], pt[2]) for pt in pts_loc]
            dB_dk[:,i] .= @. r*dY0(k*r)
            if ~isnothing(symmetries)
                for sym in symmetries
                    pts_sym = sym.sym_map.(pts)
                    pts_loc = [pt - basis.positions[i] for pt in pts_sym]
                    r = [hypot(pt[1], pt[2]) for pt in pts_loc]
                    dB_dk[:,i] .+= @. sym.parity*r*dY0(k*r)
                end
            end
        end
        if ~isnothing(symmetries)
            n = (length(symmetries)+1.0)
            return dB_dk./n
        else
            return dB_dk
        end
    end
end