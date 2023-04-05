
struct RealPlaneWaves{T,Sy,F} <: AbsBasis where  {T<:Real, Sy<:Union{AbsSymmetry,Nothing}, F<:Function}
    #cs::PolarCS{T} #not fully implemented
    dim::Int64 #using concrete type
    symmetries::Union{Vector{Sy},Nothing}
    angle_arc::T
    angle_shift::T
    angles::Vector{T}
    parity_x::Vector{Int64}
    parity_y::Vector{Int64}
    sampler::F
end

function parity_pattern(symmetries)
    if isnothing(symmetries)
        parity_x = [1,1,-1,-1]
        parity_y = [1,-1,1,-1]
    else
        parity_x = [1,-1]
        parity_y = [1,-1]
        for sym in symmetries
            if sym.axis == :y_axis
                parity_x = [sym.parity, sym.parity]
            end
            if sym.axis == :x_axis
                parity_y = [sym.parity, sym.parity]   
            end
        end
    end            
     
    return parity_x, parity_y
end

function RealPlaneWaves(dim, symmetries; angle_arc = pi, angle_shift=0.0, sampler=linear_nodes)
    par_x, par_y = parity_pattern(symmetries)
    pl = length(par_x)
    eff_dim = dim*pl
    t, dt = sampler(dim)
    angles = @. t*angle_arc + angle_shift
    angles = repeat(angles, inner=pl)
    par_x = repeat(par_x, outer=dim)
    par_y = repeat(par_y, outer=dim)
    return RealPlaneWaves(eff_dim, symmetries, angle_arc, angle_shift, angles, par_x, par_y, sampler)
end

function RealPlaneWaves(dim; angle_arc = pi, angle_shift=0.0, sampler=linear_nodes)
    symmetries = nothing
    par_x, par_y = parity_pattern(symmetries)
    pl = length(par_x)
    eff_dim = dim*pl
    t, dt = sampler(dim)
    angles = @. t*angle_arc + angle_shift
    angles = repeat(angles, inner=pl)
    par_x = repeat(par_x, outer=dim)
    par_y = repeat(par_y, outer=dim)
    return RealPlaneWaves{eltype(angles),Nothing,typeof(sampler)}(eff_dim, symmetries, angle_arc, angle_shift, angles, par_x, par_y, sampler)
end

function resize_basis(basis::Ba, billiard::Bi, dim::Int, k) where {Ba<:RealPlaneWaves,Bi<:AbsBilliard}
    return RealPlaneWaves(dim, basis.symmetries; angle_arc = basis.angle_arc, angle_shift=basis.angle_shift, sampler=basis.sampler)
end

@inline function rpw(arg, parity::Int64)
    if parity == 1
        return cos.(arg)
    else
        return sin.(arg)
    end 
end

@inline function d_rpw(arg, parity::Int64)
    if parity == 1
        return -sin.(arg)
    else
        return cos.(arg)
    end 
end

@inline function basis_fun(basis::RealPlaneWaves, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let par_x = basis.parity_x, par_y = basis.parity_y
        x = getindex.(pts,1)
        y = getindex.(pts,2)
        vx = cos(basis.angles[i])
        vy = sin(basis.angles[i])
        arg_x = k*vx.*x
        arg_y = k*vy.*y
        b = rpw(arg_x, par_x[i]).*rpw(arg_y, par_y[i])
        return b
    end
end

@inline function basis_fun(basis::RealPlaneWaves, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let par_x = basis.parity_x, par_y = basis.parity_y
        x = getindex.(pts,1)
        y = getindex.(pts,2)
        M =  length(pts)
        N = length(indices)
        B = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            vx = cos(basis.angles[i])
            vy = sin(basis.angles[i])
            arg_x = k*vx.*x
            arg_y = k*vy.*y
            B[:,i] .= rpw(arg_x, par_x[i]).*rpw(arg_y, par_y[i])
        end
        return B 
    end
end

function gradient(basis::RealPlaneWaves, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let par_x = basis.parity_x, par_y = basis.parity_y
        x = getindex.(pts,1)
        y = getindex.(pts,2)
        vx = cos(basis.angles[i])
        vy = sin(basis.angles[i])
        arg_x = k*vx.*x
        arg_y = k*vy.*y
        bx = rpw(arg_x, par_x[i])
        by = rpw(arg_y, par_y[i])
        dx = k*vx.*d_rpw(arg_x, par_x[i]).*by
        dy = bx.*k*vy.*d_rpw(arg_y, par_y[i])
        return dx, dy
    end
end

function gradient(basis::RealPlaneWaves, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let par_x = basis.parity_x, par_y = basis.parity_y
        x = getindex.(pts,1)
        y = getindex.(pts,2)
        M =  length(pts)
        N = length(indices)
        dB_dx = zeros(T,M,N)
        dB_dy = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            vx = cos(basis.angles[i])
            vy = sin(basis.angles[i])
            arg_x = k*vx.*x
            arg_y = k*vy.*y
            bx = rpw(arg_x, par_x[i])
            by = rpw(arg_y, par_y[i])
            dB_dx[:,i] .= k*vx.*d_rpw(arg_x, par_x[i]).*by
            dB_dy[:,i] .= bx.*k*vy.*d_rpw(arg_y, par_y[i])
        end
        return dB_dx, dB_dy
    end
end


function basis_and_gradient(basis::RealPlaneWaves, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let par_x = basis.parity_x, par_y = basis.parity_y
        x = getindex.(pts,1)
        y = getindex.(pts,2)
        vx = cos(basis.angles[i])
        vy = sin(basis.angles[i])
        arg_x = k*vx.*x
        arg_y = k*vy.*y
        bx = rpw(arg_x, par_x[i])
        by = rpw(arg_y, par_y[i])
        bf = bx.*by
        dx = k*vx.*d_rpw(arg_x, par_x[i]).*by
        dy = bx.*k*vy.*d_rpw(arg_y, par_y[i])
        return bf, dx, dy
    end
end


function basis_and_gradient(basis::RealPlaneWaves, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let par_x = basis.parity_x, par_y = basis.parity_y
        x = getindex.(pts,1)
        y = getindex.(pts,2)
        M =  length(pts)
        N = length(indices)
        B = zeros(T,M,N)
        dB_dx = zeros(T,M,N)
        dB_dy = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            vx = cos(basis.angles[i])
            vy = sin(basis.angles[i])
            arg_x = k*vx.*x
            arg_y = k*vy.*y
            bx = rpw(arg_x, par_x[i])
            by = rpw(arg_y, par_y[i])
            B[:,i] .= bx.*by
            dB_dx[:,i] .= k*vx.*d_rpw(arg_x, par_x[i]).*by
            dB_dy[:,i] .= bx.*k*vy.*d_rpw(arg_y, par_y[i])

        end
        return B, dB_dx, dB_dy
    end
end

@inline function dk_fun(basis::RealPlaneWaves, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let par_x = basis.parity_x, par_y = basis.parity_y
        x = getindex.(pts,1)
        y = getindex.(pts,2)
        vx = cos(basis.angles[i])
        vy = sin(basis.angles[i])
        arg_x = k*vx.*x
        arg_y = k*vy.*y
        bx = rpw(arg_x, par_x[i])
        by = rpw(arg_y, par_y[i])
        d_bx = d_rpw(arg_x, par_x[i])
        d_by = d_rpw(arg_y, par_y[i])
        dk = @. vx*x*d_bx*by + bx*vy*y*d_by
        return dk
    end
end
    

@inline function dk_fun(basis::RealPlaneWaves, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let par_x = basis.parity_x, par_y = basis.parity_y
        x = getindex.(pts,1)
        y = getindex.(pts,2)
        M =  length(pts)
        N = length(indices)
        dB_dk = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            vx = cos(basis.angles[i])
            vy = sin(basis.angles[i])
            arg_x = k*vx.*x
            arg_y = k*vy.*y
            bx = rpw(arg_x, par_x[i])
            by = rpw(arg_y, par_y[i])
            d_bx = d_rpw(arg_x, par_x[i])
            d_by = d_rpw(arg_y, par_y[i])
            dB_dk[:,i] .=  @. vx*x*d_bx*by + bx*vy*y*d_by
        end
        return dB_dk
    end
end
