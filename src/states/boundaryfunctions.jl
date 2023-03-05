#include("../abstracttypes.jl")
#include("../utils/billiardutils.jl")
#include("../utils/gridutils.jl")
#include("../solvers/matrixconstructors.jl")
using FFTW

struct BoundaryPointsU{T} <: AbsPoints where {T<:Real}
    xy::Vector{SVector{2,T}}
    normal::Vector{SVector{2,T}} #normal vectors in points
    s::Vector{T} # arc length coords
    ds::Vector{T} #integration weights
end

function boundary_coords(crv::C, N; sampler=fourier_nodes) where {C<:AbsCurve}
    L = crv.length
    t, dt = sampler(N)
    xy = curve(crv, t)
    normal = normal_vec(crv,t)
    s = arc_length(crv,t)
    ds = L.*dt #modify for different parametrizations
    return xy, normal, s, ds
end

function boundary_coords(crv::C, t, dt) where {C<:AbsCurve}
    L = crv.length
    xy = curve(crv, t)
    normal = normal_vec(crv,t)
    s = arc_length(crv,t)
    ds = L.*dt #modify for different parametrizations
    return xy, normal, s, ds
end 

#make better
function boundary_coords(billiard::Bi, N; sampler=fourier_nodes, include_virtual=true) where {Bi<:AbsBilliard}
    let boundary = billiard.boundary
        if sampler==fourier_nodes
            crv_lengths = [crv.length for crv in boundary]
            ts, dts = fourier_nodes(N, crv_lengths)

            xy_all, normal_all, s_all, ds_all = boundary_coords(boundary[1], ts[1], dts[1])
            l = boundary[1].length
            for i in 2:length(ts)
                crv = boundary[i]
                if (typeof(crv) <: AbsRealCurve || include_virtual)
                    Lc = crv.length
                    #Nc = round(Int, N*Lc/L)
                    xy,nxy,s,ds = boundary_coords(crv, ts[i], dts[i])
                    append!(xy_all, xy)
                    append!(normal_all, nxy)
                    s = s .+ l
                    append!(s_all, s)
                    append!(ds_all, ds)
                    l += Lc
                end    
            end
        else
            L = real_length(billiard)
            if include_virtual
            L += virtual_length(billiard) 
            end
            Lc = boundary[1].length
            Nc = round(Int, N*Lc/L)
            xy_all, normal_all, s_all, ds_all = boundary_coords(boundary[1], Nc; sampler=sampler)
            #println(s_all)
            l = boundary[1].length #cumulative length
            for crv in boundary[2:end]
                if (typeof(crv) <: AbsRealCurve || include_virtual)
                    Lc = crv.length
                    Nc = round(Int, N*Lc/L)
                    xy,nxy,s,ds = boundary_coords(crv, Nc; sampler=sampler)
                    append!(xy_all, xy)
                    append!(normal_all, nxy)
                    s = s .+ l
                    append!(s_all, s)
                    append!(ds_all, ds)
                    l += Lc
                end    
            end
        end
        return BoundaryPointsU(xy_all,normal_all,s_all,ds_all) 
    end
end

#this takes care of singular points
function regularize!(u)
    idx = findall(isnan, u)
    for i in idx
        u[i] = (u[i+1] + u[i-1])/2.0
    end
end

function boundary_function(state::S, basis::Ba, billiard::Bi; b=5.0, sampler=fourier_nodes, include_virtual=true) where {S<:AbsState,Ba<:AbsBasis,Bi<:AbsBilliard}
    let vec = state.vec, new_basis = resize_basis(basis,state.dim), k = state.k
        type = eltype(vec)
        L = real_length(billiard)
        N = max(round(Int, k*L*b/(2*pi)), 512)
        pts = boundary_coords(billiard, N; sampler=sampler, include_virtual=include_virtual)
        dX, dY = gradient_matrices(new_basis, k, pts.xy)
        nx = getindex.(pts.normal,1)
        ny = getindex.(pts.normal,2)
        dX = nx .* dX 
        dY = ny .* dY
        U::Array{type,2} = dX .+ dY
        u::Vector{type} = U * vec
        regularize!(u)
        w = dot.(pts.normal, pts.xy) .* pts.ds
        integrand = abs2.(u) .* w
        norm = sum(integrand)/(2*state.k^2)
        #println(norm)
        return u, pts.s::Vector{type}, norm
    end
end

function momentum_function(u,s)
    fu = rfft(u)
    sr = 1.0/diff(s)[1]
    ks = rfftfreq(length(s),sr).*(2*pi)
    return fu, ks
end

function momentum_function(state::S, basis::Ba, billiard::Bi; b=5.0, sampler=fourier_nodes, include_virtual=true) where {S<:AbsState,Ba<:AbsBasis,Bi<:AbsBilliard}
    u, s, norm = boundary_function(state, basis, billiard; b=b, sampler=sampler, include_virtual=include_virtual)
    return momentum_function(u,s)
end