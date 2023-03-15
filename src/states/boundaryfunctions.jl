#include("../abstracttypes.jl")
#include("../utils/billiardutils.jl")
#include("../utils/gridutils.jl")
#include("../solvers/matrixconstructors.jl")
using FFTW

#this takes care of singular points
function regularize!(u)
    idx = findall(isnan, u)
    for i in idx
        u[i] = (u[i+1] + u[i-1])/2.0
    end
end

function boundary_function(state::S; b=5.0, sampler=fourier_nodes, include_virtual=true) where {S<:AbsState}
    let vec = state.vec, k = state.k, k_basis = state.k_basis, new_basis = state.basis, billiard=state.billiard
        type = eltype(vec)
        L = real_length(billiard)
        N = max(round(Int, k*L*b/(2*pi)), 512)
        pts = boundary_coords(billiard, N; sampler=sampler, include_virtual=include_virtual)
        dX, dY = gradient_matrices(new_basis, k_basis, pts.xy)
        nx = getindex.(pts.normal,1)
        ny = getindex.(pts.normal,2)
        dX = nx .* dX 
        dY = ny .* dY
        U::Array{type,2} = dX .+ dY
        u::Vector{type} = U * vec
        regularize!(u)
        #compute the boundary norm
        w = dot.(pts.normal, pts.xy) .* pts.ds
        integrand = abs2.(u) .* w
        norm = sum(integrand)/(2*k^2)
        #println(norm)
        return u, pts.s::Vector{type}, norm
    end
end

function boundary_function(state_bundle::S; b=5.0, sampler=fourier_nodes, include_virtual=true) where {S<:EigenstateBundle}
    let X = state_bundle.X, k_basis = state_bundle.k_basis, ks = state_bundle.ks, new_basis = state_bundle.basis, billiard=state_bundle.billiard 
        type = eltype(X)
        L = real_length(billiard)
        N = max(round(Int, k_basis*L*b/(2*pi)), 512)
        pts = boundary_coords(billiard, N; sampler=sampler, include_virtual=include_virtual)
        dX, dY = gradient_matrices(new_basis, k_basis, pts.xy)
        nx = getindex.(pts.normal,1)
        ny = getindex.(pts.normal,2)
        dX = nx .* dX 
        dY = ny .* dY
        U::Array{type,2} = dX .+ dY
        u_bundle::Matrix{type} = U * X
        for u in eachcol(u_bundle)
            regularize!(u)
        end
        #compute the boundary norm
        w = dot.(pts.normal, pts.xy) .* pts.ds
        norms = [sum(abs2.(u_bundle[:,i]) .* w)/(2*ks[i]^2) for i in eachindex(ks)]
        #println(norm)
        us::Vector{Vector{type}} = [u for u in eachcol(u_bundle)]
        return us, pts.s::Vector{type}, norms
    end
end

function momentum_function(u,s)
    fu = rfft(u)
    sr = 1.0/diff(s)[1]
    ks = rfftfreq(length(s),sr).*(2*pi)
    return abs2.(fu)/length(fu), ks
end

function momentum_function(state::S; b=5.0, sampler=fourier_nodes, include_virtual=true) where {S<:AbsState}
    u, s, norm = boundary_function(state; b=b, sampler=sampler, include_virtual=include_virtual)
    return momentum_function(u,s)
end

#this can be optimized by usinf FFTW plans
function momentum_function(state_bundle::S; b=5.0, sampler=fourier_nodes, include_virtual=true) where {S<:EigenstateBundle}
    us, s, norms = boundary_function(state_bundle; b=b, sampler=sampler, include_virtual=include_virtual)
    mf, ks = momentum_function(us[1],s)
    type = eltype(mf)
    mfs::Vector{Vector{type}} = [mf]
    for i in 2:length(us)
        mf, ks = momentum_function(us[i],s)
        push!(mfs,mf)
    end
    return mfs, ks
end