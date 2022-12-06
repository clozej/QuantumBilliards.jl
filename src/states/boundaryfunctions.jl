include("../abstracttypes.jl")
include("../billiards/billiard.jl")
#include("../utils/gridutils.jl")
include("../solvers/matrixconstructors.jl")
#=
struct BoundaryPoints{T} <: AbsPoints
    x :: Vector{T} 
    y :: Vector{T} 
    nx :: Vector{T} 
    ny :: Vector{T}
    s :: Vector{T}
    w :: Vector{T}
end
=#
function boundary_coords(curve::AbsCurve, N; sampler=linear_nodes)
    L = curve.length
    t, dt = sampler(N)
    x, y = curve.r(t)
    nx, ny = curve.n(t)
    s = curve.s(t)
    ds = L*dt #modify for different parametrizations
    return x,y,nx,ny,s,ds
end    

function boundary_coords(billiard::AbsBilliard, N; sampler=linear_nodes, include_virtual=true)
    x_all = Float64[]
    y_all = Float64[]
    nx_all = Float64[]
    ny_all = Float64[]
    s_all = Float64[]
    ds_all = Float64[]

    L = real_length(billiard)
    if include_virtual
       L += virtual_length(billiard) 
    end
    l = 0.0 #cumulative length
    for curve in billiard.boundary
        if (typeof(curve) <: AbsRealCurve || include_virtual)
            Lc = curve.length
            Nc = round(Int, N*Lc/L)
            x,y,nx,ny,s,ds = boundary_coords(curve, Nc; sampler=sampler)
            append!(x_all, x)
            append!(y_all, y)
            append!(nx_all, nx)
            append!(ny_all, ny)
            append!(s_all, s .+ l)
            append!(ds_all, ds)
            l += Lc
        end    
    end
    return x_all,y_all,nx_all,ny_all,s_all,ds_all 
end   

function boundary_function(state::AbsState, basis::AbsBasis, billiard::AbsBilliard; b=5.0, sampler=linear_nodes, include_virtual=true)
    vec = state.vec
    #typ = eltype(vec)
    new_basis = rescale_basis(basis,state.dim)
    k = state.k
    L = real_length(billiard)
    N = max(round(Int, k*L*b/(2*pi)), 400)
    x,y,nx,ny,s,ds = boundary_coords(billiard, N; sampler=sampler, include_virtual=include_virtual)
    U = U_matrix(new_basis,k,x,y,nx,ny)
    u = U * vec
    return s, u, U
end