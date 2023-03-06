#include("../abstracttypes.jl")
#include("../utils/billiardutils.jl")
#include("../utils/gridutils.jl")
#include("../utils/benchmarkutils.jl")
#include("../solvers/matrixconstructors.jl")


#=
function basis_matrix(basis::AbsBasis,k,vec, x::Vector{T}, y::Vector{T}) where T<:Number
    M =  length(x)
    N = basis.dim
    B = Array{T}(0.0,M,N)  #basis matrix
    
    idx = collect(1:N)[abs.(vec).!=0.0]
    @inbounds Threads.@threads for i in idx
        B[:,i] = basis_fun(basis, i, k, x, y)
    end 
    return B
end
=#
#try using strided to optimize this
function compute_psi(state::S, basis::Ba, billiard::Bi, x_grid, y_grid; inside_only=true) where {S<:AbsState,Ba<:AbsBasis,Bi<:AbsBilliard}
    let vec = state.vec, k = state.k_basis, basis=basis, eps=state.eps
        #sz = length(x_grid)*length(y_grid)
        pts = collect(SVector(x,y) for y in y_grid for x in x_grid)
        #println("Points type $(eltype(pts)), $(memory_size(pts))")
        Psi = zeros(eltype(vec),length(pts))
        if inside_only
            pts_mask = is_inside(billiard,x_grid,y_grid)
            pts = pts[pts_mask]
            for i in eachindex(vec)
                if abs(vec[i]) > eps 
                    Psi[pts_mask] .+= vec[i].*basis_fun(basis,i,k,pts)
                end
            end
        else
            for i in eachindex(vec)
                if abs(vec[i]) > eps 
                    Psi .+= vec[i].*basis_fun(basis,i,k,pts)
                end
            end
        end
        #println("Psi type $(eltype(Psi)), $(memory_size(Psi))")
        return Psi
    end
end



function wavefunction(state::S, basis::Ba, billiard::Bi; b=5.0, inside_only=true) where {S<:AbsState,Ba<:AbsBasis,Bi<:AbsBilliard}
    let new_basis = resize_basis(basis, state.dim) 
        k = state.k       
        #println(new_basis.dim)
        type = eltype(state.vec)
        #try to find a lazy way to do this
        L = real_length(billiard)
        xlim,ylim = boundary_limits(billiard.boundary; grd=round(Int, k*L*b/(2*pi)), type=type)
        dx = xlim[2] - xlim[1]
        dy = ylim[2] - ylim[1]
        nx = max(round(Int, k*dx*b/(2*pi)), 512)
        ny = max(round(Int, k*dy*b/(2*pi)), 512)
        x_grid::Vector{type} = collect(type,range(xlim... , nx))
        y_grid::Vector{type} = collect(type,range(ylim... , ny))
        #dim_ind = collect(1:length(vec))[vec .!= zero(type)] #awkward
        #=
        if inside_only
            pts_mask = collect(is_inside(billiard,SVector(x,y)) for y in y_grid for x in x_grid)#collect(Iterators.flatten(gen))
            B = basis_matrix(billiard, new_basis, k, x_grid, y_grid, dim_ind)
        else
            B = basis_matrix(billiard, new_basis, k, x_grid, y_grid, dim_ind)
        end
        =#
        #B = basis_matrix(new_basis, k, x_grid, y_grid)#, dim_ind)
        #println("B type $(eltype(B)), $(memory_size(B))")
        #Psi = B * vec
        Psi::Vector{type} = compute_psi(state,new_basis,billiard,x_grid,y_grid;inside_only=inside_only) 
        #println("Psi type $(eltype(Psi)), $(memory_size(Psi))")
        Psi2d::Array{type,2} = reshape(Psi, (nx,ny))
        return Psi2d , x_grid, y_grid
    end
    
end


#not finished
function wavefunction(state::S, basis::Ba, crv::C; sampler=linear_nodes, b=5.0) where {S<:AbsState,Ba<:AbsBasis,C<:AbsCurve}
    vec = state.vec
    #typ = eltype(vec)
    k = state.k
    dim = state.dim
    new_basis = resize_basis(basis, dim)
    
    #try to find a lazy way to do this
    L = crv.length
    N = max(round(Int, k*L*b/(2*pi)), 512)
    t, dt = sampler(N)
    pts = curve(crv,t)
    s = arc_length(crv,t)
    
    dim_ind = collect(1:length(vec))[vec .!= zero(type)] #awkward
    B = basis_matrix(new_basis,k,pts, dim_ind) 
    #B = zeros(typ,length(x),basis.dim)
    #basis_matrix!(B, basis,k, x, y) 
    phi = B * vec
    return phi, s
end




#=
function wavefunction(state::AbsState, basis::AbsBasis, billiard::AbsBilliard; b=5.0, inside_only=true) 
    vec = state.vec
    typ = eltype(vec)
    k = state.k
    dim = state.dim
    #new_basis = resize_basis(basis, dim)
    
    #try to find a lazy way to do this
    L = real_length(billiard)
    xlim,ylim,dx,dy= boundary_limits(billiard.boundary; grd=round(Int, k*L*b/(2*pi)))
    nx = max(round(Int, k*dx*b/(2*pi)), 400)
    ny = max(round(Int, k*dy*b/(2*pi)), 400)
    #println((nx,ny))
    x_grid = collect(range(xlim..., length=nx))
    y_grid = collect(range(ylim..., length=ny))
    
    Psi = zeros(typ,ny,nx)
    #println(size(Psi))
    idx  = collect(1:dim)[abs.(vec).!=0.0] #only compute basis states that have sufficiently large coficients
    if inside_only
        for (i,x) in enumerate(x_grid)
            psi_col = zeros(typ,ny)#wavefunction for all y at selected x
            x_temp = x.*ones(typ,ny)
            inside = billiard.domain(x_temp, y_grid) 
            for n in idx
                bf = zeros(typ,ny)
                bf[inside] .= basis_fun(basis, n, k, x_temp[inside], y_grid[inside])
                psi_col .+= vec[n].*bf
            end
            #println(size(psi_col))
            Psi[:,i] .= psi_col
        end
    else # computing outside wavefunction
        for (i,x) in enumerate(x_grid)
            psi_col = zeros(typ,ny)#wavefunction for all y at selected x
            x_temp = x.*ones(typ,ny)
            for n in idx
                bf = basis_fun(basis, n, k, x_temp, y_grid)
                psi_col .+= vec[n].*bf
            end
            Psi[:,i] .= psi_col
        end
    end
    
    return Psi, x_grid, y_grid
end
=#
