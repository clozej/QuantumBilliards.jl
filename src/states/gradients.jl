#include("../abstracttypes.jl")
#include("../utils/billiardutils.jl")
#include("../utils/gridutils.jl")
#include("../utils/benchmarkutils.jl")
#include("../solvers/matrixconstructors.jl")

#for checking correctness
function compute_grad(state::S, basis::Ba, billiard::Bi,x_grid,y_grid ;inside_only=true) where {S<:AbsState,Ba<:AbsBasis,Bi<:AbsBilliard}
    let vec = state.vec, k = state.k, basis=basis, eps=state.eps
        #sz = length(x_grid)*length(y_grid)
        pts = collect(SVector(x,y) for y in y_grid for x in x_grid)
        #println("Points type $(eltype(pts)), $(memory_size(pts))")
        dX = zeros(eltype(vec),length(pts))
        dY = zeros(eltype(vec),length(pts))
        if inside_only
            pts_mask = is_inside(billiard,x_grid,y_grid)
            pts = pts[pts_mask]
            for i in eachindex(vec)
                if abs(vec[i]) > eps 
                    dx,dy = gradient(basis,i,k,pts)
                    dX[pts_mask] .+= vec[i].*dx
                    dY[pts_mask] .+= vec[i].*dy
                end
            end
        else
            for i in eachindex(vec)
                if abs(vec[i]) > eps 
                    dx,dy = gradient(basis,i,k,pts)
                    dX .+= vec[i].*dx
                    dY .+= vec[i].*dy
                end
            end
        end
        #println("Psi type $(eltype(Psi)), $(memory_size(Psi))")
        return dX, dY
    end
end



function wavefunction_gradient(state::S, basis::Ba, billiard::Bi; b=20.0, inside_only=true) where {S<:AbsState,Ba<:AbsBasis,Bi<:AbsBilliard}
    let new_basis = resize_basis(basis, state.dim)        
        #println(new_basis.dim)
        type = eltype(state.vec)
        k = state.k
        #try to find a lazy way to do this
        L = real_length(billiard)
        xlim,ylim = boundary_limits(billiard.boundary; grd=round(Int, k*L*b/(2*pi)), type=type)
        dx = xlim[2] - xlim[1]
        dy = ylim[2] - ylim[1]
        nx = max(round(Int, b),10)
        ny = max(round(Int, b*dy/dx),10)
        println(nx,ny)
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
        dX::Vector{type}, dY::Vector{type} = compute_grad(state,new_basis,billiard,x_grid,y_grid;inside_only=inside_only) 
        #println("Psi type $(eltype(Psi)), $(memory_size(Psi))")
        dX2d::Array{type,2} = reshape(dX, (nx,ny))
        dY2d::Array{type,2} = reshape(dY, (nx,ny))
        return dX2d, dY2d, x_grid, y_grid
    end
    
end

function wavefunction_gradient(state::S, basis::Ba; b=20.0,xlim=(-1,1),ylim=(-1,1),inside_only=true) where {S<:AbsState,Ba<:AbsBasis}
    let new_basis = resize_basis(basis, state.dim)        
        #println(new_basis.dim)
        type = eltype(state.vec)
        k = state.k
        #try to find a lazy way to do this
        L = real_length(billiard)
        #xlim,ylim = boundary_limits(billiard.boundary; grd=round(Int, k*L*b/(2*pi)), type=type)
        dx = xlim[2] - xlim[1]
        dy = ylim[2] - ylim[1]
        nx = max(round(Int, b),10)
        ny = max(round(Int, b*dy/dx),10)
        println(nx,ny)
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
        dX::Vector{type}, dY::Vector{type} = compute_grad(state,new_basis,billiard,x_grid,y_grid;inside_only=inside_only) 
        #println("Psi type $(eltype(Psi)), $(memory_size(Psi))")
        dX2d::Array{type,2} = reshape(dX, (nx,ny))
        dY2d::Array{type,2} = reshape(dY, (nx,ny))
        return dX2d, dY2d, x_grid, y_grid
    end
    
end