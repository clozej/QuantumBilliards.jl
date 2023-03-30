using StaticArrays

#try using strided to optimize this
function compute_psi(state::S, x_grid, y_grid; inside_only=true, memory_limit = 10.0e9) where {S<:AbsState}
    let vec = state.vec, k = state.k_basis, basis=state.basis, billiard=state.billiard, eps=state.eps #basis is correct size
        sz = length(x_grid)*length(y_grid)
        pts = collect(SVector(x,y) for y in y_grid for x in x_grid)
        if inside_only
            pts_mask = is_inside(billiard,pts)
            pts = pts[pts_mask]
        end
        n_pts = length(pts)
        #estimate max memory needed for the matrices
        type = eltype(vec)
        memory = sizeof(type)*basis.dim*n_pts
        Psi = zeros(type,sz)

        if memory < memory_limit
            B = basis_matrix(basis, k, pts)
            Psi_pts = B*vec
            if inside_only
                Psi[pts_mask] .= Psi_pts
            else
                Psi .= Psi_pts
            end
        else
            println("Warning: memory limit of $(Base.format_bytes(memory_limit)) exceded $(Base.format_bytes(memory)).")
            if inside_only
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
        end
        return Psi
    end
end

function wavefunction(state::S; b=5.0, inside_only=true, memory_limit = 10.0e9) where {S<:AbsState}
    let k = state.k, billiard=state.billiard       
        #println(new_basis.dim)
        type = eltype(state.vec)
        #try to find a lazy way to do this
        L = billiard.length
        xlim,ylim = boundary_limits(billiard.fundamental_boundary; grd=max(1000,round(Int, k*L*b/(2*pi))), type=type)
        dx = xlim[2] - xlim[1]
        dy = ylim[2] - ylim[1]
        nx = max(round(Int, k*dx*b/(2*pi)), 512)
        ny = max(round(Int, k*dy*b/(2*pi)), 512)
        x_grid::Vector{type} = collect(type,range(xlim... , nx))
        y_grid::Vector{type} = collect(type,range(ylim... , ny))
        Psi::Vector{type} = compute_psi(state,x_grid,y_grid;inside_only=inside_only, memory_limit = memory_limit) 
        #println("Psi type $(eltype(Psi)), $(memory_size(Psi))")
        Psi2d::Array{type,2} = reshape(Psi, (nx,ny))
        return Psi2d, x_grid, y_grid
    end
end

function wavefunction(state::BasisState; xlim =(-3.0,3.0), ylim=(-3.0,3.0), b=5.0, inside_only=true, memory_limit = 10.0e9) 
    let k = state.k, basis=state.basis      
        #println(new_basis.dim)
        type = eltype(state.vec)
        #try to find a lazy way to do this
        dx = xlim[2] - xlim[1]
        dy = ylim[2] - ylim[1]
        nx = max(round(Int, k*dx*b/(2*pi)), 512)
        ny = max(round(Int, k*dy*b/(2*pi)), 512)
        x_grid::Vector{type} = collect(type,range(xlim... , nx))
        y_grid::Vector{type} = collect(type,range(ylim... , ny))
        pts_grid = [SVector(x,y) for y in y_grid for x in x_grid]
        Psi::Vector{type} = basis_fun(basis,state.idx,k,pts_grid) 
        #println("Psi type $(eltype(Psi)), $(memory_size(Psi))")
        Psi2d::Array{type,2} = reshape(Psi, (nx,ny))
        return Psi2d, x_grid, y_grid
    end
end

#this can be optimized
function compute_psi(state_bundle::S, x_grid, y_grid; inside_only=true, memory_limit = 10.0e9) where {S<:EigenstateBundle}
    let k = state_bundle.k_basis, basis=state_bundle.basis, billiard=state_bundle.billiard, X=state_bundle.X #basis is correct size
        sz = length(x_grid)*length(y_grid)
        pts = collect(SVector(x,y) for y in y_grid for x in x_grid)
        if inside_only
            pts_mask = is_inside(billiard,pts)
            pts = pts[pts_mask]
        end
        n_pts = length(pts)
        n_states = length(state_bundle.ks)
        #estimate max memory needed for the matrices
        type = eltype(state_bundle.X)
        memory = sizeof(type)*basis.dim*n_pts
        #Vector of results
        Psi_bundle = zeros(type,(sz,n_states))    
        if memory < memory_limit
            #Psi = zeros(type,sz)
            B = basis_matrix(basis, k, pts)
            Psi_pts = B*X
            Psi_bundle[pts_mask,:] .= Psi_pts
        else
            println("Warning: memory limit of $(Base.format_bytes(memory_limit)) exceded $(Base.format_bytes(memory)).")
            
            Psi_pts = zeros(type,(n_pts,n_states))
            for i in 1:basis.dim
                bf = basis_fun(basis,i,k,pts) #vector of length n_pts
                for j in 1:n_states
                    Psi_pts[:,j] .+= X[i,j].*bf
                end
            end
            if inside_only
                Psi_bundle[pts_mask,:] = Psi_pts
            else
                Psi_bundle = Psi_pts
            end
            #println("Psi type $(eltype(Psi)), $(memory_size(Psi))")
        end
        return Psi_bundle #this is a matrix 
    end
end

function wavefunction(state_bundle::S; b=5.0, inside_only=true, memory_limit = 10.0e9) where {S<:EigenstateBundle}
    let k = state_bundle.k_basis, billiard=state_bundle.billiard       
        #println(new_basis.dim)
        type = eltype(state_bundle.X)
        #try to find a lazy way to do this
        L = billiard.length
        xlim,ylim = boundary_limits(billiard.fundamental_boundary; grd=max(1000,round(Int, k*L*b/(2*pi))), type=type)
        dx = xlim[2] - xlim[1]
        dy = ylim[2] - ylim[1]
        nx = max(round(Int, k*dx*b/(2*pi)), 512)
        ny = max(round(Int, k*dy*b/(2*pi)), 512)
        x_grid::Vector{type} = collect(type,range(xlim... , nx))
        y_grid::Vector{type} = collect(type,range(ylim... , ny))
        Psi_bundle::Matrix{type} = compute_psi(state_bundle,x_grid,y_grid;inside_only=inside_only, memory_limit = memory_limit) 
        #println("Psi type $(eltype(Psi)), $(memory_size(Psi))")
        Psi2d::Vector{Array{type,2}} = [reshape(Psi, (nx,ny)) for Psi in eachcol(Psi_bundle)]
        return Psi2d, x_grid, y_grid
    end
end




#not finished
#=
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

=#


#=
function wavefunction(state::AbsState, basis::AbsBasis, billiard::AbsBilliard; b=5.0, inside_only=true) 
    vec = state.vec
    typ = eltype(vec)
    k = state.k
    dim = state.dim
    #new_basis = resize_basis(basis, dim)
    
    #try to find a lazy way to do this
    L = billiard.length
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
