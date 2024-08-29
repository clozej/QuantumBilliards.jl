using StaticArrays

#try using strided to optimize this
"""
Compute the wavefunction values on a grid for a given state.

# Description
This function computes the wavefunction values on a specified grid of points based on the provided state. It efficiently handles large matrices by checking memory usage and adjusts computation accordingly.

# Arguments
- `state::S`: The quantum state for which the wavefunction is computed, where `S` is a subtype of `AbsState`.
- `x_grid`: A vector of x-coordinates.
- `y_grid`: A vector of y-coordinates.
- `inside_only::Bool=true`: If `true`, computes the wavefunction only for points inside the billiard.
- `memory_limit`: The maximum memory allowed for computation in bytes. The default value is 10.0e9 bytes.

# Returns
- `Psi`: A vector containing the computed wavefunction values for the grid points.
"""
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

"""
Compute the wavefunction over a 2D grid for a given state.

# Description
This function computes the wavefunction on a specified grid within the billiard. The grid resolution is automatically determined based on the wavenumber `k` and the geometry of the billiard. The function also supports reflection symmetries to extend the computed wavefunction to the full domain.

# Arguments
- `state::S`: The quantum state for which the wavefunction is computed, where `S` is a subtype of `AbsState`.
- `b::Float64=5.0`: A parameter controlling the grid resolution.
- `inside_only::Bool=true`: If `true`, computes the wavefunction only for points inside the billiard.
- `fundamental_domain::Bool=true`: If `true`, computes the wavefunction only in the fundamental domain; otherwise, it includes reflections.
- `memory_limit`: The maximum memory allowed for computation in bytes. The default value is 10.0e9 bytes.

# Returns
- `Psi2d`: A 2D array containing the computed wavefunction values.
- `x_grid`: The x-coordinates of the grid.
- `y_grid`: The y-coordinates of the grid.
"""
function wavefunction(state::S; b=5.0, inside_only=true, fundamental_domain = true, memory_limit = 10.0e9) where {S<:AbsState}
    let k = state.k, billiard=state.billiard, symmetries=state.basis.symmetries       
        #println(new_basis.dim)
        type = eltype(state.vec)
        #try to find a lazy way to do this
        L = billiard.length
        
        xlim,ylim = boundary_limits(billiard.fundamental_boundary; grd=max(1000,round(Int, k*L*b/(2*pi))))
        dx = xlim[2] - xlim[1]
        dy = ylim[2] - ylim[1]
        nx = max(round(Int, k*dx*b/(2*pi)), 512)
        ny = max(round(Int, k*dy*b/(2*pi)), 512)
        x_grid::Vector{type} = collect(type,range(xlim... , nx))
        y_grid::Vector{type} = collect(type,range(ylim... , ny))
        Psi::Vector{type} = compute_psi(state,x_grid,y_grid;inside_only=inside_only, memory_limit = memory_limit) 
        #println("Psi type $(eltype(Psi)), $(memory_size(Psi))")
        Psi2d::Array{type,2} = reshape(Psi, (nx,ny))
        if ~fundamental_domain 
            if ~isnothing(symmetries)
                Psi2d, x_grid, y_grid = reflect_wavefunction(Psi2d,x_grid,y_grid,symmetries)
            end
        end
        return Psi2d, x_grid, y_grid
    end
end

"""
Compute the wavefunction over a 2D grid for a basis state.

# Description
This function computes the wavefunction on a specified grid within a predefined Cartesian grid. The grid resolution is automatically determined based on the wavenumber `k` and the grid limits.

# Arguments
- `state::BasisState`: The basis state for which the wavefunction is computed.
- `xlim`: A tuple specifying the limits of the x-axis.
- `ylim`: A tuple specifying the limits of the y-axis.
- `b::Float64=5.0`: A parameter controlling the grid resolution.

# Returns
- `Psi2d`: A 2D array containing the computed wavefunction values.
- `x_grid`: The x-coordinates of the grid.
- `y_grid`: The y-coordinates of the grid.
"""
function wavefunction(state::BasisState; xlim =(-2.0,2.0), ylim=(-2.0,2.0), b=5.0) 
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
"""
Compute the wavefunction values on a grid for a bundle of eigenstates.

# Description
This function computes the wavefunction values for a bundle of eigenstates on a specified grid of points. It efficiently handles large matrices by checking memory usage and adjusts computation accordingly.

# Arguments
- `state_bundle::S`: The eigenstate bundle for which the wavefunction is computed, where `S` is a subtype of `EigenstateBundle`.
- `x_grid`: A vector of x-coordinates.
- `y_grid`: A vector of y-coordinates.
- `inside_only::Bool=true`: If `true`, computes the wavefunction only for points inside the billiard.
- `memory_limit`: The maximum memory allowed for computation in bytes. The default value is 10.0e9 bytes.

# Returns
- `Psi_bundle`: A matrix where each column represents the computed wavefunction values for a different eigenstate in the bundle.
"""
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

"""
Compute the wavefunction over a 2D grid for a bundle of eigenstates.

# Description
This function computes the wavefunctions for a bundle of eigenstates on a specified grid within the billiard. The grid resolution is automatically determined based on the wavenumber `k` and the geometry of the billiard. The function also supports reflection symmetries to extend the computed wavefunction to the full domain.

# Arguments
- `state_bundle::S`: The eigenstate bundle for which the wavefunction is computed, where `S` is a subtype of `EigenstateBundle`.
- `b::Float64=5.0`: A parameter controlling the grid resolution.
- `inside_only::Bool=true`: If `true`, computes the wavefunction only for points inside the billiard.
- `fundamental_domain::Bool=true`: If `true`, computes the wavefunction only in the fundamental domain; otherwise, it includes reflections.
- `memory_limit`: The maximum memory allowed for computation in bytes. The default value is 10.0e9 bytes.

# Returns
- `Psi2d`: A vector of 2D arrays, each representing the wavefunction for a different eigenstate.
- `x_grid`: The x-coordinates of the grid.
- `y_grid`: The y-coordinates of the grid.
"""
function wavefunction(state_bundle::S; b=5.0, inside_only=true, fundamental_domain = true, memory_limit = 10.0e9) where {S<:EigenstateBundle}
    let k = state_bundle.k_basis, billiard=state_bundle.billiard       
        #println(new_basis.dim)
        type = eltype(state_bundle.X)
        #try to find a lazy way to do this
        L = billiard.length
        xlim,ylim = boundary_limits(billiard.fundamental_boundary; grd=max(1000,round(Int, k*L*b/(2*pi))))
        dx = xlim[2] - xlim[1]
        dy = ylim[2] - ylim[1]
        nx = max(round(Int, k*dx*b/(2*pi)), 512)
        ny = max(round(Int, k*dy*b/(2*pi)), 512)
        x_grid::Vector{type} = collect(type,range(xlim... , nx))
        y_grid::Vector{type} = collect(type,range(ylim... , ny))
        Psi_bundle::Matrix{type} = compute_psi(state_bundle,x_grid,y_grid;inside_only=inside_only, memory_limit = memory_limit) 
        #println("Psi type $(eltype(Psi)), $(memory_size(Psi))")
        Psi2d::Vector{Array{type,2}} = [reshape(Psi, (nx,ny)) for Psi in eachcol(Psi_bundle)]
        if ~fundamental_domain 
            if ~isnothing(symmetries)
                for i in eachindex(Psi2d)
                    Psi_new, x_grid, y_grid = reflect_wavefunction(Psi2d[i],x_grid,y_grid,symmetries)
                    Psi2d[i] = Psi_new
                end
            end
        end
        return Psi2d, x_grid, y_grid
    end
end


