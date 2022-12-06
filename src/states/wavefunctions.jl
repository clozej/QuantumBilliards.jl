include("../abstracttypes.jl")
include("../billiards/billiard.jl")
include("../utils/gridutils.jl")
include("../solvers/matrixconstructors.jl")


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

function wavefunction(state::AbsState, basis::AbsBasis, billiard::AbsBilliard; b=5.0, inside_only=true) 
    vec = state.vec
    typ = eltype(vec)
    k = state.k
    dim = state.dim
    new_basis = rescale_basis(basis, dim)
    
    #try to find a lazy way to do this
    L = real_length(billiard)
    xlim,ylim,dx,dy= boundary_limits(billiard.boundary; grd=round(Int, k*L*b/(2*pi)))
    nx = max(round(Int, k*dx*b/(2*pi)), 400)
    ny = max(round(Int, k*dy*b/(2*pi)), 400)
    #println((nx,ny))
    #x, y, x_grid, y_grid = lazy_grid(xlim, ylim, (nx,ny))
    x_plot = range(xlim... , nx)
    y_plot = range(ylim... , ny)
    if typ == Float32
        k = Float32(k)
        x_plot = Float32.(x_plot)
        y_plot = Float32.(y_plot)
    end
    x = repeat(x_plot , outer = length(y_plot))
    y = repeat(y_plot , inner = length(x_plot))
    
    dim_ind = collect(1:length(vec))[vec .!= 0.0] #awkward
    if inside_only
        inside_ind  = collect(1:length(x))[billiard.domain(x,y)]
        #println(inside)
        B = basis_matrix(new_basis,k, x, y, dim_ind, inside_ind) 
    else
        B = basis_matrix(new_basis,k, x, y) 
        #B = zeros(typ,length(x),basis.dim)
        #basis_matrix!(B, basis,k, x, y) 
    end
    
    phi = B * vec
    Psi = reshape(phi, (nx,ny))

    return Psi, x_plot, y_plot
end

function wavefunction(state::AbsState, basis::AbsBasis, curve::AbsCurve; sampler=linear_nodes, b=5.0) 
    vec = state.vec
    typ = eltype(vec)
    k = state.k
    dim = state.dim
    new_basis = rescale_basis(basis, dim)
    
    #try to find a lazy way to do this
    L = curve.length
    N = max(round(Int, k*L*b/(2*pi)), 400)
    t, dt = sampler(N)
    x, y = curve.r(t)
    s = curve.s(t)

    if typ == Float32
        k = Float32(k)
        x = Float32.(x)
        y = Float32.(y)
        s = Float32.(s)
    end
    
    dim_ind = collect(1:length(vec))[vec .!= 0.0] #awkward
    B = basis_matrix(new_basis,k, x, y, dim_ind) 
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
    #new_basis = rescale_basis(basis, dim)
    
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
