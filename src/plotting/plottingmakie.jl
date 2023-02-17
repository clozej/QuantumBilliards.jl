include("../abstracttypes.jl")
include("../utils/gridutils.jl")
using Makie
#helper functions
function plot_heatmap!(f,x,y,Z ;vmax = 1.0,log=(false,-5.0), cmap=Reverse(:gist_heat),hmargs=Dict(),axargs=Dict())
    if log[1]
        X = log10.(Z)
        ax = Axis(f[1,1],axargs...)        
        m = findmax(X)[1]
        range_val = (log[2],m*vmax)
        hmap = heatmap!(ax,x, y, X, colormap = cmap, colorrange=range_val, hmargs...)
        ax.aspect=DataAspect()
        Colorbar(f[1,2], colormap = cmap, limits = Float64.(range_val),tellheight=true)
        rowsize!(f.layout, 1, ax.scene.px_area[].widths[2])
    else
        ax = Axis(f[1,1],axargs...)        
        m = findmax(Z)[1]
        range_val = (0,m*vmax)
        hmap = heatmap!(ax,x, y, Z, colormap = cmap, colorrange=range_val, hmargs...)
        ax.aspect=DataAspect()
        Colorbar(f[1,2], colormap = cmap, limits = Float64.(range_val),tellheight=true)
        rowsize!(f.layout, 1, ax.scene.px_area[].widths[2])
    end
    return hmap, ax
end

function plot_heatmap_balaced!(f,x,y,Z ;vmax = 1.0, cmap=Reverse(:balance),hmargs=Dict(),axargs=Dict())
    ax = Axis(f[1,1],axargs...)        
    m = findmax(abs.(Z))[1]
    range_val = (-m*vmax,m*vmax)
    hmap = heatmap!(ax,x, y, Z, colormap = cmap, colorrange=range_val, hmargs...)
    ax.aspect=DataAspect()
    Colorbar(f[1,2], colormap = cmap, limits = Float64.(range_val),tellheight=true)
    rowsize!(f.layout, 1, ax.scene.px_area[].widths[2])
    return hmap, ax
end

#curve and billiard ploting
function plot_curve!(ax, crv::AbsRealCurve; plot_normal=true, dens = 20.0)
    L = crv.length
    grid = round(Int, L*dens)
    t = range(0.0,1.0, grid)
    pts = curve(crv,t)
    lines!(ax,pts, color = :black )
    if plot_normal
        ns = normal_vec(crv,t)
        arrows!(ax,getindex.(pts,1),getindex.(pts,2), getindex.(ns,1),getindex.(ns,2), color = :black, lengthscale = 0.1)
    end
    ax.aspect=DataAspect()
end

function plot_curve!(ax, crv::AbsVirtualCurve; plot_normal=false, dens = 10.0)
    L = crv.length
    grid = round(Int, L*dens)
    t = range(0.0,1.0, grid)
    pts = curve(crv,t)
    lines!(ax,pts, color = :black, linestyle = :dash)
    if plot_normal
        ns = normal_vec(crv,t)
        arrows!(ax,getindex.(pts,1),getindex.(pts,2), getindex.(ns,1),getindex.(ns,2), color = :black, lengthscale = 0.1)
    end
    ax.aspect=DataAspect()
end

function plot_boundary!(ax, billiard::AbsBilliard; dens = 100.0, plot_normal=true)
    for curve in billiard.boundary
        plot_curve!(ax, curve; dens = dens, plot_normal = plot_normal)
    end
end

function plot_domain_fun!(f, curve::C; xlim=(-1.0,1.0),ylim=(-1.0,1.0), dens=100.0, hmargs=Dict(),cmap=:binary) where {C<:AbsCurve}
    d = one(dens)/dens
    x_grid = range(xlim... ; step=d)
    y_grid = range(ylim... ; step=d)
    Z = reshape(domain(curve,x_grid,y_grid),length(x_grid),length(y_grid))
    hmap, ax = plot_heatmap_balaced!(f,x_grid,y_grid,Z) 
    ax.aspect=DataAspect()
    return ax, hmap
end

function plot_domain!(ax, curve::AbsCurve;xlim=(-1.0,1.0),ylim=(-1.0,1.0), dens=100.0, hmargs=Dict(),cmap=Reverse(:binary))
    d = one(dens)/dens
    x_grid = range(xlim... ; step=d)
    y_grid = range(ylim... ; step=d)
    Z = reshape(is_inside(curve,x_grid,y_grid),length(x_grid),length(y_grid))
    hmap = heatmap!(ax, x_grid, y_grid, Z, colormap = cmap, colorrange=(-1,1) ,hmargs...)
    ax.aspect=DataAspect()
end
#modify for consistency
function plot_domain!(ax, billiard::AbsBilliard; dens=100.0, hmargs=Dict(),cmap=Reverse(:binary))
    d = one(dens)/dens
    #sz = (d,d)
    xlim, ylim = boundary_limits(billiard.boundary; grd=1000, type=Float32) 
    x_grid = range(xlim... ; step=d)
    y_grid = range(ylim... ; step=d)
    Z = reshape(is_inside(billiard,x_grid,y_grid),length(x_grid),length(y_grid))
    hmap = heatmap!(ax, x_grid, y_grid, Z, colormap = cmap, colorrange=(-1,1) ,hmargs...)
    ax.aspect=DataAspect()
end

function plot_lattice!(ax, billiard::AbsBilliard; dens=50.0, scargs=Dict())
    d = one(dens)/dens
    sz = (d,d)
    x_plot, y_plot, gen = interior_grid(billiard,sz)
    X = [pt.xy  for pt in gen if pt.inside]
    hmap = scatter!(ax,X)
    ax.aspect=DataAspect()
end



#wavefunction plotting

function plot_wavefunction!(f,state::AbsState, basis::AbsBasis, billiard::AbsBilliard; 
    b=5.0,dens = 10.0, inside_only=true, plot_normal=false, vmax = 1.0, cmap=Reverse(:balance),hmargs=Dict(),axargs=Dict())
    Psi, x, y = wavefunction(state,basis,billiard;b=b, inside_only=inside_only)
    #Psi[Psi .== zero(eltype(Psi))] .= NaN

    hmap, ax = plot_heatmap_balaced!(f,x,y,Psi ;vmax = vmax, cmap=cmap,hmargs=hmargs,axargs=axargs)
    plot_boundary!(ax, billiard; dens = dens, plot_normal=plot_normal)
end

function plot_wavefunction_gradient!(f,state::AbsState, basis::AbsBasis, billiard::AbsBilliard; 
    b=5.0,dens = 10.0, inside_only=true, plot_normal=false, lengthscale = 0.001, cmap=Reverse(:balance),hmargs=Dict(),axargs=Dict())
    #Psi[Psi .== zero(eltype(Psi))] .= NaN
    ax = Axis(f[1,1])  
    dX, dY, x_grid, y_grid =  wavefunction_gradient(state,basis,billiard;b=b, inside_only=inside_only)
    arrows!(ax,x_grid,y_grid, dX,dY, color = :black, lengthscale = lengthscale)
    plot_boundary!(ax, billiard; dens = dens, plot_normal=plot_normal)
    ax.aspect=DataAspect()
end

function plot_probability!(f,state::AbsState, basis::AbsBasis, billiard::AbsBilliard; 
    b=5.0,dens = 100.0,log=false, inside_only=true, plot_normal=false, vmax = 1.0, cmap=Reverse(:gist_heat),hmargs=Dict(),axargs=Dict())
    Psi, x, y = wavefunction(state,basis,billiard;b=b, inside_only=inside_only)
    Psi = abs2.(Psi)
    println("Psi type $(eltype(Psi)), $(memory_size(Psi))")
    
    hmap, ax = plot_heatmap!(f,x,y,Psi ;vmax = vmax, cmap=cmap,hmargs=hmargs,axargs=axargs,log=log)
    plot_boundary!(ax, billiard; dens = dens, plot_normal=plot_normal)
end

function plot_boundary_function!(ax,state::AbsState, basis::AbsBasis, billiard::AbsBilliard; 
    b=5.0,log=false,linesargs=Dict(),axargs=Dict())
    u, s = boundary_function(state, basis, billiard; b= b)
    if log
        lines!(ax, s, log10.(abs.(u)); linesargs...)
    else
        lines!(ax, s, u; linesargs...)
    end
end



#=
function plot_basis_function!(ax, basis::AbsBasis, i, k; xlim=(-1,1), ylim=(-1,1), grid::Tuple = (200,200))
    x_plot = LinRange(xlim... , grid[1])
    y_plot = LinRange(ylim... , grid[2])
    x = repeat(x_plot , outer = length(y_plot))
    y = repeat(y_plot , inner = length(x_plot))
    phi = basis_fun(basis, i, k, x, y) 
    Z = reshape(phi, grid)
    heatmap!(ax, x_plot,y_plot,Z, colormap = :balance)
    ax.aspect=DataAspect()
end

function plot_basis_function!(ax, basis::AbsBasis, i, k, curve::AbsCurve, sampler;  grid::Int = 200)
    t, dt = sampler(grid)
    x, y = curve.r(t)
    phi = basis_fun(basis, i, k, x, y) 
    
    lines!(ax,t,phi, color = :black )
end
=#
