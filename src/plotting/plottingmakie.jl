include("../abstracttypes.jl")
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
function plot_curve!(ax, curve::AbsRealCurve; plot_normal=true, dens = 20.0)
    L = curve.length
    grid = round(Int, L*dens)
    t = LinRange(0,1, grid)
    x,y = curve.r(t)
    lines!(ax,x,y, color = :black )
    if plot_normal
        nx, ny = curve.n(t)
        arrows!(ax,x,y,nx,ny, color = :black, lengthscale = 0.1)
    end
    ax.aspect=DataAspect()
end

function plot_curve!(ax, curve::AbsVirtualCurve; plot_normal=false, dens = 10.0)
    L = curve.length
    grid = round(Int, L*dens)
    t = LinRange(0,1, grid)
    x,y = curve.r(t)
    
    lines!(ax,x,y, color = :black, linestyle = :dash)
    
    if plot_normal
        nx, ny = curve.n(t)
        arrows!(ax,x,y,nx,ny, color = :black, lengthscale = 0.1)
    end
    ax.aspect=DataAspect()
end

function plot_billiard!(ax, billiard::AbsBilliard; dens = 10.0, plot_normal=true)
    for curve in billiard.boundary
        plot_curve!(ax, curve; dens = dens, plot_normal = plot_normal)
    end
end

function plot_wavefunction!(f,state::AbsState, basis::AbsBasis, billiard::AbsBilliard; 
    b=5.0,dens = 10.0, inside_only=true, plot_normal=false, vmax = 1.0, cmap=Reverse(:balance),hmargs=Dict(),axargs=Dict())
    Psi, x, y = wavefunction(state,basis,billiard;b=b, inside_only=inside_only)
    hmap, ax = plot_heatmap_balaced!(f,x,y,Psi ;vmax = vmax, cmap=cmap,hmargs=hmargs,axargs=axargs)
    plot_billiard!(ax, billiard; dens = dens, plot_normal=plot_normal)
end

function plot_probability!(f,state::AbsState, basis::AbsBasis, billiard::AbsBilliard; 
    b=5.0,dens = 100.0,log=false, inside_only=true, plot_normal=false, vmax = 1.0, cmap=Reverse(:gist_heat),hmargs=Dict(),axargs=Dict())
    Psi, x, y = wavefunction(state,basis,billiard;b=b, inside_only=inside_only)
    P = abs2.(Psi)
    hmap, ax = plot_heatmap!(f,x,y,P ;vmax = vmax, cmap=cmap,hmargs=hmargs,axargs=axargs,log=log)
    plot_billiard!(ax, billiard; dens = dens, plot_normal=plot_normal)
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
