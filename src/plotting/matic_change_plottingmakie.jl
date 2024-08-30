using QuantumBilliards
using GLMakie

# CHANGE: Add nrows and ncols to plot the last wavefunctions or all if smaller than specified least
function plot_probability!(f::Figure,state_bundle::EigenstateBundle; 
    b=5.0,dens = 100.0, nrows::Int, ncols::Int, start_idx::Int, log=false, inside_only=true, fundamental_domain = true, plot_normal=false, 
    vmax = 1.0, cmap=Reverse(:gist_heat),hmargs=Dict(),axargs=Dict(), 
    memory_limit = 10.0e9)
    Psi_bundle, x, y = wavefunction(state_bundle;b=b,fundamental_domain=fundamental_domain, inside_only=inside_only, memory_limit=memory_limit)
    billiard = state_bundle.billiard
    
    if start_idx > length(Psi_bundle)
        start_idx = length(Psi_bundle)
    end
    end_idx = start_idx + nrows*ncols
    if end_idx > length(Psi_bundle)
        end_idx = length(Psi_bundle)
    end
    indices = start_idx:end_idx
    l_idx = length(indices)
    Psi_to_plot = Psi_bundle[indices]
    if length(indices) != nrows*ncols
        k, _ = divrem(length(indices), ncols)
        nrows = k+1
    end

    for m in 1:l_idx
        i, j = divrem(m - 1, ncols)
        i += 1
        j += 1
        P = abs2.(Psi_to_plot[m])
        hmap, ax = plot_heatmap!(f[i, j],x,y,P ;vmax = vmax, cmap=cmap,hmargs=hmargs,axargs=axargs,log=log)
        plot_boundary!(ax, billiard; dens=dens, plot_normal=plot_normal)
    end
end


# Modified plot_boundary_function! function to include rows and cols with the option to choose start_idx
function plot_boundary_function!(f::Figure, state_bundle::EigenstateBundle, nrows::Int, 
    ncols::Int, 
    start_idx::Int; b=5.0, 
    log=false, 
    linesargs=Dict(), 
    axargs=Dict(), 
)
    us, s, norms = boundary_function(state_bundle; b=b)
    billiard = state_bundle.billiard
    edges = curve_edge_lengths(billiard)
    
    if start_idx > length(us)
        start_idx = length(us)
    end
    end_idx = start_idx + nrows*ncols
    if end_idx > length(us)
        end_idx = length(us)
    end
    indices = start_idx:end_idx
    l_idx = length(indices)
    us_to_plot = us[indices]
    if length(indices) != nrows*ncols
        k, _ = divrem(length(indices), ncols)
        nrows = k+1
    end

    for m in 1:l_idx
        i, j = divrem(m - 1, ncols)
        i += 1
        j += 1
        if log
            lines!(f[i, j], s, log10.(abs.(us_to_plot[m])); linesargs...)
        else
            lines!(f[i, j], s, us_to_plot[m]; linesargs...)
        end
        # Add vertical lines for billiard edges
        vlines!(f[i, j], edges; color=:black, linewidth=0.5)
    end
end

# Modified plot_momentum_function! function to include rows and cols with the option to choose start_idx
function plot_momentum_function!(f,state_bundle::EigenstateBundle, nrows::Int, ncols::Int, start_idx::Int; 
    b=5.0, log=false, linesargs=Dict(),axargs=Dict())
    mfs, k_range = momentum_function(state_bundle; b=b)
    ks = state_bundle.ks

    if start_idx > length(mfs)
        start_idx = length(mfs)
    end
    end_idx = start_idx + nrows*ncols
    if end_idx > length(mfs)
        end_idx = length(mfs)
    end
    indices = start_idx:end_idx
    l_idx = length(indices)
    mfs_to_plot = mfs[indices]
    ks_to_plot = ks[indices]
    if length(indices) != nrows*ncols
        k, _ = divrem(length(indices), ncols)
        nrows = k+1
    end

    for m in 1:l_idx
        i, j = divrem(m - 1, ncols)
        i += 1
        j += 1
        if log
            lines!(f[i, j], k_range, log10.(abs.(mfs_to_plot[m])); linesargs...)
            vlines!(f[i, j], [ks_to_plot[m]]; color=:black, linewidth=0.5)
        else
            lines!(f[i, j], k_range, mfs[i]; linesargs...)
            vlines!(f[i, j], [ks_to_plot[m]]; color=:black, linewidth=0.5)
            xlims!(f[i, j], 0.0, 1.2*ks_to_plot[m])
        end
    end
end

# Modified plot_husimi_function! function to include rows and cols with the option to choose start_idx
function plot_husimi_function!(f,state_bundle::EigenstateBundle, nrows::Int, ncols::Int, start_idx::Int; 
    b=5.0,log=false, vmax = 1.0, cmap=Reverse(:gist_heat),hmargs=Dict(),axargs=Dict())
    billiard = state_bundle.billiard
    L = billiard.length
    us, s, norms = boundary_function(state_bundle; b=b)
    ks = state_bundle.ks
    edges = curve_edge_lengths(billiard)
    
    if start_idx > length(us)
        start_idx = length(us)
    end
    end_idx = start_idx + nrows*ncols
    if end_idx > length(us)
        end_idx = length(us)
    end
    indices = start_idx:end_idx
    l_idx = length(indices)
    us_to_plot = us[indices]
    ks_to_plot = ks[indices]
    if length(indices) != nrows*ncols
        k, _ = divrem(length(indices), ncols)
        nrows = k+1
    end

    for m in 1:l_idx
        i, j = divrem(m - 1, ncols)
        i += 1
        j += 1
        H, qs, ps = husimi_function(ks_to_plot[m],us_to_plot[m],s,L; w = 7.0)    
        hmap, ax = plot_heatmap!(f[i,j],qs,ps,H; vmax = vmax, cmap=cmap,hmargs=hmargs,axargs=axargs,log=log)
        vlines!(ax, edges; color=:black, linewidth=0.5)
    end
end



