include("../abstracttypes.jl")
using Makie
#not exported
function plot_matrix!(f, A; log=false)
    Z = deepcopy(A)
    ax = Axis(f[1,1])
    if log
        m1 = log10(findmax(Z)[1])
        m2 = log10(abs(findmin(Z)[1]))
        S = sign.(Z)
        Z = abs.(Z)
        Z[Z.<eps()] .= NaN
        Z = log10.(Z)
        N = deepcopy(Z)
        N[S .>= 0.0] .= NaN #mask positive amplitudes
        Z[S .< 0.0] .= NaN #mask negative amplitudes

        println("$m1, $m2")
        range_val1 = (-15,m1)
        range_val2 = (-15,m2)

        hmap1 = heatmap!(ax, N, colormap = :algae, colorrange=range_val1)
        hmap2 = heatmap!(ax, Z, colormap = :amp, colorrange=range_val2)
        ax.yreversed=true
        ax.aspect=DataAspect()
        Colorbar(f[1,2], colormap = :algae, limits = Float64.(range_val1),tellheight=true)
        rowsize!(f.layout, 1, ax.scene.px_area[].widths[2])
        Colorbar(f[1,3], colormap = :amp, limits = Float64.(range_val2),tellheight=true)
        rowsize!(f.layout, 1, ax.scene.px_area[].widths[2])
        #Colorbar(f[1, 2], colormap = :balance, limits = Float64.(range_val))
    else
        m = findmax(abs.(Z))[1]
        Z[abs.(Z).<eps()] .= NaN
        println("$m")
        range_val = (-m,m) 
        hmap = heatmap!(ax, Z, colormap = :balance, colorrange=range_val)
        ax.yreversed=true
        ax.aspect=DataAspect()
        Colorbar(f[1,2], colormap = :balance, limits = Float64.(range_val),tellheight=true)
        rowsize!(f.layout, 1, ax.scene.px_area[].widths[2])
    end
    return f
end