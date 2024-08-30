using GLMakie
using LinearAlgebra

# Test to see how to implement only plotting part of the results
if !isdir("PlottingTest")
    mkdir("PlottingTest")
end

# Function to plot the probability and save the figure
function plot_test!(
    fig::Figure, 
    Psi_bundle::Vector{Matrix{Float64}}; 
    nrows::Int, 
    ncols::Int, 
    log=false, 
    vmax=1.0, 
    cmap=:viridis, 
    hmargs=Dict(), 
    axargs=Dict(), 
    start_idx::Int = length(Psi_bundle)  # Default value indicating the last nrows * ncols elements
)
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
        P = abs2.(Psi_to_plot[m])  # Calculate the probability density
        ax = Axis(fig[i, j], title="Heatmap ($i,$j)")
        ax.aspect = DataAspect()
        heatmap!(ax, P; colormap=cmap, colorrange=(0, vmax))
    end
    println("Start index: $start_idx")
    println("End index: $end_idx")
    println("Actual elements to plot: $l_idx")
end

# Testing function
function test_plot()
    # Define the number of rows and columns
    nrows = 2
    ncols = 3

    n = 7
    matrix_size = 5
    state_bundle = [zeros(matrix_size, matrix_size) for _ in 1:n]
    for i in 1:n
        row = div(i - 1, matrix_size) + 1 
        col = mod(i - 1, matrix_size) + 1
        state_bundle[i][col, end-row+1] = 1
    end

    # Display the result
    for matrix in state_bundle
        println(matrix)
    end

    # start_idx = length(state_bundle)
    println("default")
    fig_default = Figure(resolution = (800, 600))
    plot_test!(fig_default, state_bundle, nrows=nrows, ncols=ncols)
    save("PlottingTest/plot_default.png", fig_default)

    # start_idx in middle of length(vec)
    println("middle")
    fig_middle = Figure(resolution = (800, 600))
    plot_test!(fig_middle, state_bundle, nrows=nrows, ncols=ncols, start_idx=3)
    save("PlottingTest/plot_middle.png", fig_middle)
    
    # start_idx=1
    println("start")
    fig_start = Figure(resolution = (800, 600))
    plot_test!(fig_start, state_bundle, nrows=nrows, ncols=ncols, start_idx=1)
    save("PlottingTest/plot_start.png", fig_start)
    
    # start_idx>lenght(vec)
    println("almost ending")
    fig_large_start = Figure(resolution = (800, 600))
    plot_test!(fig_large_start, state_bundle, nrows=nrows, ncols=ncols, start_idx=10)
    save("PlottingTest/plot_ending_almost.png", fig_large_start)

    #  Start index = length(vec) - 1
    println("second last")
    fig_second_end = Figure(resolution = (800, 600))
    plot_test!(fig_second_end, state_bundle, nrows=nrows, ncols=ncols, start_idx=length(state_bundle) - 1)
    save("PlottingTest/plot_second_last.png", fig_second_end)
end

test_plot()