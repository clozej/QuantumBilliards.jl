#include("../abstracttypes.jl")
#include("../solvers/particularsolutionsmethod.jl")
#include("../plotting/matrixplotting.jl")
using Makie

"""
Get the memory size of object a` in a human-readable format using `Base.format_bytes` to convert the size into bytes, kilobytes, megabytes, etc.
"""
memory_size(a) = Base.format_bytes(Base.summarysize(a)) 

"""
A struct to store benchmarking information for a solver.

# Description
`BenchmarkInfo` stores the results and details of benchmarking a solver. This includes the solver used, matrix dimensions and memory usage, times for matrix construction, decomposition, and solution, as well as the final results.

# Fields
- `solver::AbsSolver`: The solver used in the benchmarking.
- `matrix_dimensions::Vector`: A vector containing the dimensions of the matrices involved.
- `matrix_memory::Vector`: A vector containing the memory usage of the matrices.
- `construction_time::Float64`: The time taken to construct the matrices.
- `decomposition_time::Float64`: The time taken to decompose the matrices.
- `solution_time::Float64`: The time taken to solve the problem using the matrices.
- `results::Tuple`: The results of the solution, typically the wavenumber and tension or something else.
"""
struct BenchmarkInfo
    solver::AbsSolver
    matrix_dimensions::Vector
    matrix_memory::Vector
    construction_time::Float64
    decomposition_time::Float64
    solution_time::Float64
    results::Tuple
end

"""
Print the benchmarking information stored in a `BenchmarkInfo` object.

# Description
This function prints detailed information about the benchmarking process, including the solver used, matrix dimensions, memory usage, construction time, decomposition time, solution time, and final results.

# Arguments
- `info::BenchmarkInfo`: The `BenchmarkInfo` object containing the benchmarking details to be printed.
"""
function print_benchmark_info(info::BenchmarkInfo)
    println("Solver: $(info.solver)")
    for i in eachindex(info.matrix_dimensions)
        println("Matrix $i, $(info.matrix_dimensions[i]), memory: $(info.matrix_memory[i])")
    end
    println("Construction time: $(info.construction_time) s")
    println("Decomposition time: $(info.decomposition_time) s")
    println("Solution time: $(info.solution_time) s")
    println("Results: $(info.results)")
end

"""
Benchmark a solver with a given basis and billiard domain.

# Description
This function benchmarks a solver by constructing matrices, decomposing them, and solving the problem for given wavenumbers. It measures the time taken for each of these steps and optionally prints the results and plots the matrices.

# Logic
- The solver is used to resize the basis and evaluate points within the billiard domain.
- The matrices are constructed and benchmarked for time.
- Depending on the solver type (AcceleratedSolver vs. SweepSolver), the matrices are either solved directly or in a sweep, with times recorded for each step.
- If `print_info` is true, the benchmark information is printed via `print_benchmark_info`.
- If `plot_matrix` is true, the matrices are plotted using Makie.
- The function returns a `BenchmarkInfo` object containing all the benchmarking details.

# Arguments
- `solver::AbsSolver`: The solver to be benchmarked.
- `basis::AbsBasis`: The basis functions used in the benchmarking.
- `billiard::AbsBilliard`: The billiard domain in which the problem is solved.
- `k`: The reference wavenumber for which the problem is solved.
- `dk`: The range of wavenumbers around `k`.
- `btimes`: The number of times to repeat the benchmark for averaging. Default is 1.
- `print_info`: Whether to print the benchmark information. Default is true.
- `plot_matrix`: Whether to plot the matrices. Default is false.
- `log`: Whether to plot the matrices on a logarithmic scale. Default is false.

# Returns
- A `BenchmarkInfo` object containing the benchmarking details.
"""
function benchmark_solver(solver::AbsSolver, basis::AbsBasis, billiard::AbsBilliard, k, dk; btimes = 1, print_info=true, plot_matrix=false,log=false, kwargs...) 
    let L = billiard.length, dim = round(Int, L*k*solver.dim_scaling_factor/(2*pi))
        basis_new = resize_basis(basis,billiard, dim, k)
        pts = evaluate_points(solver, billiard, k)

        function run(fun, args...; kwargs...)
            times = Float64[]
            values = []
            for i in 1:btimes
                res = @timed fun(args...; kwargs...)
                push!(times,res.time)
                push!(values, res.value)
            end
            return values[1], mean(times)
        end
        if print_info
            mat_res, mat_time = run(construct_matrices_benchmark, solver, basis_new, pts, k)
        #returns also inforamtion about the basis matrix in the last element
        else
            mat_res, mat_time = run(construct_matrices, solver, basis_new, pts, k)
    
        end
        mat_dim = []
        mat_mem = []
        for A in mat_res
            push!(mat_dim, size(A))
            push!(mat_mem, memory_size(A))
        end
        
        if typeof(solver) <:AcceleratedSolver
            sol, decomp_time = run(solve, solver, mat_res[1],mat_res[2],k,dk)
            ks, ten = sol
            idx = findmin(abs.(ks.-k))[2]
            res = ks[idx], ten[idx]
            info = BenchmarkInfo(solver,mat_dim,mat_mem,mat_time,decomp_time,mat_time+decomp_time,res)
        end
        
        if typeof(solver) <:SweepSolver
            t, decomp_time = run(solve, solver, mat_res[1],mat_res[2])
            res, sol_time = run(solve_wavenumber,solver, basis, billiard,k,dk)
            info = BenchmarkInfo(solver,mat_dim,mat_mem,mat_time,decomp_time,sol_time,res) 
        end
                
        if print_info
            print_benchmark_info(info)
        end
        
        if plot_matrix
            mat_n = length(mat_res)
            f = Figure(resolution = (500*mat_n,500))
            for i in 1:mat_n
                plot_matrix!(f[1,i],mat_res[i],log=log)
            end
            display(f)
        end

        return info

    end
end

"""
Benchmark a solver over a range of scaling factors.

# Description
This function benchmarks a solver across different combinations of dimensional (`d_range`) and point (`b_range`) scaling factors. It returns a matrix of `BenchmarkInfo` objects, each containing the results for a specific combination.

# Arguments
- `solver`: The solver to benchmark.
- `basis`: The basis functions used.
- `billiard`: The billiard domain.
- `k`: The wavenumber to solve for.
- `dk`: The wavenumber range.
- `d_range`: A range of dimensional scaling factors (default `[2.0]`).
- `b_range`: A range of point scaling factors (default `[2.0]`).
- `btimes`: Number of benchmark repetitions for averaging (default `1`).

# Returns
- A matrix of `BenchmarkInfo` objects for each combination of scaling factors.
"""
function compute_benchmarks(solver, basis, billiard, k, dk; d_range = [2.0], b_range=[2.0],btimes=1)
    #eps = solver.eps
    make_solver(solver,d,b) = typeof(solver)(d,b) 
    grid_indices = CartesianIndices((length(d_range), length(b_range)))
    info_matrix = [benchmark_solver(make_solver(solver,d_range[i[1]],b_range[i[2]]),basis,billiard,k,dk;print_info=false,btimes=btimes) for i in grid_indices]
    return info_matrix
end
