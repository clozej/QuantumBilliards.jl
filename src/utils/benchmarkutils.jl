#include("../abstracttypes.jl")
#include("../solvers/particularsolutionsmethod.jl")
#include("../plotting/matrixplotting.jl")
using Makie


memory_size(a) = Base.format_bytes(Base.summarysize(a)) 

struct BenchmarkInfo
    solver::AbsSolver
    matrix_dimensions::Vector
    matrix_memory::Vector
    construction_time::Float64
    decomposition_time::Float64
    solution_time::Float64
    results::Tuple
end

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

function benchmark_solver(solver::AbsSolver, basis::AbsBasis, billiard::AbsBilliard, sampler::Function, k, dk; btimes = 1, print_info=true, plot_matrix=false,log=false, kwargs...) 
    let L = real_length(billiard), dim = round(Int, L*k*solver.dim_scaling_factor/(2*pi))
        basis_new = resize_basis(basis,billiard, dim, k)
        pts = evaluate_points(solver, billiard, sampler, k)

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
            res, sol_time = run(solve_wavenumber,solver, basis, billiard,k,dk;sampler=sampler)
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


function compute_benchmarks(solver, basis, billiard, sampler, k, dk; d_range = [2.0], b_range=[2.0],btimes=1)
    eps = solver.eps
    make_solver(solver,d,b) = typeof(solver)(d,b,eps) 
    grid_indices = CartesianIndices((length(d_range), length(b_range)))
    info_matrix = [benchmark_solver(make_solver(solver,d_range[i[1]],b_range[i[2]]),basis,billiard,sampler,k,dk;print_info=false,btimes=btimes) for i in grid_indices]
    return info_matrix
end
