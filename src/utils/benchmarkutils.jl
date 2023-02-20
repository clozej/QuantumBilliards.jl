include("../abstracttypes.jl")
#include("../solvers/particularsolutionsmethod.jl")

memory_size(a) = Base.format_bytes(Base.summarysize(a)) 

#=
function benchmark_solver(solver::ParticularSolutionsMethod, basis::AbsBasis, billiard::AbsBilliard, sampler::Function, k, dk; btimes = 5, print_info=true, kwargs...)
    L = real_length(billiard)
    dim = round(Int, L*k*solver.dim_scaling_factor/(2*pi))
    basis_new = resize_basis(basis, dim)
    
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
    
    mat_res, mat_time = run(construct_matrices, solver, basis_new, pts, k)
    
    B, B_int = mat_res
    svd_res, decomp_time = run(svdvals, B, B_int)
    k0, t0 = solve_wavenumber(solver, basis, pts, k, dk)

    if print_info
        println("Basis dimension: $dim")
        println("Boundary point number: $(length(pts.x))")
        println("Interior point number: $(length(pts.x_int))")
        println("Matrix construction time: $mat_time s")
        println("B matrix size: $(size(B))")
        println("B_int size: $(size(B_int))")
        println("Decomposition time: $decomp_time s")
        println("Solver time: $(decomp_time + mat_time) s")
        println("Wavenumber: k = $k0")
        println("Minimum: t = $t0")
    end
    
    return k0, t0, decomp_time, mat_time
end

function benchmark_solver(solver::DistributedParticularSolutionsMethod, basis_collection, billiard::AbsBilliard, sampler::Function, k, dk; btimes = 5, print_info=true, kwargs...)
    basis_vec = similar(basis_collection)
    for i in 1:length(basis_vec)
        L = billiard.boundary[i].length
        dim = round(Int, L*k*solver.dim_scaling_factor[i]/(2*pi))
        basis_vec[i] = resize_basis(basis_collection[i], dim)
    end
    
    pts_vec = evaluate_points(solver, billiard, sampler, k)
       
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
    
    mat_res, mat_time = run(construct_matrices, solver, basis_vec, pts_vec, k)
    
    B, B_int = mat_res
    svd_res, decomp_time = run(svdvals, B, B_int)
    k0, t0 = solve_wavenumber(solver, basis_vec, pts_vec, k, dk)

    if print_info
        println("Basis dimension: $([basis_vec[i].dim for i in eachindex(basis_vec)])")
        println("Boundary point number: $([length(pts.x) for pts in pts_vec])")
        println("Interior point number: $([length(pts.x_int) for pts in pts_vec])")
        println("Matrix construction time: $mat_time s")
        println("B matrix size: $(size(B))")
        println("B_int size: $(size(B_int))")
        println("Decomposition time: $decomp_time s")
        println("Solver time: $(decomp_time + mat_time) s")
        println("Wavenumber: k = $k0")
        println("Minimum: t = $t0")
    end
    
    return k0, t0, decomp_time, mat_time
end

function benchmark_solver(solver::ScalingMethod, basis::AbsBasis, billiard::AbsBilliard, sampler::Function, k, dk; 
    btimes = 1, print_info=true, plot_info=false, fig_res=(1000,2000), log_mat=false, limits=[(-2.0,2.0),(-2.0,2.0)], kwargs...)
    L = real_length(billiard)
    dim = round(Int, L*k*solver.dim_scaling_factor/(2*pi))
    basis_new = resize_basis(basis, dim)
    
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
    println("Constructing matrix.")
    mat_res, mat_time = run(construct_matrices, solver, basis_new, pts, k)
    
    F, Fk   = mat_res
    println("Solving decomposition.")
    sm_res, decomp_time = run(solve, solver,F,Fk, k, dk)
    ks, ten = sm_res

    if print_info
        println("Basis dimension: $dim")
        println("Boundary point number: $(length(pts.x))")
        println("Matrix construction time: $mat_time s")
        println("F matrix size: $(size(F))")
        println("Fk size: $(size(Fk))")
        println("Decomposition time: $decomp_time s")
        println("Solver time: $(decomp_time + mat_time) s")
        #println("Wavenumber: k = $k0")
        println("Minimum: t = $(minimum(ten))")
    end
    if plot_info
        f = Figure(resolution = fig_res);
        axis = Axis(f[1,1:2])
        xlim, ylim = limits
        plot_basis_function!(axis,basis_new, 1, k; xlim=xlim, ylim=ylim)
        plot_billiard!(axis, billiard)
        plot_matrix!(f[2,1], F; log=log_mat)
        plot_matrix!(f[2,2], Fk; log=log_mat)
        ax = Axis(f[3,1:2])
        scatter!(ax, ks, log10.(ten), color = :black)
        display(f)
    end
    return ks, ten, decomp_time, mat_time
end
=#