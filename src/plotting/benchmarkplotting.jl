using Makie

function plot_benchmarks!(f, solver, basis, billiard, k, dk, d_range::AbstractArray, b)
    d_range = collect(d_range)
    b_range = [b]
    info_matrix = compute_benchmarks(solver, basis, billiard, k, dk; d_range = d_range, b_range = b_range)
    info = info_matrix[:,1]
    ks = [i.results[1] for i in info]
    ten = [i.results[2] for i in info]
    const_time = [i.construction_time for i in info]
    dec_time = [i.decomposition_time for i in info]
    ax_ten = Axis(f[1,1], xlabel=L"d", ylabel=L"\mathrm{min(SV)}")
    #println(ks)
    #println(ten)
    scatter!(ax_ten, d_range, ten)
    ax_k = Axis(f[2,1], xlabel=L"d", ylabel=L"k")
    scatter!(ax_k, d_range, ks)
    #ma = maximum(ks)
    #mi = minimum(ks)
    ylims!(ax_k, k-dk, k+dk)
    ax_tim = Axis(f[:,2], xlabel=L"d", ylabel=L"\mathrm{time}")
    scatter!(ax_tim, d_range, const_time) 
    scatter!(ax_tim, d_range, dec_time) 
end

function plot_benchmarks!(f, solver, basis, billiard, k, dk, d, b_range::AbstractArray)
    b_range = collect(b_range)
    d_range = [d]
    info_matrix = compute_benchmarks(solver, basis, billiard, k, dk; d_range = d_range, b_range = b_range)
    info = info_matrix[1,:]
    ks = [i.results[1] for i in info]
    ten = [i.results[2] for i in info]
    const_time = [i.construction_time for i in info]
    dec_time = [i.decomposition_time for i in info]
    ax_ten = Axis(f[1,1], xlabel=L"b", ylabel=L"\mathrm{min(SV)}")
    #println(ks)
    #println(ten)
    scatter!(ax_ten, b_range, ten)
    ax_k = Axis(f[2,1], xlabel=L"b", ylabel=L"k")
    scatter!(ax_k, b_range, ks)
    #ma = maximum(ks)
    #mi = minimum(ks)
    ylims!(ax_k, k-dk, k+dk)
    ax_tim = Axis(f[:,2], xlabel=L"b", ylabel=L"\mathrm{time}")
    scatter!(ax_tim, b_range, const_time) 
    scatter!(ax_tim, b_range, dec_time) 
end