include("../abstracttypes.jl")

using IntervalArithmetic
function is_equal(x,dx,y,dy)
    #check if numbers are equal within tolerances
    X = x ± dx
    Y = y ± dy 
    Z = X ∩ Y
    return  ~(Z == ∅)
end

function match_wavenumbers(ks_l,ts_l,ks_r,ts_r)
    #vectors ks_l and_ks_r must be sorted
    i = j = 1 #counting index
    control = Vector{Bool}()#control bits
    ks = Vector{eltype(ks_l)}()#final wavenumbers
    ts = Vector{eltype(ts_l)}()#final tensions
    while i <= length(ks_l) && j <= length(ks_r)
        x, dx = ks_l[i], ts_l[i]
        y, dy = ks_r[j], ts_r[j]
        if  is_equal(x,dx,y,dy) #check equality with errorbars
            i += 1 
            j += 1
            if dx < dy
                push!(ks, x)
                push!(ts, dx)
                push!(control, true)
            else
                push!(ks, y)
                push!(ts, dy)
                push!(control, true)
            end
        elseif x < y
            i += 1
            push!(ks, x)
            push!(ts, dx)
            push!(control, false)
        else 
            j += 1
            push!(ks, y)
            push!(ts, dy)
            push!(control, false)
        end
    end
    return ks, ts, control 
end

function overlap_and_merge!(k_left, ten_left, k_right, ten_right, control_left, kl, kr; tol=1e-3)
    #find overlaps in interval [k1,k2]
    idx_l = k_left .> (kl-tol) .&& k_left .< (kr+tol)
    idx_r = k_right .> (kl-tol) .&& k_right .< (kr+tol)
    
    ks_l,ts_l,ks_r,ts_r = k_left[idx_l], ten_left[idx_l], k_right[idx_r], ten_right[idx_r]
    #check if wavnumbers match in overlap interval
    ks, ts, control = match_wavenumbers(ks_l,ts_l,ks_r,ts_r)
    #println("left: $ks_l")
    #println("right: $ks_r")
    #println("overlaping: $ks")
    #i_l = idx_l[1]
    #i_r = idx_r[end]+1
    deleteat!(k_left, idx_l)
    append!(k_left, ks)
    deleteat!(ten_left, idx_l)
    append!(ten_left, ts)
    deleteat!(control_left, idx_l)
    append!(control_left, control)

    idx_last = findlast(idx_r) + 1
    append!(k_left, k_right[idx_last:end])
    append!(ten_left, ten_right[idx_last:end])
    append!(control_left, [false for i in idx_last:length(k_right)])

end

function compute_spectrum(solver::AbsSolver, basis::AbsBasis, pts::AbsPoints,k1,k2,dk;tol=1e-4, plot_info=false)
    k0 = k1
    #initial computation
    k_res, ten_res = solve(solver, basis, pts, k0, dk+tol)
    control = [false for i in 1:length(k_res)]
    while k0 < k2
        k0 += dk
        k_new, ten_new = solve(solver, basis, pts, k0, dk+tol)
        overlap_and_merge!(k_res, ten_res, k_new, ten_new, control, k0-dk, k0; tol=tol)

    end
    return k_res, ten_res, control
end

using Makie
function compute_spectrum_test(solver::AbsSolver, basis::AbsBasis, pts::AbsPoints,k1,k2,dk;tol=1e-4, plot_info=false)
    k0 = k1
    #initial computation
    k_res, ten_res = solve(solver, basis, pts, k0, dk+tol)
    control = [false for i in 1:length(k_res)]
    cycle = Makie.wong_colors()[1:6]
    f = Figure(resolution = (1000,1000));
    ax = Axis(f[1,1])
    scatter!(ax, k_res, log10.(ten_res),color=(cycle[1], 0.5))
    scatter!(ax, k_res, zeros(length(k_res)),color=(cycle[1], 0.5))
    #println("iteration 0")
    #println("merged: $k_res")
    
    #println("overlaping: $ks")
    i=1
    while k0 < k2
        #println("iteration $i")
        k0 += dk
        k_new, ten_new = solve(solver, basis, pts, k0, dk+tol)
        scatter!(ax, k_new, log10.(ten_new),color=(cycle[mod1(i+1,6)], 0.5))
        scatter!(ax, k_new, zeros(length(k_new)),color=(cycle[mod1(i+1,6)], 0.5))
        i+=1
        #println("new: $k_new")
        #println("overlap: $(k0-dk), $k0")
        overlap_and_merge!(k_res, ten_res, k_new, ten_new, control, k0-dk, k0; tol=tol)
        #println("merged: $k_res")
        #println("control: $control")
    end
    scatter!(ax, k_res, log10.(ten_res), color=(:black, 1.0), marker=:x,  ms = 100)
    scatter!(ax, k_res, zeros(length(k_res)), color=(:black, 1.0), marker=:x, ms = 100)
    display(f)
    return k_res, ten_res, control
end