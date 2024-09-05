#include("../abstracttypes.jl")


function is_equal(x::T, dx::T, y::T, dy::T) :: Bool where {T<:Real}
    # Define the intervals
    x_lower = x - dx
    x_upper = x + dx
    y_lower = y - dy
    y_upper = y + dy
    # Check if the intervals overlap
    return max(x_lower, y_lower) <= min(x_upper, y_upper)
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
    #check if intervals are empty 
    if isempty(k_left)
        #println("Left interval is empty.")
        append!(k_left, k_right)
        append!(ten_left, ten_right)
        append!(control_left, [false for i in 1:length(k_right)])
        return nothing #return short circuits further evaluation
    end

    #if right is empty just skip the mergeing
    if isempty(k_right)
        #println("Right interval is empty.")
        return nothing
    end
    
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

    fl = findlast(idx_r)
    idx_last = isnothing(fl) ? 1 : fl + 1
    append!(k_left, k_right[idx_last:end])
    append!(ten_left, ten_right[idx_last:end])
    append!(control_left, [false for i in idx_last:length(k_right)])
end

function compute_spectrum(solver::AbsSolver, basis::AbsBasis, billiard::AbsBilliard,k1,k2,dk; tol=1e-4, parallel_matrix = true)
    k0 = k1
    #initial computation
    k_res, ten_res = solve_spectrum(solver, basis, billiard, k0, dk+tol)
    control = [false for i in 1:length(k_res)]
    while k0 < k2
        k0 += dk
        k_new, ten_new = solve_spectrum(solver, basis, billiard, k0, dk+tol)
        overlap_and_merge!(k_res, ten_res, k_new, ten_new, control, k0-dk, k0; tol=tol)

    end
    return k_res, ten_res, control
end

struct SpectralData{T} 
    k::Vector{T}
    ten::Vector{T}
    control::Vector{Bool}
    k_min::T
    k_max::T
end

function SpectralData(k,ten,control)
    k_min = minimum(k)
    k_max = maximum(k)
    return SpectralData(k,ten,control,k_min,k_max)
end

function merge_spectra(s1::SpectralData, s2::SpectralData; tol=1e-4)
    # Define the overlap interval manually
    overlap_start = max(s1.k_min - tol / 2, s2.k_min - tol / 2)
    overlap_end = min(s1.k_max + tol / 2, s2.k_max + tol / 2)
    
    # Identify the indices of wavenumbers within the overlap interval
    idx_1 = [(k >= overlap_start) && (k <= overlap_end) for k in s1.k]
    idx_2 = [(k >= overlap_start) && (k <= overlap_end) for k in s2.k]

    # Extract wavenumbers and tensions in the overlap interval
    ks1 = s1.k[idx_1]
    ts1 = s1.ten[idx_1]
    ks2 = s2.k[idx_2]
    ts2 = s2.ten[idx_2]

    # Match wavenumbers and combine the results in the overlap region
    ks_ov, ts_ov, cont_ov = match_wavenumbers(ks1, ts1, ks2, ts2)
    
    # Combine the non-overlapping regions and the matched overlap region
    ks = append!(s1.k[.~idx_1], ks_ov)
    ts = append!(s1.ten[.~idx_1], ts_ov)
    control = append!(s1.control[.~idx_1], cont_ov)

    append!(ks, s2.k[.~idx_2])
    append!(ts, s2.ten[.~idx_2])
    append!(control, s2.control[.~idx_2])

    # Sort the resulting wavenumbers
    p = sortperm(ks) 
    return SpectralData(ks[p], ts[p], control[p])
end

function compute_spectrum(solver::AbsSolver,basis::AbsBasis,billiard::AbsBilliard,N1::Int,N2::Int,dN::Int; N_expect = 2.0, tol=1e-4, parallel_matrix = false)
    let solver=solver, basis=basis, billiard=billiard
        N_intervals = range(N1-dN/2,N2+dN/2,step=dN)
        #println(N_intervals)
        if hasproperty(billiard,:angles)
            k_intervals = [k_at_state(n, billiard.area, billiard.length, billiard.angles) for n in N_intervals]
        else
            k_intervals = [k_at_state(n, billiard.area, billiard.length) for n in N_intervals]
        end

        results = Vector{SpectralData}(undef,length(k_intervals)-1)
        for i in 1:(length(k_intervals)-1)
            k1 = k_intervals[i]
            k2 = k_intervals[i+1]
            dk = N_expect * 2.0*pi / (billiard.area * k1) #fix this
            #println(k1)
            #println(k2)
            #println(dk)
            k_res, ten_res, control = compute_spectrum(solver,basis,billiard,k1,k2,dk; parallel_matrix, tol)
            #println(k_res)
            results[i] = SpectralData(k_res, ten_res, control)
        end

        return reduce(merge_spectra, results)
    end
end