using StaticArrays

struct BoundaryPoints{T} <: AbsPoints where {T<:Real}
    xy::Vector{SVector{2,T}}
    normal::Vector{SVector{2,T}} #normal vectors in points
    s::Vector{T} # arc length coords
    ds::Vector{T} #integration weights
end

function boundary_coords(crv::C, N; sampler=fourier_nodes) where {C<:AbsCurve}
    L = crv.length
    t, dt = sampler(N)
    xy = curve(crv, t)
    normal = normal_vec(crv,t)
    s = arc_length(crv,t)
    ds = L.*dt #modify for different parametrizations
    return xy, normal, s, ds
end

function boundary_coords(crv::C, t, dt) where {C<:AbsCurve}
    L = crv.length
    xy = curve(crv, t)
    normal = normal_vec(crv,t)
    s = arc_length(crv,t)
    ds = L.*dt #modify for different parametrizations
    return xy, normal, s, ds
end 

#make better watch out for primes parameter
function boundary_coords(billiard::Bi, N; sampler=fourier_nodes, primes=true) where {Bi<:AbsBilliard}
    let boundary = billiard.full_boundary
        if sampler==fourier_nodes
            crv_lengths = [crv.length for crv in boundary]
            if primes
                ts, dts = fourier_nodes(N, crv_lengths)
            else
                ts, dts = fourier_nodes(N, crv_lengths;primes=false)
            end

            xy_all, normal_all, s_all, ds_all = boundary_coords(boundary[1], ts[1], dts[1])
            l = boundary[1].length
            for i in 2:length(ts)
                crv = boundary[i]
                if (typeof(crv) <: AbsRealCurve)
                    Lc = crv.length
                    #Nc = round(Int, N*Lc/L)
                    xy,nxy,s,ds = boundary_coords(crv, ts[i], dts[i])
                    append!(xy_all, xy)
                    append!(normal_all, nxy)
                    s = s .+ l
                    append!(s_all, s)
                    append!(ds_all, ds)
                    l += Lc
                end    
            end
        else
            L = billiard.length

            Lc = boundary[1].length
            Nc = round(Int, N*Lc/L)
            xy_all, normal_all, s_all, ds_all = boundary_coords(boundary[1], Nc; sampler=sampler)
            #println(s_all)
            l = boundary[1].length #cumulative length
            for crv in boundary[2:end]
                if (typeof(crv) <: AbsRealCurve )
                    Lc = crv.length
                    Nc = round(Int, N*Lc/L)
                    xy,nxy,s,ds = boundary_coords(crv, Nc; sampler=sampler)
                    append!(xy_all, xy)
                    append!(normal_all, nxy)
                    s = s .+ l
                    append!(s_all, s)
                    append!(ds_all, ds)
                    l += Lc
                end    
            end
        end
        return BoundaryPoints(xy_all,normal_all,s_all,ds_all) 
    end
end

function dilated_boundary_points(billiard::Bi, N, k; sampler=fourier_nodes, primes=false) where {Bi<:AbsBilliard}
    lam = 2*pi/k #wavelength
    #M = max(N,100)
    pts = boundary_coords(billiard, N; sampler=sampler, primes = primes)
    xy = pts.xy 
    n = pts.normal
    return xy .+ lam .* n
end