using StaticArrays

struct BoundaryPoints{T} <: AbsPoints where {T<:Real}
    xy::Vector{SVector{2,T}}
    normal::Vector{SVector{2,T}} #normal vectors in points
    s::Vector{T} # arc length coords
    ds::Vector{T} #integration weights
end

function boundary_coords(crv::C,sampler::S, N) where {C<:AbsCurve, S<:AbsSampler}
    L = crv.length
    t, dt = sample_points(sampler, N)
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
function boundary_coords(billiard::Bi, sampler::FourierNodes, N) where {Bi<:AbsBilliard}
    let boundary = billiard.full_boundary
        #crv_lengths = [crv.length for crv in boundary]
        ts, dts = sample_points(sampler, N)

        xy_all, normal_all, s_all, ds_all = boundary_coords(boundary[1], ts[1], dts[1])
        l = boundary[1].length
        for i in 2:length(ts)
            crv = boundary[i]
            if (typeof(crv) <: AbsRealCurve)
                Lc = crv.length
                xy,nxy,s,ds = boundary_coords(crv, ts[i], dts[i])
                append!(xy_all, xy)
                append!(normal_all, nxy)
                s = s .+ l
                append!(s_all, s)
                append!(ds_all, ds)
                l += Lc
            end    
        end
        return BoundaryPoints(xy_all,normal_all,s_all,ds_all) 
    end
end


#make better watch out for primes parameter
function boundary_coords(billiard::Bi, sampler::S, N) where {Bi<:AbsBilliard, S<:AbsSampler}
    let boundary = billiard.full_boundary
            L = billiard.length
            Lc = boundary[1].length
            Nc = round(Int, N*Lc/L)
            xy_all, normal_all, s_all, ds_all = boundary_coords(boundary[1],sampler,Nc)
            #println(s_all)
            l = boundary[1].length #cumulative length
            for crv in boundary[2:end]
                if (typeof(crv) <: AbsRealCurve )
                    Lc = crv.length
                    Nc = round(Int, N*Lc/L)
                    xy,nxy,s,ds = boundary_coords(crv,sampler,Nc)
                    append!(xy_all, xy)
                    append!(normal_all, nxy)
                    s = s .+ l
                    append!(s_all, s)
                    append!(ds_all, ds)
                    l += Lc
                end    
            end
        return BoundaryPoints(xy_all,normal_all,s_all,ds_all) 
    end
end

function dilated_boundary_points(billiard::Bi, sampler::S, N, k) where {Bi<:AbsBilliard, S<:AbsSampler}
    lam = 2*pi/k #wavelength
    #M = max(N,100)
    pts = boundary_coords(billiard, sampler, N)
    xy = pts.xy 
    n = pts.normal
    return xy .+ lam .* n
end