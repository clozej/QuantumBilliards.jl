using StaticArrays

struct BoundaryPoints{T} <: AbsPoints where {T<:Real}
    xy::Vector{SVector{2,T}}
    normal::Vector{SVector{2,T}} #normal vectors in points
    s::Vector{T} # arc length coords
    ds::Vector{T} #integration weights
end

"""
This function creates the boundary coordinates for a given `curve`. The logic behind the program is this:

# Arguments:
- `curve`: An `AbsCurve` object - needs to implement a `curve`, `normal_vec` and `arc_length` functions 
- `sampler`: An `AbsSampler` object - for sampling the `ts`
- `N`: The number of disretization points

# Logic:
- We get the lenght `L` of the boundary from the curve struct field
- We sample points from the curve using the `sampler` given (LinearNodes, Legendre...). This gives the disretization parameters (`ts`) and the differences between them (`dts`)
- For the sample points we calculate the corresponding `xy` point on the curve using the curve equation
- We calculate the `normal` vector to the curve using the curve's derivative for each disretization point (each `t` that defines an `xy` point in the curve)
- We calculate the arc length `s` using the curve's `arc_length` field
- !The weigths (arclength differences) are gotten from the parametrization of the curve -> For nonhomogenous `dts` we get nonhomogenous `ds` arclength weights

# Returns
- `xy`: Vector of 2D points on the curve
- `normal`: Vector of 2D normal vectors to the curve
- `s`: Vector of arc length coordinates
- `ds`: Vector of integration weights
"""
function boundary_coords(crv::C,sampler::S, N) where {C<:AbsCurve, S<:AbsSampler}
    L = crv.length
    t, dt = sample_points(sampler, N)
    xy = curve(crv, t)
    normal = normal_vec(crv,t)
    s = arc_length(crv,t)
    ds = L.*dt #modify for different parametrizations
    return xy, normal, s, ds
end

"""
This function creates the boundary coordinates for a given `curve`. The logic behind the program is this:

# Arguments:
- `curve`: An `AbsCurve` object - needs to implement a `curve`, `normal_vec` and `arc_length` functions 
- `t::AbstractArray{T, 1}`: A parametrizatiomn for the `curve` object.
- `dt::AbstractArray{T, 1}`: The difference between parametrizationsThe number of disretization points

# Logic:
- We get the lenght `L` of the boundary from the curve struct field
- From the `ts` we calculate the corresponding `xy` point on the curve using the curve equation. Here we do not use the sampler but supply the points ourselves
- We calculate the `normal` vector to the curve using the curve's derivative for each disretization point (each `t` that defines an `xy` point in the curve)
- We calculate the arc length `s` using the curve's `arc_length` field
- The weigths (arclength differences) are gotten from the parametrization of the curve -> For nonhomogenous `dts` we get nonhomogenous `ds` arclength weights

# Returns
- `xy`: Vector of 2D points on the curve
- `normal`: Vector of 2D normal vectors to the curve
- `s`: Vector of arc length coordinates
- `ds`: Vector of integration weights
"""
function boundary_coords(crv::C, t, dt) where {C<:AbsCurve}
    L = crv.length
    xy = curve(crv, t)
    normal = normal_vec(crv,t)
    s = arc_length(crv,t)
    ds = L.*dt #modify for different parametrizations
    return xy, normal, s, ds
end 

#make better watch out for primes parameter
"""
Compute boundary coordinates for a billiard domain using Fourier nodes as the sampler.

# Arguments:
- `billiard`: instance of subtype of `AbsBilliard`
- `sampler`: instance of subtype of `AbsSampler`, default is `FourierNodes`
- `N`: Number of disretization points for each boundary segment

# Logic:
- Retrieve the boundary segments of the billiard domain.
- Sample each boundary segment using Fourier nodes, generating the sample points (`ts`) and intervals (`dts`).
- Compute the boundary coordinates (`xy`, `normal`, `s`, and `ds`) for the first segment.
- For subsequent segments, adjust the arc length coordinates `s` to account for the cumulative length of previous segments.
- Concatenate the results across all segments.

# Returns:
- `BoundaryPoints{T}`: A struct containing:
- `xy`: Vector of 2D points on the boundary.
- `normal`: Vector of 2D normal vectors to the boundary at each point.
- `s`: Vector of arc length coordinates for the sampled points.
- `ds`: Vector of integration weights (arc length differences).
"""
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
"""
Compute boundary coordinates for a billiard domain using a specified sampler.

# Arguments:
- `billiard`: instance of subtype of `AbsBilliard`
- `sampler`: instance of subtype of `AbsSampler`, has no default
- `N`: Number of disretization points for each boundary segment

# Logic:
- Retrieve the total length `L` of the billiard boundary.
- For each boundary segment, calculate the number of sample points `Nc` based on its length relative to the total boundary length.
- Compute the boundary coordinates (`xy`, `normal`, `s`, and `ds`) for each segment using the specified sampler.
- Adjust the arc length coordinates `s` for each segment to account for the cumulative length of previous segments.
- Concatenate the results across all segments.

# Returns:
- `BoundaryPoints{T}`: A struct containing:
- `xy`: Vector of 2D points on the boundary.
- `normal`: Vector of 2D normal vectors to the boundary at each point.
- `s`: Vector of arc length coordinates for the sampled points.
- `ds`: Vector of integration weights (arc length differences).
"""
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

"""
Compute dilated boundary points by moving each point on the billiard boundary outward by a wavelength.

# Description:
This function calculates dilated boundary points for a billiard domain by shifting each boundary point outward along its normal vector by a distance equal to the wavelength corresponding to the wavenumber `k`.

# Arguments:
- `billiard`: An instance of subtype of `AbsBilliard`.
- `sampler`: An instance of subtype of `AbsSampler`, default is unavailable.
- `N`: The number of sample points to generate for each boundary segment.
- `k`: The wavelenght

# Logic:
- Calculate the wavelength `λ` as `λ = 2π/k`.
- Compute the boundary coordinates (`xy`, `normal`) using the specified sampler and number of sample points.
- Dilate the boundary by adding the product of the normal vectors and the wavelength `λ` to the original boundary points.

# Returns:
- A vector of 2D points representing the dilated boundary.
"""
function dilated_boundary_points(billiard::Bi, sampler::S, N, k) where {Bi<:AbsBilliard, S<:AbsSampler}
    lam = 2*pi/k #wavelength
    #M = max(N,100)
    pts = boundary_coords(billiard, sampler, N)
    xy = pts.xy 
    n = pts.normal
    return xy .+ lam .* n
end