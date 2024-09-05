using ForwardDiff
using StaticArrays

#TODO Needs testing
"""
Calculates the corner corrections to Weyl's law with the boundary ones

#Arguments
- The corners of the boundary

Returns:
- The corner corrections
"""
function corner_correction(corner_angles::Vector{T}) where T<:Real
    return sum([(pi^2 - c^2)/(24*pi*c) for c in corner_angles])
end    

"""
The first derivative of the curve (x(t),y(t)) for each t in ts

#Arguments
- Each curve that constructs the boundary
- The parametrizations that are constructed from the boundary and sampler

Returns:
A vector of SVector's containing the first derivatives computed at the t in ts
"""
function first_derivative(crv::QuantumBilliards.AbsRealCurve, ts::Vector{T}) where T<:Real
    if hasmethod(tangent, Tuple{typeof(crv), AbstractArray{T,1}}) # Check if the curve method exists for the AbsRealCurve subtype
        return tangent(crv, ts)
    else
        return nothing
    end
end

"""
The second derivative of the curve (x(t),y(t)) for each t in ts

#Arguments
- Each curve that constructs the boundary
- The parametrizations that are constructed from the boundary and sampler

Returns:
A vector of SVector's containing the second derivatives computed at the t in ts
"""
function second_derivative(crv::QuantumBilliards.AbsRealCurve, ts::Vector{T}) where T<:Real
    if hasmethod(tangent, Tuple{typeof(crv), AbstractArray{T,1}}) # Check if the curve method exists for the AbsRealCurve subtype
        # Initialize an array to store the results
        result = Vector{SVector{2, T}}(undef, length(ts))
        h = T(1e-8)
        # Iterate over each t in ts
        for i in 1:length(ts)
            t = ts[i]
            # Clamp the t+h and t-h within [0, 1]
            t_forward = min(t + h, T(1))
            t_backward = max(t - h, T(0))
            # Evaluate the curve at t+h and t-h
            pt_forward = curve(crv, [t_forward])[1]
            pt_backward = curve(crv, [t_backward])[1]
            # Compute the numerical second derivative using finite differences
            derivative_at_t = (pt_forward - pt_backward) / (2 * h)
            # Ensure that the derivative result is stored as SVector{2, Float64}
            result[i] = derivative_at_t
        end
        
        return result
    else
        return nothing
    end  
end

"""
Computes the curvature of the boundary

# Arguments
- Each curve that constructs the boundary
- The parametrizations that are constructed from the boundary and sampler

Returns:
A vector of real numbers that are the curvatures computed at the t in ts
"""
function curvature(crv::M, ts::Vector{T}) where {T<:Real, M<:QuantumBilliards.AbsRealCurve}
    if crv isa QuantumBilliards.LineSegment
        println("line segment")
        return [0 for _ in ts]  # Curvature is 0 for a straight line
    elseif crv isa QuantumBilliards.CircleSegment || crv isa QuantumBilliards.DispersingCircleSegment
        println("circle segment")
        return [1/crv.radius for _ in ts]  # Constant curvature for circular segments
    else
        println("other curve")
        # Compute the curvature using numerical differentiation
        first_derivs = first_derivative(crv, ts)
        second_derivs = second_derivative(crv, ts)
        curvature_result = Vector{T}(undef, length(ts))
        for i in 1:length(ts)
            dx_dt, dy_dt = first_derivs[i]
            d2x_dt2, d2y_dt2 = second_derivs[i]
            num = abs(dx_dt * d2y_dt2 - dy_dt * d2x_dt2)
            denom = (dx_dt^2 + dy_dt^2)^(3/2)
            curvature_result[i] = num / denom
        end
        return curvature_result
    end

end

"""
Computes the total curvature of the boundary

# Arguments
- The scaling method to use
- The billiard/boundary to compute the curvarure on
- The wavenumber k for curvature computation as the discretizations is dependant on it


Returns:
- The total curvature of the boundary
"""
function curvature_correction(solver::QuantumBilliards.AbsScalingMethod, billiard::QuantumBilliards.AbsBilliard, k::T) where T<:Real
    bs, samplers = QuantumBilliards.adjust_scaling_and_samplers(solver, billiard)
    curves = billiard.fundamental_boundary
    curvature_corrections = 0.0
    for i in eachindex(curves)
        println("Num of curves: ", i)
        crv = curves[i]
        L = crv.length
        N = max(solver.min_pts,round(Int, k*L*bs[i]/(2*pi)))
        sampler = samplers[i]
        ts, dts = QuantumBilliards.sample_points(sampler, N)
        curvatures = curvature(crv, ts)
        curvature_boundary_integral = sum([curvatures[i] * 2*pi*dts[i] for i in eachindex(curvatures)])
        curvature_corrections += curvature_boundary_integral
    end
    return 1/(12*pi)*curvature_corrections
end

"""
The total correction to Weyl's law, incorporating both the curvature corrections and the corner corrections

# Arguments
- The scaling method to use
- The billiard/boundary to compute the curvature on
- The wavenumber k for curvature computation as the discretizations is dependant on
- The corners of the boundary as a vector of values in radians

Returns:
- The total correction to Weyl's law
"""
function C(solver::QuantumBilliards.AbsScalingMethod, billiard::QuantumBilliards.AbsBilliard, k::T, corner_angles::Vector{T}) where T<:Real
    return curvature_correction(solver, billiard, k) + corner_correction(corner_angles)
end

"""
The unfolding function in billiards

# Arguments
- The billiard/boundary to compute the curvature on
- The wavenumber k for curvature computation as the discretizations is dependant on
- The scaling method to use
- The corners of the boundary as a vector of values in radians

# Returns:
- The unfolding function
"""
function weyl_law(A::T,L::T,k::T, solver::QuantumBilliards.AbsScalingMethod, billiard::QuantumBilliards.AbsBilliard, corner_angles::Vector{T}) where T<:Real
    A * k^2 - L * k/(4*pi) - C(solver, billiard, k, corner_angles)
end

# TESTING CIRCLE

full_circle(cap_radius) = [CircleSegment(cap_radius, 0.0, Float64(2*pi), 0.0, 0.0; origin=SVector(0.0, 0.0), rot_angle=0.0)]
quarter_circle(cap_radius) = [CircleSegment(cap_radius, 0.5pi, 0.5pi, 0.0, 0.0; origin=SVector(0.0, 0.0), rot_angle=0.0)]

struct FullCircle{T} <: QuantumBilliards.AbsBilliard where {T<:Real}
    fundamental_boundary::Vector{CircleSegment{T}}
    full_boundary::Vector{CircleSegment{T}}
    length::T
    area::T
    cap_radius::T
end

function FullCircle(radius::T) where T<:Real
    fundamental_boundary = quarter_circle(radius)
    full_boundary = full_circle(radius)
    area = pi*radius^2
    length = 2*pi*radius
    cap_radius = radius
    return FullCircle{T}(fundamental_boundary, full_boundary, length, area, cap_radius)
end

circle = FullCircle(1.0)

d = 1.5
b = [20.0] # only 1 curve, not composite
solver = QuantumBilliards.ScalingMethodA(d,b)

println(corner_correction(Float64[])) # Check if we add an empty vector it does not crash
println(C(solver, circle, 10.0, Float64[])) # Should be 1/6 ≈ 0.16666...
println(curvature(quarter_circle(2.0)[1], [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])) # Should be all 0.5 -> 1/(r=2)
