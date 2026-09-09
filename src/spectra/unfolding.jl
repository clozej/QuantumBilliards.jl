

corner_correction(corner_angles) =  sum([(pi^2 - c^2)/(24*pi*c) for c in corner_angles])

weyl_law(k,A,L) =  @. (A * k^2 - L * k)/(4*pi)
weyl_law(k,A,L,corner_angles) =  weyl_law(k,A,L) .+ corner_correction(corner_angles)


function k_at_state(state, A, L)
    a = A
    b = -L
    c = -state*4*pi
    dis = sqrt(b^2-4*a*c)
    return (-b+dis)/(2*a)
end

function k_at_state(state, A, L, corner_angles)
    a = A
    b = -L
    c = (corner_correction(corner_angles)-state)*4*pi 
    dis = sqrt(b^2-4*a*c)
    return (-b+dis)/(2*a)
end

################################################################################
######################## AREA (GREEN'S THEOREM) ###############################
################################################################################

# Twice the oriented area contribution of the planar curve `crv`,
# ∫(x y' - y x')dt, used by `area` via Green's theorem.
@inline function _area_integral(crv::AbsCurve; rtol=1e-10)
    T = typeof(crv.length)
    f(t) = begin
        r = curve(crv, t)
        dr = tangent(crv, t)
        r[1]*dr[2] - r[2]*dr[1]
    end
    I, _ = quadgk(f, zero(T), one(T); rtol=rtol)
    return I
end

"""
    area(crv::AbsCurve; rtol=1e-10) → A::Real

Returns the geometric area enclosed by the planar boundary curve `crv`, via
Green's theorem `A = (1/2)|∫(x y' - y x')dt|`. The absolute value makes the
result independent of the boundary parametrization's orientation.
"""
@inline area(crv::AbsCurve; rtol=1e-10) = abs(_area_integral(crv; rtol=rtol))/2

"""
    area(crv::CompositeCurve; rtol=1e-10) → A::Real

Returns the geometric area enclosed by the composite boundary `crv`, summing
the Green-theorem contributions of every constituent curve before taking the
absolute value (preserving cancellation between oppositely oriented boundary
components).
"""
function area(crv::CompositeCurve; rtol=1e-10)
    T = typeof(crv.length)
    I = zero(T)
    @inbounds for subcrv in crv.subcurves
        I += _area_integral(subcrv; rtol=rtol)
    end
    return abs(I)/2
end

"""
    area(curves::AbstractVector{<:AbsCurve}; rtol=1e-10) → A::Real

Returns the geometric area enclosed by a collection of planar boundary
curves, summing oriented Green-theorem contributions before taking the
absolute value.
"""
function area(curves::AbstractVector{<:AbsCurve}; rtol=1e-10)
    isempty(curves) && return 0.0
    T = typeof(first(curves).length)
    I = zero(T)
    @inbounds for crv in curves
        I += _area_integral(crv; rtol=rtol)
    end
    return abs(I)/2
end

"""
    area(billiard::AbsBilliard; kwargs...) → A::Real

Returns the area of the complete physical billiard, using
[`BilliardGeometry.full_boundary`](@ref) so the result is independent of any
symmetry reduction represented by the billiard's fundamental domain.
"""
@inline area(billiard::AbsBilliard; kwargs...) = area(full_boundary(billiard); kwargs...)

# How many-fold a discrete symmetry divides the physical area onto the
# fundamental domain (used to estimate the fundamental-domain area below,
# without needing a separate field on the billiard struct).
@inline symmetry_reduction_factor(::AbsSymmetry) = 1
@inline symmetry_reduction_factor(::BilliardGeometry.XAxisReflection) = 2
@inline symmetry_reduction_factor(::BilliardGeometry.YAxisReflection) = 2
@inline symmetry_reduction_factor(::BilliardGeometry.XYAxisReflection) = 4
@inline symmetry_reduction_factor(::BilliardGeometry.DiagonalReflection) = 2
@inline symmetry_reduction_factor(::BilliardGeometry.AntiDiagonalReflection) = 2
@inline symmetry_reduction_factor(sym::BilliardGeometry.NFoldRotation) = sym.order
# Composite reflection's fundamental-domain reduction factor is not tracked
# per-irrep here; only used as a Weyl-window size estimate in BeynSolver, so
# a conservative (unreduced) factor of 1 never invalidates the windowing.
@inline symmetry_reduction_factor(::BilliardGeometry.CompositeReflection) = 1

function maximal_symmetry(billiard::AbsBilliard)
    isempty(billiard.symmetries) && return nothing
    _, i = findmax(symmetry_reduction_factor, billiard.symmetries)
    return billiard.symmetries[i]
end

@inline function symmetry_reduction_factor(billiard::AbsBilliard)
    sym = maximal_symmetry(billiard)
    return isnothing(sym) ? 1 : symmetry_reduction_factor(sym)
end

"""
    fundamental_area(billiard::AbsBilliard) → A::Real

Returns the area of the billiard's fundamental domain, obtained by dividing
the full physical [`area`](@ref) by the highest-order discrete symmetry's
[`symmetry_reduction_factor`](@ref).
"""
@inline function fundamental_area(billiard::AbsBilliard)
    return area(billiard)/symmetry_reduction_factor(billiard)
end