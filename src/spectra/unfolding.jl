"""
    corner_correction(corner_angles)

Return the constant corner correction entering the two-dimensional Weyl law,

    C = Σᵢ (π²-αᵢ²)/(24παᵢ),

where `αᵢ` are the interior corner angles.
"""
@inline function corner_correction(corner_angles)
    return sum(c -> (pi^2-c^2)/(24*pi*c),corner_angles)
end

"""
    weyl_law(k,A,L)

Return the leading two-dimensional Weyl counting function

    N(k) = (A k² - L k)/(4π),

where `A` is the billiard area and `L` its perimeter.
"""
@inline weyl_law(k,A,L)=@. (A*k^2-L*k)/(4*pi)

"""
    weyl_law(k,A,L,corner_angles)

Return the two-dimensional Weyl counting function including the constant corner
correction,

    N(k) = (A k² - L k)/(4π) + C,

where

    C = Σᵢ (π²-αᵢ²)/(24παᵢ).
"""
@inline function weyl_law(k,A,L,corner_angles)
    return weyl_law(k,A,L).+corner_correction(corner_angles)
end

"""
    k_at_state(state,A,L)

Invert the leading two-dimensional Weyl law and return the positive wavenumber
corresponding to the counting index `state`.

The inversion is

    k = (L + √(L² + 16πA state))/(2A).
"""
@inline function k_at_state(state,A,L)
    disc=L^2+16*pi*A*state
    return (L+sqrt(disc))/(2*A)
end

"""
    k_at_state(state,A,L,corner_angles)

Invert the two-dimensional Weyl law including the constant corner correction
and return the positive wavenumber corresponding to the counting index `state`.

If `C = corner_correction(corner_angles)`, the inversion is

    k = (L + √(L² + 16πA(state-C)))/(2A).
"""
@inline function k_at_state(state,A,L,corner_angles)
    C=corner_correction(corner_angles)
    disc=L^2+16*pi*A*(state-C)
    return (L+sqrt(disc))/(2*A)
end

"""
    _area_integral(crv::AbsCurve;rtol=1e-10)

Return twice the oriented area contribution of the planar curve `crv`,

    ∫ (x y' - y x') dt,

as used in Green's theorem.

This is an internal helper used by [`area`](@ref).
"""
@inline function _area_integral(crv::AbsCurve;rtol=1e-10)
    T=typeof(crv.length)
    f(t)=begin
        r=BilliardGeometry.curve(crv,t)
        dr=tangent(crv,t)
        r[1]*dr[2]-r[2]*dr[1]
    end
    I,_=quadgk(f,zero(T),one(T);rtol=rtol)
    return I
end

"""
    area(crv::AbsCurve;rtol=1e-10)

Return the geometric area enclosed by the planar boundary curve `crv`.

The area is evaluated using Green's theorem,

    A = 1/2 |∫ (x y' - y x') dt|.

The absolute value makes the result independent of the orientation of the
boundary parametrization.
"""
@inline function area(crv::BilliardGeometry.AbsCurve;rtol=1e-10)
    return abs(_area_integral(crv;rtol=rtol))/2
end

"""
    area(crv::BilliardGeometry.CompositeCurve;rtol=1e-10)

Return the geometric area enclosed by the composite boundary `crv`.

The Green-theorem contributions of all constituent curves are summed before the
absolute value is taken,

    A = 1/2 |Σⱼ ∫ (xⱼ yⱼ' - yⱼ xⱼ') dt|.

Summing before taking the absolute value preserves the correct cancellation
between oppositely oriented boundary components.
"""
function area(crv::BilliardGeometry.CompositeCurve;rtol=1e-10)
    T=typeof(crv.length)
    I=zero(T)
    @inbounds for subcrv in crv.subcurves
        I+=_area_integral(subcrv;rtol=rtol)
    end
    return abs(I)/2
end

"""
    area(curves::AbstractVector{<:BilliardGeometry.AbsCurve};rtol=1e-10)

Return the geometric area enclosed by a collection of planar boundary curves.
The oriented Green-theorem contributions are summed before taking the absolute
value so that oppositely oriented boundary components cancel correctly.
"""
function area(curves::AbstractVector{<:BilliardGeometry.AbsCurve};rtol=1e-10)
    isempty(curves)&&return 0.0
    T=typeof(first(curves).length)
    I=zero(T)
    @inbounds for crv in curves
        I+=_area_integral(crv;rtol=rtol)
    end
    return abs(I)/2
end

"""
    area(components::AbstractVector{<:AbstractVector{<:BilliardGeometry.AbsCurve}};rtol=1e-10)

Return the geometric area of a multiply connected planar domain.

The first connected boundary component is interpreted as the outer boundary and
all subsequent components as holes. The area of each component is evaluated
independently using Green's theorem and the hole areas are subtracted,

    A = A_outer - Σₕ A_h.

## Arguments
* `components::AbstractVector{<:AbstractVector{<:BilliardGeometry.AbsCurve}}`: Connected physical boundary components, with the outer boundary first and holes following.

## Keyword Arguments
* `rtol::Real=1e-10`: Relative tolerance used by the Green-theorem quadrature.

## Returns
* `A`: Physical area enclosed by the outer boundary with all hole areas removed.
"""
function area(components::AbstractVector{<:AbstractVector{<:BilliardGeometry.AbsCurve}};rtol=1e-10)
    isempty(components)&&return 0.0
    A=area(first(components);rtol=rtol)
    @inbounds for comp in @view components[2:end]
        A-=area(comp;rtol=rtol)
    end
    return A
end

"""
    area(billiard::BilliardGeometry.AbsBilliard;kwargs...)

Return the area of the complete physical billiard.

The calculation uses `billiard.full_boundary`, so the returned area is
independent of any symmetry reduction represented by the billiard's fundamental
domain.
"""
@inline area(billiard::BilliardGeometry.AbsBilliard;kwargs...)=area(billiard.full_boundary;kwargs...)

# Nice ways to get the fundamental domain size without needing a separate field for it in the billiard struct. Assign to each symmetry in the vector of classical dynamics symmetries and find the one which is the highest
# Then just divide the area by that number to get the fundamental domain area. This is used in Beyn (and potentially in EBIM if one wants constant level per interval). 
@inline symmetry_reduction_factor(::BilliardGeometry.AbsSymmetry)=1
@inline symmetry_reduction_factor(::BilliardGeometry.XAxisReflection)=2
@inline symmetry_reduction_factor(::BilliardGeometry.YAxisReflection)=2
@inline symmetry_reduction_factor(::BilliardGeometry.XYAxisReflection)=4
@inline symmetry_reduction_factor(::DiagonalReflection)=2
@inline symmetry_reduction_factor(::AntiDiagonalReflection)=2
@inline symmetry_reduction_factor(sym::BilliardGeometry.NFoldRotation)=sym.order
#TODO symmetry reduction factor for composite reflection symmetry. No matter, since this always reduces A to the max and since 
# in Beyn this is just an estimate for the number of eigvals inside the contour does not change anything
@inline symmetry_reduction_factor(::CompositeReflection)=1 #FIX THIS

function maximal_symmetry(billiard::BilliardGeometry.AbsBilliard)
    isempty(billiard.symmetries)&&return nothing
    _,i=findmax(symmetry_reduction_factor,billiard.symmetries)
    return billiard.symmetries[i]
end

@inline function symmetry_reduction_factor(billiard::BilliardGeometry.AbsBilliard)
    sym=maximal_symmetry(billiard)
    return isnothing(sym) ? 1 : symmetry_reduction_factor(sym)
end

@inline function fundamental_area(billiard::BilliardGeometry.AbsBilliard)
    return area(billiard)/symmetry_reduction_factor(billiard)
end