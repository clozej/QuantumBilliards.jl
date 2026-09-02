struct CartesianCS{T} <:CoordinateSystem where {T<:Number}
    origin::SVector{2,T}
    rot_angle::T
    affine_map::CoordinateTransformations.AffineMap{Rotations.Angle2d{T},SVector{2,T}}
    local_map::CoordinateTransformations.AffineMap{Rotations.Angle2d{T},SVector{2,T}}
end

"""
    CartesianCS(origin::SVector{2,T},rot_angle::T) where {T<:Number}

Create a `CartesianCS` struct that contains information about transformations possible in that coordinate transformation based on the usage choice. One can choose between `PolarCS` and `CartesianCS`. It contains maps that handle translations and rotations of points in the respective coordinate system along with their inverses.

# Arguments
- `origin::SVector{2,T}`: Origin of the coordinate system.
- `rot_angle::T`: The given rotation of the coordinate system.

# Returns
- `CartesianCS`: The struct containing the `AffineMap` for the rotations/translations and inverses for this coordiante system choice.
"""
function CartesianCS(origin::SVector{2,T},rot_angle::T) where {T<:Number}
    Rot=CoordinateTransformations.LinearMap(Rotations.Angle2d(rot_angle))
    Tran=CoordinateTransformations.Translation(origin[1],origin[2])
    Tran_inv=CoordinateTransformations.Translation(-origin[1],-origin[2])
    Rot_inv=CoordinateTransformations.LinearMap(Rotations.Angle2d(-rot_angle))
    affine_map=CoordinateTransformations.compose(Tran,Rot)
    local_map=CoordinateTransformations.compose(Rot_inv,Tran_inv)
    return CartesianCS(origin,rot_angle,affine_map,local_map)
end

struct PolarCS{T} <:CoordinateSystem  where {T<:Number}
    origin::SVector{2,T}
    rot_angle::T
    affine_map::CoordinateTransformations.AffineMap{Rotations.Angle2d{T},SVector{2, T}}  #maps carthesian coordinates
    local_map::CoordinateTransformations.AffineMap{Rotations.Angle2d{T},SVector{2, T}} #transform carthesian into local polar coords
end

"""
    PolarCS(origin::SVector{2,T},rot_angle::T) where {T<:Number}

Create a `PolarCS` struct that contains information about transformations possible in that coordinate transformation based on the usage choice. One can choose between `PolarCS` and `CartesianCS`. It contains maps that handle translations and rotations of points in the respective coordinate system along with their inverses.

# Arguments
- `origin::SVector{2,T}`: Origin of the coordinate system.
- `rot_angle::T`: The given rotation of the coordinate system.

# Returns
- `PolarCS`: The struct containing the `AffineMap` for the rotations/translations and inverses for this coordiante system choice.
"""
function PolarCS(origin::SVector{2,T},rot_angle::T) where {T<:Number}
    Rot=CoordinateTransformations.LinearMap(Rotations.Angle2d(rot_angle))
    Tran=CoordinateTransformations.Translation(origin[1],origin[2])
    Tran_inv=CoordinateTransformations.Translation(-origin[1],-origin[2])
    Rot_inv=CoordinateTransformations.LinearMap(Rotations.Angle2d(-rot_angle))
    affine_map=CoordinateTransformations.compose(Tran,Rot) #already in cartesian coordinates
    local_map=CoordinateTransformations.compose(Rot_inv,Tran_inv) # rotate in local polar coordinates
    return PolarCS(origin,rot_angle,affine_map,local_map)
end

"""
    linear_map(cs, v) -> SVector{2,T}

Apply only the linear (rotation) part of a coordinate-system transformation
to a vector `v`, discarding any translational component.

# Purpose
In the geometry framework, `AffineMap`s stored in `CartesianCS` and `PolarCS`
represent transformations of the form

    x ↦ R*x + t

where `R` is a rotation and `t` is a translation.

While this is correct for points, it is incorrect for vectors such as:
- tangents,
- normals,
- higher derivatives,

because vectors must transform only under the linear part `R`.

This helper extracts that linear action by removing the translation.

# Definition
The linear map is computed as

    linear_map(cs, v) = cs.affine_map(v) - cs.affine_map(0)

which is equivalent to applying the rotation `R` alone.

# Arguments
- `cs`: A coordinate system (`CartesianCS` or `PolarCS`)
- `v`: A vector in ℝ² (e.g. tangent, normal, derivative)

# Returns
- `SVector{2,T}`: The rotated vector

# Important
- This function must be used for transforming derivatives.
- Using `affine_map` directly on tangents introduces a translation offset,
  which leads to incorrect geometry (e.g. twisted normals for off-center curves).
"""
@inline linear_map(cs,v)=cs.affine_map(v)-cs.affine_map(zero(v))
@inline function linear_map(cs::Union{CartesianCS,PolarCS},v::SVector{2})
    s,c=sincos(cs.rot_angle)
    return SVector(c*v[1]-s*v[2],s*v[1]+c*v[2])
end

"""
    polar_to_cartesian(pt::SVector{2,T}) where {T<:Number}

Standard `(r,ϕ) -> (x,y)` transformation.

# Arguments
- `pt::SVector{2,T}`: A 2D point represented as a static vector `(r,ϕ)`.

# Returns
- `SVector{2,T}`: A 2D point represented as a static vector `(x,y)`.
"""
function polar_to_cartesian(pt::SVector{2,T}) where {T<:Number}
    s,c=sincos(pt[2]) # pt[2] = ϕ
    return SVector(pt[1]*c,pt[1]*s)
end
    
"""
    cartesian_to_polar(pt::SVector{2,T}; rotation_angle_discontinuity::T = zero(T)) where T<:Number

Convert a 2D Cartesian point to polar coordinates, with an optional rotation applied before conversion. This is neccesery when we have a discontinutiy in the basis functions since e.g. `atan` has a jump at `ϕ=π`. This means that we will need to shift all the evaluation points of the basis by a given angle (rotation_angle_discontinuity) such that the discontinutiy is outside of our domain.

# Arguments
- `pt::SVector{2,T}`: A 2D point represented as a static vector `(x, y)`.
- `rotation_angle_discontinuity::T`: (optional) An angle in radians. If non-zero, the input point is first rotated by this angle before conversion. The resulting polar angle `θ` is then adjusted by subtracting this rotation to maintain continuity. Default is `0`.

# Returns
- `SVector{2,T}`: A static vector `(r, θ)` where:
- `r`: radial distance from the origin.
- `θ`: polar angle in radians from the positive x-axis.
"""
function cartesian_to_polar(pt::SVector{2,T};rotation_angle_discontinuity::T=zero(T)) where {T<:Number}
    if rotation_angle_discontinuity!=zero(T)
        s,c=sincos(rotation_angle_discontinuity)
        xr=pt[1]*c-pt[2]*s
        yr=pt[1]*s+pt[2]*c
        r=hypot(xr,yr)
        θ=atan(yr,xr)
        return SVector(r,θ-rotation_angle_discontinuity)
    else
        r=hypot(pt[1],pt[2])
        θ=atan(pt[2],pt[1])
        return SVector(r,θ)
    end
end

# in-place polar (zero allocs)
@inline function _polar_coords!(r::AbstractVector,φ::AbstractVector,pm,pts,rotation_angle_discontinuity)
    θ0=rotation_angle_discontinuity
    s0,c0=sincos(θ0)
    @inbounds for j in eachindex(pts)
        x=pm(pts[j])
        if θ0==zero(θ0)
            xr=x[1]
            yr=x[2]
        else
            xr=x[1]*c0-x[2]*s0
            yr=x[1]*s0+x[2]*c0
        end
        r[j]=hypot(xr,yr)
        φ[j]=atan(yr,xr)-θ0
    end
    return r,φ
end