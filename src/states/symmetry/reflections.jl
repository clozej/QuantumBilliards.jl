@inline function apply_symmetries_to_wavefunction(Psi::AbstractMatrix{U},x_grid::AbstractVector{T},y_grid::AbstractVector{T},::Nothing) where {U<:Number,T<:Real};return Psi,x_grid,y_grid;end

"""
    apply_symmetries_to_wavefunction(Psi::AbstractMatrix{U},x_grid::AbstractVector{T},y_grid::AbstractVector{T},symmetry::BilliardGeometry.XAxisReflection) where {U<:Number,T<:Real} → Tuple

Reflect a wavefunction matrix across the `x` axis.

## Arguments
* `Psi::AbstractMatrix{U}`: Wavefunction on the fundamental-domain grid.
* `x_grid::AbstractVector{T}`: Fundamental-domain `x` coordinates.
* `y_grid::AbstractVector{T}`: Fundamental-domain `y` coordinates.
* `symmetry::BilliardGeometry.XAxisReflection`: Active reflection symmetry.

## Returns
* `full_Psi`: Wavefunction reflected across the `x` axis.
* `full_x`: Unchanged `x` grid.
* `full_y`: Symmetry-expanded `y` grid.
"""
function apply_symmetries_to_wavefunction(Psi::AbstractMatrix{U},x_grid::AbstractVector{T},y_grid::AbstractVector{T},symmetry::BilliardGeometry.XAxisReflection) where {U<:Number,T<:Real}
    p=symmetry.parity_y
    Psi_ref=reverse(p.*Psi,dims=2)
    full_y=vcat(-reverse(y_grid),y_grid)
    return hcat(Psi_ref,Psi),x_grid,full_y
end

"""
    apply_symmetries_to_wavefunction(Psi::AbstractMatrix{U},x_grid::AbstractVector{T},y_grid::AbstractVector{T},symmetry::YAxisReflection) where {U<:Number,T<:Real} → Tuple

Reflect a wavefunction matrix across the `y` axis.

## Arguments
* `Psi::AbstractMatrix{U}`: Wavefunction on the fundamental-domain grid.
* `x_grid::AbstractVector{T}`: Fundamental-domain `x` coordinates.
* `y_grid::AbstractVector{T}`: Fundamental-domain `y` coordinates.
* `symmetry::BilliardGeometry.YAxisReflection`: Active reflection symmetry.

## Returns
* `full_Psi`: Wavefunction reflected across the `y` axis.
* `full_x`: Symmetry-expanded `x` grid.
* `full_y`: Unchanged `y` grid.
"""
function apply_symmetries_to_wavefunction(Psi::AbstractMatrix{U},x_grid::AbstractVector{T},y_grid::AbstractVector{T},symmetry::BilliardGeometry.YAxisReflection) where {U<:Number,T<:Real}
    p=symmetry.parity_x
    Psi_ref=reverse(p.*Psi,dims=1)
    full_x=vcat(-reverse(x_grid),x_grid)
    return vcat(Psi_ref,Psi),full_x,y_grid
end

"""
    apply_symmetries_to_wavefunction(Psi::AbstractMatrix{U},x_grid::AbstractVector{T},y_grid::AbstractVector{T},symmetry::BilliardGeometry.XYAxisReflection) where {U<:Number,T<:Real} → Tuple

Reflect a wavefunction from the first-quadrant fundamental domain onto all four
quadrants.

The stored parities are

```text
parity_x : x → -x, reflection across the y axis,
parity_y : y → -y, reflection across the x axis.
```

With matrix dimension `1` corresponding to `x` and dimension `2` to `y`, the
four quadrants are assembled as

```text
Q3 | Q2
-------
Q4 | Q1
```

with `Q1` the supplied wavefunction.

## Arguments
* `Psi::AbstractMatrix{U}`: Wavefunction on the first-quadrant fundamental grid.
* `x_grid::AbstractVector{T}`: Positive-side `x` coordinates.
* `y_grid::AbstractVector{T}`: Positive-side `y` coordinates.
* `symmetry::BilliardGeometry.XYAxisReflection`: Two-axis reflection symmetry.

## Returns
* `full_Psi`: Wavefunction unfolded onto all four quadrants.
* `full_x`: Symmetry-expanded `x` grid.
* `full_y`: Symmetry-expanded `y` grid.
"""
function apply_symmetries_to_wavefunction(Psi::AbstractMatrix{U},x_grid::AbstractVector{T},y_grid::AbstractVector{T},symmetry::BilliardGeometry.XYAxisReflection) where {U<:Number,T<:Real}
    px=symmetry.parity_x
    py=symmetry.parity_y
    Psi_Y=reverse(px.*Psi,dims=1)
    Psi_X=reverse(py.*Psi,dims=2)
    Psi_XY=reverse((px*py).*Psi,dims=(1,2))
    full_Psi=hcat(vcat(Psi_XY,Psi_X),vcat(Psi_Y,Psi))
    full_x=vcat(-reverse(x_grid),x_grid)
    full_y=vcat(-reverse(y_grid),y_grid)
    return full_Psi,full_x,full_y
end

@inline function apply_symmetries_to_boundary_function(u::AbstractVector{U},::Nothing) where {U<:Number};return u;end

"""
    apply_symmetries_to_boundary_function(u::AbstractVector{U},symmetry::BilliardGeometry.XAxisReflection) where {U<:Number} → Vector

Unfold boundary data across the `x` axis.

The reflected copy is order-reversed to preserve the physical boundary
traversal and multiplied by `symmetry.parity_y`.

## Arguments
* `u::AbstractVector{U}`: Boundary data on the fundamental arc.
* `symmetry::BilliardGeometry.XAxisReflection`: Active reflection symmetry.

## Returns
* `full_u::Vector`: Symmetry-expanded boundary data.
"""
function apply_symmetries_to_boundary_function(u::AbstractVector{U},symmetry::BilliardGeometry.XAxisReflection) where {U<:Number}
    return vcat(u,symmetry.parity_y.*reverse(u))
end

"""
    apply_symmetries_to_boundary_function(u::AbstractVector{U},symmetry::BilliardGeometry.YAxisReflection) where {U<:Number} → Vector

Unfold boundary data across the `y` axis.

The reflected copy is order-reversed to preserve the physical boundary
traversal and multiplied by `symmetry.parity_x`.

## Arguments
* `u::AbstractVector{U}`: Boundary data on the fundamental arc.
* `symmetry::BilliardGeometry.YAxisReflection`: Active reflection symmetry.

## Returns
* `full_u::Vector`: Symmetry-expanded boundary data.
"""
function apply_symmetries_to_boundary_function(u::AbstractVector{U},symmetry::BilliardGeometry.YAxisReflection) where {U<:Number}
    return vcat(u,symmetry.parity_x.*reverse(u))
end

"""
    apply_symmetries_to_boundary_function(u::AbstractVector{U},symmetry::BilliardGeometry.XYAxisReflection) where {U<:Number} → Vector

Unfold boundary data from one quadrant onto the complete four-quadrant physical
boundary.

The copies are appended in continuous boundary order

```text
Q1 → Q2 → Q3 → Q4
```

where `Q2` and `Q4` are orientation-reversing reflected copies and `Q3` is the
orientation-preserving two-axis image. The factors are read directly from
`symmetry.parity_x` and `symmetry.parity_y`.

## Arguments
* `u::AbstractVector{U}`: Boundary data on the fundamental arc.
* `symmetry::BilliardGeometry.XYAxisReflection`: Two-axis reflection symmetry.

## Returns
* `full_u::Vector`: Boundary data unfolded onto the complete physical boundary.
"""
function apply_symmetries_to_boundary_function(u::AbstractVector{U},symmetry::BilliardGeometry.XYAxisReflection) where {U<:Number}
    px=symmetry.parity_x
    py=symmetry.parity_y
    return vcat(u,px.*reverse(u),(px*py).*u,py.*reverse(u))
end

"""
    apply_symmetries_to_boundary_points(pts::BoundaryPoints{T},symmetry::Nothing) where {T<:Real} → BoundaryPoints{T}

Return fundamental-domain boundary points unchanged when no reflection symmetry
is active.

## Arguments
* `pts::BoundaryPoints{T}`: Fundamental-domain boundary discretization.
* `symmetry::Nothing`: No active reflection symmetry.

## Returns
* `pts::BoundaryPoints{T}`: Input boundary discretization.
"""
@inline function apply_symmetries_to_boundary_points(pts::BoundaryPoints{T},::Nothing) where {T<:Real}
    return pts
end

"""
    apply_symmetries_to_boundary_points(pts::BoundaryPoints{T},symmetry::BilliardGeometry.XAxisReflection) where {T<:Real} → BoundaryPoints{T}

Unfold boundary points across the `x` axis.

Positions and outward normals are reflected geometrically and the reflected
copy is order-reversed to preserve the physical boundary traversal. The
quadrature weights `ds` are reversed with the corresponding nodes.

The returned [`BoundaryPoints`](@ref) is intended for basis-solver
post-processing and therefore populates only the full physical positions,
normals, arclength coordinates and physical quadrature weights.

## Arguments
* `pts::BoundaryPoints{T}`: Boundary points on the fundamental arc.
* `symmetry::BilliardGeometry.XAxisReflection`: Active reflection symmetry.

## Returns
* `full_pts::BoundaryPoints{T}`: Symmetry-expanded physical boundary points.
"""
function apply_symmetries_to_boundary_points(pts::BoundaryPoints{T},symmetry::BilliardGeometry.XAxisReflection) where {T<:Real}
    xy_ref=reverse(BilliardGeometry.apply_symmetry(symmetry,pts.xy))
    normal_ref=reverse(BilliardGeometry.apply_symmetry(symmetry,pts.normal))
    full_xy=vcat(pts.xy,xy_ref)
    full_normal=vcat(pts.normal,normal_ref)
    full_ds=vcat(pts.ds,reverse(pts.ds))
    L=sum(pts.ds)
    s0=pts.s.-pts.s[1]
    full_s=vcat(s0,L.+L.-reverse(s0))
    return BoundaryPoints(full_xy;normal=full_normal,s=full_s,ds=full_ds)
end

"""
    apply_symmetries_to_boundary_points(pts::BoundaryPoints{T},symmetry::BilliardGeometry.YAxisReflection) where {T<:Real} → BoundaryPoints{T}

Unfold boundary points across the `y` axis.

Positions and outward normals are reflected geometrically and the reflected
copy is order-reversed to preserve the physical boundary traversal. The
quadrature weights `ds` are reversed with the corresponding nodes.

The returned [`BoundaryPoints`](@ref) is intended for basis-solver
post-processing and therefore populates only the full physical positions,
normals, arclength coordinates and physical quadrature weights.

## Arguments
* `pts::BoundaryPoints{T}`: Boundary points on the fundamental arc.
* `symmetry::BilliardGeometry.YAxisReflection`: Active reflection symmetry.

## Returns
* `full_pts::BoundaryPoints{T}`: Symmetry-expanded physical boundary points.
"""
function apply_symmetries_to_boundary_points(pts::BoundaryPoints{T},symmetry::BilliardGeometry.YAxisReflection) where {T<:Real}
    xy_ref=reverse(BilliardGeometry.apply_symmetry(symmetry,pts.xy))
    normal_ref=reverse(BilliardGeometry.apply_symmetry(symmetry,pts.normal))
    full_xy=vcat(pts.xy,xy_ref)
    full_normal=vcat(pts.normal,normal_ref)
    full_ds=vcat(pts.ds,reverse(pts.ds))
    L=sum(pts.ds)
    s0=pts.s.-pts.s[1]
    full_s=vcat(s0,L.+L.-reverse(s0))
    return BoundaryPoints(full_xy;normal=full_normal,s=full_s,ds=full_ds)
end

"""
    apply_symmetries_to_boundary_points(pts::BoundaryPoints{T},symmetry::BilliardGeometry.XYAxisReflection) where {T<:Real} → BoundaryPoints{T}

Unfold boundary points from one quadrant onto the complete four-quadrant
physical boundary.

The physical copies are appended in continuous boundary order

```text
Q1 → Q2 → Q3 → Q4
```

`Q2` is obtained by reflection across the `y` axis and order reversal, `Q3`
by the combined two-axis reflection without order reversal, and `Q4` by
reflection across the `x` axis and order reversal. The same ordering is used
by [`apply_symmetries_to_boundary_function`](@ref).

The returned [`BoundaryPoints`](@ref) is intended for basis-solver
post-processing and therefore populates only the full physical positions,
normals, arclength coordinates and physical quadrature weights.

## Arguments
* `pts::BoundaryPoints{T}`: Boundary points on the fundamental arc.
* `symmetry::BilliardGeometry.XYAxisReflection`: Two-axis reflection symmetry.

## Returns
* `full_pts::BoundaryPoints{T}`: Boundary points unfolded onto the complete physical boundary.
"""
function apply_symmetries_to_boundary_points(pts::BoundaryPoints{T},symmetry::BilliardGeometry.XYAxisReflection) where {T<:Real}
    sym_y=BilliardGeometry.YAxisReflection(symmetry.parity_x)
    sym_x=BilliardGeometry.XAxisReflection(symmetry.parity_y)
    xyY=reverse(BilliardGeometry.apply_symmetry(sym_y,pts.xy))
    normalY=reverse(BilliardGeometry.apply_symmetry(sym_y,pts.normal))
    xyXY=BilliardGeometry.apply_symmetry(symmetry,pts.xy)
    normalXY=BilliardGeometry.apply_symmetry(symmetry,pts.normal)
    xyX=reverse(BilliardGeometry.apply_symmetry(sym_x,pts.xy))
    normalX=reverse(BilliardGeometry.apply_symmetry(sym_x,pts.normal))
    full_xy=vcat(pts.xy,xyY,xyXY,xyX)
    full_normal=vcat(pts.normal,normalY,normalXY,normalX)
    full_ds=vcat(pts.ds,reverse(pts.ds),pts.ds,reverse(pts.ds))
    L=sum(pts.ds)
    s0=pts.s.-pts.s[1]
    full_s=vcat(s0,L.+L.-reverse(s0),2*L.+s0,3*L.+L.-reverse(s0))
    return BoundaryPoints(full_xy;normal=full_normal,s=full_s,ds=full_ds)
end
