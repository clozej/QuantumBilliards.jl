"""
    RealPlaneWaves{T,Sa}<:AbsBasis

Concrete real plane-wave basis with an optional reflection symmetries.

Each basis function has the separable form

```math
f(x,y)=F_x(kv_xx)F_y(kv_yy),
```

where `F_x,F_y ∈ {cos,sin}` are selected by `parity_x` and `parity_y`
(`+1` selects `cos`, `-1` selects `sin`), and `(v_x,v_y)` is determined by
the corresponding entry of `angles`.

At most one reflection symmetry is stored:
- `nothing`: no symmetry; all four parity combinations are used.
- `XAxisReflection`: reflection across the x axis (`y→-y`), fixing y parity.
- `YAxisReflection`: reflection across the y axis (`x→-x`), fixing x parity.
- `XYAxisReflection`: both reflections, fixing both parities.

The symmetry object itself stores the relevant parity quantum number(s).

## Attributes
* `dim::Int`: Effective number of basis functions.
* `symmetries`: Reflection symmetry, or `nothing`.
* `angle_arc::T`: Angular range over which directions are sampled.
* `angle_shift::T`: Angular offset applied to the sampled directions.
* `angles::Vector{T}`: Propagation angles.
* `parity_x::Vector{Int}`: x-factor parity (`+1` = `cos`, `-1` = `sin`).
* `parity_y::Vector{Int}`: y-factor parity (`+1` = `cos`, `-1` = `sin`).
* `sampler::Sa`: Sampling strategy used to sample the angles.

## API
The following functions can be evaluated for this type:
- [`resize_basis`](@ref)
- [`basis_fun`](@ref)
- [`gradient`](@ref)
- [`basis_and_gradient`](@ref)
- [`dk_fun`](@ref)
"""
struct RealPlaneWaves{T<:Real,Sa<:AbsSampler}<:AbsBasis
    dim::Int
    symmetries::Union{BilliardGeometry.AbsReflection,Nothing}
    angle_arc::T
    angle_shift::T
    angles::Vector{T}
    parity_x::Vector{Int}
    parity_y::Vector{Int}
    sampler::Sa
end

"""
    parity_pattern(symmetry) → (parity_x,parity_y)

Return the cos/sin parity patterns selected by a single reflection symmetry.

`+1` selects `cos` and `-1` selects `sin`.
"""
@inline parity_pattern(::Nothing)=Int[1,1,-1,-1],Int[1,-1,1,-1]

@inline function parity_pattern(symmetry::BilliardGeometry.XAxisReflection)
    p=Int(symmetry.parity_y)
    return Int[1,-1],Int[p,p]
end

@inline function parity_pattern(symmetry::BilliardGeometry.YAxisReflection)
    p=Int(symmetry.parity_x)
    return Int[p,p],Int[1,-1]
end

@inline function parity_pattern(symmetry::BilliardGeometry.XYAxisReflection)
    return Int[Int(symmetry.parity_x)],Int[Int(symmetry.parity_y)]
end

"""
    RealPlaneWaves(dim::Int,symmetry::Union{BilliardGeometry.AbsReflection,Nothing}=nothing;angle_arc::Union{Real,Nothing}=nothing,angle_shift::Union{Real,Nothing}=nothing,sampler=LinearNodes()) → RealPlaneWaves

Construct a real plane-wave basis with at most one reflection symmetry.

`dim` is the number of sampled propagation angles before parity expansion.
The effective basis dimension is `dim*length(first(parity_pattern(symmetry)))`.

When `angle_arc` or `angle_shift` is omitted, defaults are chosen from the
symmetry:
- no symmetry: `angle_arc=π`, `angle_shift=0`;
- `XAxisReflection`: `angle_arc=π`, `angle_shift=0`;
- `YAxisReflection`: `angle_arc=π`, `angle_shift=-π/2`;
- `XYAxisReflection`: `angle_arc=π/2`, `angle_shift=0`.

## Arguments
* `dim::Int`: Number of propagation angles to sample.
* `symmetry`: Reflection symmetry, or `nothing`.

## Keyword Arguments
* `angle_arc::Union{Real,Nothing}=nothing`: Angular sampling range.
* `angle_shift::Union{Real,Nothing}=nothing`: Angular sampling offset.
* `sampler::AbsSampler=LinearNodes()`: Angular sampler.

## Returns
* `basis::RealPlaneWaves`: Constructed basis.
"""
function RealPlaneWaves(dim::Int,symmetry::Union{BilliardGeometry.AbsReflection,Nothing};angle_arc::Union{Real,Nothing}=nothing,angle_shift::Union{Real,Nothing}=nothing,sampler=BilliardGeometry.LinearNodes())
    dim>0||throw(ArgumentError("dim must be positive"))
    par_x,par_y=parity_pattern(symmetry)
    pl=length(par_x)
    eff_dim=dim*pl
    default_arc=symmetry isa BilliardGeometry.XYAxisReflection ? π/2 : π
    default_shift=symmetry isa BilliardGeometry.YAxisReflection ? -π/2 : 0.0
    arc=isnothing(angle_arc) ? default_arc : angle_arc
    shift=isnothing(angle_shift) ? default_shift : angle_shift
    t,_=BilliardGeometry.sample_points(sampler,dim)
    T=eltype(t)
    arcT=T(arc);shiftT=T(shift)
    angles=Vector{T}(undef,eff_dim)
    parity_x=Vector{Int}(undef,eff_dim)
    parity_y=Vector{Int}(undef,eff_dim)
    @inbounds for i in 1:dim
        angle=t[i]*arcT+shiftT
        base=(i-1)*pl
        for j in 1:pl
            idx=base+j
            angles[idx]=angle
            parity_x[idx]=par_x[j]
            parity_y[idx]=par_y[j]
        end
    end
    return RealPlaneWaves{T,typeof(sampler)}(eff_dim,symmetry,arcT,shiftT,angles,parity_x,parity_y,sampler)
end

"""
    RealPlaneWaves(dim::Int;sym_x::Union{Int,Nothing}=nothing,sym_y::Union{Int,Nothing}=nothing,angle_arc::Union{Real,Nothing}=nothing,angle_shift::Union{Real,Nothing}=nothing,sampler=LinearNodes()) → RealPlaneWaves

Construct a real plane-wave basis by specifying reflection parities.

`sym_x` is the parity under reflection across the y axis (`x→-x`), while
`sym_y` is the parity under reflection across the x axis (`y→-y`).

The two parity values are represented by a single symmetry object:
- neither specified → `nothing`;
- only `sym_y` → `XAxisReflection(sym_y)`;
- only `sym_x` → `YAxisReflection(sym_x)`;
- both → `XYAxisReflection(sym_x,sym_y)`.

## Arguments
* `dim::Int`: Number of propagation angles to sample before parity expansion.

## Keyword Arguments
* `sym_x::Union{Int,Nothing}=nothing`: x parity (`±1`).
* `sym_y::Union{Int,Nothing}=nothing`: y parity (`±1`).
* `angle_arc::Union{Real,Nothing}=nothing`: Angular sampling range.
* `angle_shift::Union{Real,Nothing}=nothing`: Angular sampling offset.
* `sampler::AbsSampler=LinearNodes()`: Angular sampler.

## Returns
* `basis::RealPlaneWaves`: Constructed basis.
"""
function RealPlaneWaves(dim::Int;sym_x::Union{Int,Nothing}=nothing,sym_y::Union{Int,Nothing}=nothing,angle_arc::Union{Real,Nothing}=nothing,angle_shift::Union{Real,Nothing}=nothing,sampler=BilliardGeometry.LinearNodes())
    isnothing(sym_x)||(sym_x==1||sym_x==-1)||throw(ArgumentError("sym_x must be ±1 or nothing"))
    isnothing(sym_y)||(sym_y==1||sym_y==-1)||throw(ArgumentError("sym_y must be ±1 or nothing"))
    symmetry=if isnothing(sym_x)&&isnothing(sym_y)
        nothing
    elseif isnothing(sym_x)
        BilliardGeometry.XAxisReflection(sym_y)
    elseif isnothing(sym_y)
        BilliardGeometry.YAxisReflection(sym_x)
    else
        BilliardGeometry.XYAxisReflection(sym_x,sym_y)
    end
    return RealPlaneWaves(dim,symmetry;angle_arc=angle_arc,angle_shift=angle_shift,sampler=sampler)
end

"""
    resize_basis(basis::RealPlaneWaves,billiard::AbsBilliard,dim::Int,k) → RealPlaneWaves

Resize a [`RealPlaneWaves`](@ref) basis while preserving its symmetry and
angular sampling parameters.

## Arguments
* `basis::RealPlaneWaves`: Basis to resize.
* `billiard::AbsBilliard`: Billiard associated with the basis.
* `dim::Int`: Number of sampled propagation angles before parity expansion.
* `k`: Wavenumber.

## Returns
* `basis_new::RealPlaneWaves`: Resized basis.
"""
@inline function resize_basis(basis::RealPlaneWaves,billiard::AbsBilliard,dim::Int,k)
    return RealPlaneWaves(dim,basis.symmetries;angle_arc=basis.angle_arc,angle_shift=basis.angle_shift,sampler=basis.sampler)
end

"""
    rescale_dimension(basis::Ba,dim::Integer) where {Ba<:AbsBasis} → Int

Convert an effective basis dimension to the number of sampled directions
required by basis constructors that internally expand parity sectors.

For [`RealPlaneWaves`](@ref), the multiplicity is four without symmetry, two
for a single-axis reflection, and one for `XYAxisReflection`.

## Arguments
* `basis::Ba`: Basis.
* `dim::Integer`: Effective basis dimension.

## Returns
* `dim::Int`: Dimension to pass to `resize_basis`.
"""
@inline function rescale_dimension(basis::Ba,dim::Integer) where {Ba<:AbsBasis}
    basis isa RealPlaneWaves||return Int(dim)
    symmetries=basis.symmetries
    multiplicity=isnothing(symmetries) ? 4 : symmetries isa BilliardGeometry.XYAxisReflection ? 1 : 2
    return div(Int(dim),multiplicity)
end

# Helper functions for cos/sin pattern
# parity = 1 → cos, parity = -1 → sin
@inline _cos(arg)=cos(arg)
@inline _sin(arg)=sin(arg)
@inline _rpw_fun(par::Int)=par==1 ? _cos : _sin
@inline _drpw_fun(par::Int)=par==1 ? (x->-sin(x)) : _cos

"""
    basis_fun(basis::RealPlaneWaves,i::Int,k::T,pts::AbstractArray) where {T<:Real} → out::Vector{T}

Evaluate the `i`-th real plane-wave basis function at wavenumber `k` on
`pts`.

## Arguments
* `basis::RealPlaneWaves`: Basis.
* `i::Int`: Basis-function index.
* `k::T`: Wavenumber.
* `pts::AbstractArray`: Evaluation points.

## Returns
* `out::Vector{T}`: Values of basis function `i`.
"""
@inline function basis_fun(basis::RealPlaneWaves,i::Int,k::T,pts::AbstractArray) where {T<:Real}
    parx=basis.parity_x[i]
    pary=basis.parity_y[i]
    vx=cos(basis.angles[i])
    vy=sin(basis.angles[i])
    fx=_rpw_fun(parx)
    fy=_rpw_fun(pary)
    M=length(pts)
    out=Vector{T}(undef,M)
    @inbounds @simd for j=1:M
        x=pts[j][1]
        y=pts[j][2]
        out[j]=fx(k*vx*x)*fy(k*vy*y)
    end
    return out
end

"""
    basis_fun(basis::RealPlaneWaves,indices::AbstractArray,k::T,pts::AbstractArray;multithreaded::Bool=true) where {T<:Real} → B::Matrix{T}

Evaluate selected real plane-wave basis functions at wavenumber `k`.

## Arguments
* `basis::RealPlaneWaves`: Basis.
* `indices::AbstractArray`: Basis-function indices.
* `k::T`: Wavenumber.
* `pts::AbstractArray`: Evaluation points.

## Keyword Arguments
* `multithreaded::Bool=true`: Construct columns in parallel.

## Returns
* `B::Matrix{T}`: Basis matrix.
"""
@inline function basis_fun(basis::RealPlaneWaves,indices::AbstractArray,k::T,pts::AbstractArray;multithreaded::Bool=true) where {T<:Real}
    M=length(pts)
    N=length(indices)
    B=Matrix{T}(undef,M,N)
    @use_threads multithreading=multithreaded for c in 1:N
        idx=indices[c]
        parx=basis.parity_x[idx]
        pary=basis.parity_y[idx]
        vx=cos(basis.angles[idx])
        vy=sin(basis.angles[idx])
        fx=_rpw_fun(parx)
        fy=_rpw_fun(pary)
        col=@view B[:,c]
        @inbounds @simd for j=1:M
            x=pts[j][1]
            y=pts[j][2]
            col[j]=fx(k*vx*x)*fy(k*vy*y)
        end
    end
    return B
end

"""
    gradient(basis::RealPlaneWaves,i::Int,k::T,pts::AbstractArray) where {T<:Real} → (dx,dy)

Evaluate the spatial gradient of the `i`-th real plane-wave basis function.

## Arguments
* `basis::RealPlaneWaves`: Basis.
* `i::Int`: Basis-function index.
* `k::T`: Wavenumber.
* `pts::AbstractArray`: Evaluation points.

## Returns
* `(dx,dy)`: x and y derivatives.
"""
function gradient(basis::RealPlaneWaves,i::Int,k::T,pts::AbstractArray) where {T<:Real}
    parx=basis.parity_x[i]
    pary=basis.parity_y[i]
    vx=cos(basis.angles[i])
    vy=sin(basis.angles[i])
    fx=_rpw_fun(parx)
    fy=_rpw_fun(pary)
    dfx=_drpw_fun(parx)
    dfy=_drpw_fun(pary)
    M=length(pts)
    dx=Vector{T}(undef,M)
    dy=Vector{T}(undef,M)
    @inbounds @simd for j=1:M
        x=pts[j][1]
        y=pts[j][2]
        ax=k*vx*x
        ay=k*vy*y
        bx=fx(ax)
        by=fy(ay)
        dx[j]=k*vx*dfx(ax)*by
        dy[j]=bx*k*vy*dfy(ay)
    end
    return dx,dy
end

"""
    gradient(basis::RealPlaneWaves,indices::AbstractArray,k::T,pts::AbstractArray;multithreaded::Bool=true) where {T<:Real} → (dB_dx,dB_dy)

Evaluate spatial gradients of selected real plane-wave basis functions.

## Arguments
* `basis::RealPlaneWaves`: Basis.
* `indices::AbstractArray`: Basis-function indices.
* `k::T`: Wavenumber.
* `pts::AbstractArray`: Evaluation points.

## Keyword Arguments
* `multithreaded::Bool=true`: Construct columns in parallel.

## Returns
* `(dB_dx,dB_dy)`: x and y derivative matrices.
"""
function gradient(basis::RealPlaneWaves,indices::AbstractArray,k::T,pts::AbstractArray;multithreaded::Bool=true) where {T<:Real}
    M=length(pts)
    N=length(indices)
    dBdx=Matrix{T}(undef,M,N)
    dBdy=Matrix{T}(undef,M,N)
    @use_threads multithreading=multithreaded for c in 1:N
        idx=indices[c]
        parx=basis.parity_x[idx]
        pary=basis.parity_y[idx]
        vx=cos(basis.angles[idx])
        vy=sin(basis.angles[idx])
        fx=_rpw_fun(parx)
        fy=_rpw_fun(pary)
        dfx=_drpw_fun(parx)
        dfy=_drpw_fun(pary)
        cx=@view dBdx[:,c]
        cy=@view dBdy[:,c]
        @inbounds @simd for j=1:M
            x=pts[j][1]
            y=pts[j][2]
            ax=k*vx*x
            ay=k*vy*y
            bx=fx(ax)
            by=fy(ay)
            cx[j]=k*vx*dfx(ax)*by
            cy[j]=bx*k*vy*dfy(ay)
        end
    end
    return dBdx,dBdy
end

"""
    basis_and_gradient(basis::RealPlaneWaves,i::Int,k::T,pts::AbstractArray) where {T<:Real} → (bf,dx,dy)

Evaluate the `i`-th basis function and its spatial gradient in one pass.

## Arguments
* `basis::RealPlaneWaves`: Basis.
* `i::Int`: Basis-function index.
* `k::T`: Wavenumber.
* `pts::AbstractArray`: Evaluation points.

## Returns
* `(bf,dx,dy)`: Basis values and x/y derivatives.
"""
function basis_and_gradient(basis::RealPlaneWaves,i::Int,k::T,pts::AbstractArray) where {T<:Real}
    parx=basis.parity_x[i]
    pary=basis.parity_y[i]
    vx=cos(basis.angles[i])
    vy=sin(basis.angles[i])
    fx=_rpw_fun(parx)
    fy=_rpw_fun(pary)
    dfx=_drpw_fun(parx)
    dfy=_drpw_fun(pary)
    M=length(pts)
    bf=Vector{T}(undef,M)
    dx=Vector{T}(undef,M)
    dy=Vector{T}(undef,M)
    @inbounds @simd for j=1:M
        x=pts[j][1]
        y=pts[j][2]
        ax=k*vx*x
        ay=k*vy*y
        bx=fx(ax)
        by=fy(ay)
        bf[j]=bx*by
        dx[j]=k*vx*dfx(ax)*by
        dy[j]=bx*k*vy*dfy(ay)
    end
    return bf,dx,dy
end

"""
    basis_and_gradient(basis::RealPlaneWaves,indices::AbstractArray,k::T,pts::AbstractArray;multithreaded::Bool=true) where {T<:Real} → (B,dB_dx,dB_dy)

Evaluate selected basis functions and their spatial gradients in one pass.

## Arguments
* `basis::RealPlaneWaves`: Basis.
* `indices::AbstractArray`: Basis-function indices.
* `k::T`: Wavenumber.
* `pts::AbstractArray`: Evaluation points.

## Keyword Arguments
* `multithreaded::Bool=true`: Construct columns in parallel.

## Returns
* `(B,dB_dx,dB_dy)`: Basis matrix and x/y derivative matrices.
"""
function basis_and_gradient(basis::RealPlaneWaves,indices::AbstractArray,k::T,pts::AbstractArray;multithreaded::Bool=true) where {T<:Real}
    M=length(pts)
    N=length(indices)
    B=Matrix{T}(undef,M,N)
    dBdx=Matrix{T}(undef,M,N)
    dBdy=Matrix{T}(undef,M,N)
    @use_threads multithreading=multithreaded for c in 1:N
        idx=indices[c]
        parx=basis.parity_x[idx]
        pary=basis.parity_y[idx]
        vx=cos(basis.angles[idx])
        vy=sin(basis.angles[idx])
        fx=_rpw_fun(parx)
        fy=_rpw_fun(pary)
        dfx=_drpw_fun(parx)
        dfy=_drpw_fun(pary)
        col=@view B[:,c]
        cx=@view dBdx[:,c]
        cy=@view dBdy[:,c]
        @inbounds @simd for j=1:M
            x=pts[j][1]
            y=pts[j][2]
            ax=k*vx*x
            ay=k*vy*y
            bx=fx(ax)
            by=fy(ay)
            col[j]=bx*by
            cx[j]=k*vx*dfx(ax)*by
            cy[j]=bx*k*vy*dfy(ay)
        end
    end
    return B,dBdx,dBdy
end

"""
    dk_fun(basis::RealPlaneWaves,i::Int,k::T,pts::AbstractArray) where {T<:Real} → dk::Vector{T}

Evaluate the derivative with respect to `k` of the `i`-th basis function.

## Arguments
* `basis::RealPlaneWaves`: Basis.
* `i::Int`: Basis-function index.
* `k::T`: Wavenumber.
* `pts::AbstractArray`: Evaluation points.

## Returns
* `dk::Vector{T}`: Derivative with respect to `k`.
"""
@inline function dk_fun(basis::RealPlaneWaves,i::Int,k::T,pts::AbstractArray) where {T<:Real}
    parx=basis.parity_x[i]
    pary=basis.parity_y[i]
    vx=cos(basis.angles[i])
    vy=sin(basis.angles[i])
    fx=_rpw_fun(parx)
    fy=_rpw_fun(pary)
    dfx=_drpw_fun(parx)
    dfy=_drpw_fun(pary)
    M=length(pts)
    dk=Vector{T}(undef,M)
    @inbounds @simd for j=1:M
        x=pts[j][1]
        y=pts[j][2]
        ax=k*vx*x
        ay=k*vy*y
        bx=fx(ax)
        by=fy(ay)
        dk[j]=vx*x*dfx(ax)*by+bx*vy*y*dfy(ay)
    end
    return dk
end

"""
    dk_fun(basis::RealPlaneWaves,indices::AbstractArray,k::T,pts::AbstractArray;multithreaded::Bool=true) where {T<:Real} → dB_dk::Matrix{T}

Evaluate derivatives with respect to `k` of selected basis functions.

## Arguments
* `basis::RealPlaneWaves`: Basis.
* `indices::AbstractArray`: Basis-function indices.
* `k::T`: Wavenumber.
* `pts::AbstractArray`: Evaluation points.

## Keyword Arguments
* `multithreaded::Bool=true`: Construct columns in parallel.

## Returns
* `dB_dk::Matrix{T}`: Derivative matrix with respect to `k`.
"""
@inline function dk_fun(basis::RealPlaneWaves,indices::AbstractArray,k::T,pts::AbstractArray;multithreaded::Bool=true) where {T<:Real}
    M=length(pts)
    N=length(indices)
    dBdk=Matrix{T}(undef,M,N)
    @use_threads multithreading=multithreaded for c in 1:N
        idx=indices[c]
        parx=basis.parity_x[idx]
        pary=basis.parity_y[idx]
        vx=cos(basis.angles[idx])
        vy=sin(basis.angles[idx])
        fx=_rpw_fun(parx)
        fy=_rpw_fun(pary)
        dfx=_drpw_fun(parx)
        dfy=_drpw_fun(pary)
        col=@view dBdk[:,c]
        @inbounds @simd for j=1:M
            x=pts[j][1]
            y=pts[j][2]
            ax=k*vx*x
            ay=k*vy*y
            bx=fx(ax)
            by=fy(ay)
            col[j]=vx*x*dfx(ax)*by+bx*vy*y*dfy(ay)
        end
    end
    return dBdk
end
