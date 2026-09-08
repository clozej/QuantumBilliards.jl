"""
    pad_limits(xlim, ylim; padding::Real = 0.01) → (xlim_padded::Tuple, ylim_padded::Tuple)

Pads a pair of `(min, max)` limits `xlim` and `ylim` symmetrically by
`padding` on each side.

## Arguments
* `xlim`: The `(xmin, xmax)` limits to pad.
* `ylim`: The `(ymin, ymax)` limits to pad.

## Keyword arguments
*  `padding::Real = 0.01` : Amount subtracted from the lower limit and added to the upper limit of each pair.

## Returns
*  `xlim_padded` : `(xlim[1] - padding, xlim[2] + padding)`.
*  `ylim_padded` : `(ylim[1] - padding, ylim[2] + padding)`.
"""
function pad_limits(xlim, ylim; padding=0.01)
    return (xlim[1] - padding, xlim[2] + padding), (ylim[1] - padding, ylim[2] + padding)
end

"""
    rectify_grid(grid::AbstractVector) → new_grid::AbstractVector

Shifts `grid` so that its entry closest to zero is exactly `0`, then keeps
only the strictly positive half of the shifted grid; returns `grid`
unchanged if it does not straddle zero.

## Description
This is used to build a half-grid for wavefunction evaluation on the
fundamental domain of a reflection-symmetric billiard, where only the `x > 0`
(or `y > 0`) half-plane needs to be sampled.

## Arguments
* `grid`: The coordinate grid to rectify.

## Returns
*  `new_grid` : The shifted, strictly-positive half of `grid` if `grid` straddles zero (i.e. `grid[1] <= 0 <= grid[end]`), otherwise `grid` unchanged.
"""
function rectify_grid(grid)
    type = eltype(grid)
    if grid[1] <= zero(type) <= grid[end]
        idx = argmin(abs.(grid))
        new_grid = grid .- grid[idx] #.- ds/2.0
        return new_grid[new_grid .> zero(type)]
    else
        return grid
    end
end

"""
    boundary_limits(curves; grd::Int = 1000, padding::Real = 0.01) → (xlim::Tuple, ylim::Tuple)

Computes padded bounding-box limits `(xlim, ylim)` enclosing a collection of
boundary `curves`.

## Description
Each curve is sampled at `N_bnd = max(512, round(Int, grd/L))` equally spaced
parameter values (with `L` the curve length), the sampled points from all
curves are pooled, and the extrema of their `x` and `y` coordinates are
padded with [`pad_limits`](@ref).

## Arguments
* `curves`: A collection of boundary curves to sample, e.g. from `get_boundary_curves`.

## Keyword arguments
*  `grd::Int = 1000` : Target total sampling density (points per unit curve length) used to set the number of points sampled per curve.
*  `padding::Real = 0.01` : Padding added to the bounding box, passed to [`pad_limits`](@ref).

## Returns
*  `xlim` : Padded `(xmin, xmax)` limits enclosing all sampled boundary points.
*  `ylim` : Padded `(ymin, ymax)` limits enclosing all sampled boundary points.
"""
function boundary_limits(curves; grd=1000, padding=0.01) 
    x_bnd = Vector{Any}()
    y_bnd = Vector{Any}()
    for crv in curves #names of variables not very nice
        L = crv.length
        N_bnd = max(512,round(Int, grd/L))
        t = range(0.0,1.0, N_bnd)[1:end-1]
        pts = curve(crv,t)
        append!(x_bnd, getindex.(pts,1))
        append!(y_bnd, getindex.(pts,2))
    end
    x_bnd[end] = x_bnd[1]
    y_bnd[end] = y_bnd[1]
    xlim = extrema(x_bnd)
    #dx =  xlim[2] - xlim[1]
    ylim = extrema(y_bnd)
    #dy =  ylim[2] - ylim[1]
    return pad_limits(xlim, ylim; padding=padding)
end


"""
    ϕ_slp(x::T, y::T, k::T, pts::BoundaryPoints{T}, u::AbstractVector; float32_bessel::Bool = true, use_chebyshev::Bool = false) where {T<:Real} → ψ::Number

Evaluates the single-layer-potential (SLP) Green's-function reconstruction of
a Dirichlet eigenfunction at `(x,y)`,

```math
\\psi(x,y) = \\frac{1}{4}\\int_{\\partial\\Omega} Y_0(k|x-q|)\\,u(q)\\,ds_q,
```

where `u = ∂ₙψ` is the boundary normal derivative (e.g. from
[`boundary_function`](@ref)). The overall sign is irrelevant for an
eigenfunction. This is the reconstruction kernel behind
[`wavefunction(state::BIMEigenstate)`](@ref).

!!! note "Chebyshev acceleration not yet implemented"
    `use_chebyshev = true` is reserved for a future Chebyshev-interpolated
    `Y₀` evaluation (see the BIM-solver migration plan's Chebyshev
    acceleration step); passing it currently raises an error.

## Arguments
* `x`,`y`: Evaluation coordinates.
* `k`: Wavenumber.
* `pts`: Boundary discretization providing `pts.xy` and `pts.ds`.
* `u`: Boundary normal derivative values at `pts`.

## Keyword arguments
*  `float32_bessel::Bool = true` : Evaluate `Y₀` in `Float32` arithmetic for speed, converting back to `T`.
*  `use_chebyshev::Bool = false` : Reserved for future Chebyshev-interpolated `Y₀` evaluation; must currently be `false`.

## Returns
*  `ψ` : The reconstructed wavefunction value at `(x,y)`.
"""
@inline function ϕ_slp(x::T, y::T, k::T, pts::BoundaryPoints{T}, u::AbstractVector; float32_bessel::Bool=true, use_chebyshev::Bool=false) where {T<:Real}
    use_chebyshev && error("Chebyshev SLP wavefunction reconstruction not yet implemented; see migration plan step 10")
    xy = pts.xy
    ds = pts.ds
    S = eltype(u)
    acc = zero(S)
    @inbounds @fastmath for j in eachindex(u)
        p = xy[j]
        dx = x - p[1]
        dy = y - p[2]
        r2 = muladd(dx, dx, dy*dy)
        r2 == zero(T) && continue # only guard exact coincidence with a source node
        r = sqrt(r2)
        y0 = float32_bessel ? T(Bessels.bessely0(Float32(k*r))) : Bessels.bessely0(k*r)
        acc += (y0*ds[j])*u[j]
    end
    return acc*T(0.25)
end

"""
    compute_psi(state::S, x_grid::AbstractVector, y_grid::AbstractVector; inside_only::Bool = true, memory_limit::Real = 10.0e9, multithreaded::Bool = true) where {S<:AbsState} → Psi::Vector

Evaluates the wavefunction of `state` on the Cartesian grid formed by
`x_grid` and `y_grid`, returning it as a flattened vector.

## Description
The evaluation points are the tensor-product grid `(x,y)` for `y in y_grid`,
`x in x_grid`, optionally restricted to points inside `state.billiard` when
`inside_only = true` (via `is_inside`). If the estimated memory required to
build the full basis matrix (`sizeof(eltype(vec)) * basis.dim * n_pts`) is
below `memory_limit`, the basis matrix is built in one shot with
`basis_matrix` and multiplied by the coefficient vector; otherwise the
wavefunction is accumulated basis function by basis function with
`basis_fun`, skipping coefficients smaller than `state.eps` in magnitude,
trading memory for compute time. Points outside the billiard (when
`inside_only = true`) are set to `NaN`.

## Arguments
* `state`: The eigenstate (`S<:AbsState`) whose wavefunction is evaluated.
* `x_grid`: Grid of `x` coordinates.
* `y_grid`: Grid of `y` coordinates.

## Keyword arguments
*  `inside_only::Bool = true` : Whether to evaluate only at points inside `state.billiard`, setting the wavefunction to `NaN` elsewhere.
*  `memory_limit::Real = 10.0e9` : Memory threshold (in bytes) above which the wavefunction is accumulated basis function by basis function instead of via a full basis matrix.
*  `multithreaded::Bool = true` : Whether the basis matrix construction is multithreaded.

## Returns
*  `Psi` : The wavefunction values on the flattened grid `(x_grid, y_grid)`, ordered as `x` varying fastest.
"""
function compute_psi(state::S, x_grid, y_grid; inside_only=true, memory_limit = 10.0e9, multithreaded = true) where {S<:AbsState}
    let vec = state.vec, k = state.k_basis, basis=state.basis, billiard=state.billiard, eps=state.eps #basis is correct size
        sz = length(x_grid)*length(y_grid)
        pts = collect(SVector(x,y) for y in y_grid for x in x_grid)
        if inside_only
            pts_mask = is_inside(billiard,pts)
            pts = pts[pts_mask]
        end
        n_pts = length(pts)
        #estimate max memory needed for the matrices
        type = eltype(vec)
        memory = sizeof(type)*basis.dim*n_pts
        Psi = zeros(type,sz)

        if memory < memory_limit
            B = basis_matrix(basis, k, pts; multithreaded)
            Psi_pts = B*vec
            if inside_only
                Psi[pts_mask] .= Psi_pts
            else
                Psi .= Psi_pts
            end

        else
            println("Warning: memory limit of $(Base.format_bytes(memory_limit)) exceded $(Base.format_bytes(memory)).")
            if inside_only
                for i in eachindex(vec)
                    if abs(vec[i]) > eps 
                        Psi[pts_mask] .+= vec[i].*basis_fun(basis,i,k,pts)
                    end
                end
            else
                for i in eachindex(vec)
                    if abs(vec[i]) > eps 
                        Psi .+= vec[i].*basis_fun(basis,i,k,pts)
                    end
                end
            end
        end
        if inside_only
            Psi[.!pts_mask] .= convert(type, NaN)
        end
        return Psi
    end
end

"""
    wavefunction(state::S; b::Real = 5.0, inside_only::Bool = true, fundamental_domain::Bool = true, memory_limit::Real = 10.0e9, multithreaded::Bool = true) where {S<:AbsState} → (Psi2d::Matrix, x_grid::Vector, y_grid::Vector)

Computes the wavefunction of an eigenstate `state` on a regular grid covering
its billiard, optionally unfolding it from the fundamental domain onto the
full billiard.

## Description
A bounding box for `state.billiard` is computed with [`boundary_limits`](@ref)
at a sampling density of `max(1000, round(Int, k*L*b/(2*pi)))` (with `L` the
boundary length), and grids `x_grid`, `y_grid` are built with
`max(round(Int, k*d*b/(2*pi)), 512)` points along each dimension `d` (`dx` or
`dy`), giving roughly `b` grid points per de Broglie wavelength. If the basis
carries reflection symmetries, the corresponding grid(s) are restricted to
the fundamental domain with [`rectify_grid`](@ref). The wavefunction is then
evaluated with [`compute_psi`](@ref) and reshaped into a 2D array `Psi2d`. If
`fundamental_domain = false` and the basis has symmetries, `Psi2d` and the
grids are unfolded onto the full billiard with
[`apply_symmetries_to_wavefunction`](@ref).

## Arguments
* `state`: The eigenstate (`S<:AbsState`) for which the wavefunction is computed.

## Keyword arguments
*  `b::Real = 5.0` : Oversampling factor controlling the grid resolution; roughly `b` grid points per de Broglie wavelength.
*  `inside_only::Bool = true` : Whether to evaluate only at points inside `state.billiard`, passed to [`compute_psi`](@ref).
*  `fundamental_domain::Bool = true` : Whether to return the wavefunction restricted to the symmetry-reduced fundamental domain (`true`) or unfolded onto the full billiard (`false`).
*  `memory_limit::Real = 10.0e9` : Memory threshold (in bytes) passed to [`compute_psi`](@ref).
*  `multithreaded::Bool = true` : Whether the underlying matrix construction is multithreaded.

## Returns
*  `Psi2d` : The wavefunction values on the grid `(x_grid, y_grid)`.
*  `x_grid` : The `x` coordinates of the grid.
*  `y_grid` : The `y` coordinates of the grid.
"""
function wavefunction(state::S; b=5.0, inside_only=true, fundamental_domain = true, memory_limit = 10.0e9, multithreaded = true) where {S<:AbsState}
    let k = state.k, billiard=state.billiard, symmetries=state.basis.symmetries     
        #println(new_basis.dim)
        type = eltype(state.vec)
        #try to find a lazy way to do this
        L = CompositeCurve(get_boundary_curves(billiard)).length
        
        xlim,ylim = boundary_limits(get_boundary_curves(billiard); grd=max(1000,round(Int, k*L*b/(2*pi))))
        dx = xlim[2] - xlim[1]
        dy = ylim[2] - ylim[1]
        nx = max(round(Int, k*dx*b/(2*pi)), 512)
        ny = max(round(Int, k*dy*b/(2*pi)), 512)
        x_grid::Vector{type} = collect(type,range(xlim... , nx))
        y_grid::Vector{type} = collect(type,range(ylim... , ny))

        if ~isnothing(symmetries)
            has_x = any(s -> s isa BilliardGeometry.XAxisReflection, symmetries)
            has_y = any(s -> s isa BilliardGeometry.YAxisReflection, symmetries)
            if has_x
                x_grid = rectify_grid(x_grid)
                nx = length(x_grid)
            end
            if has_y
                y_grid = rectify_grid(y_grid)
                ny = length(y_grid)
            end
        end

        Psi::Vector{type} = compute_psi(state,x_grid,y_grid; inside_only, memory_limit, multithreaded) 
        #println("Psi type $(eltype(Psi)), $(memory_size(Psi))")
        Psi2d::Array{type,2} = reshape(Psi, (nx,ny))
        if ~fundamental_domain 
            if ~isnothing(symmetries)
                Psi2d, x_grid, y_grid = apply_symmetries_to_wavefunction(Psi2d,x_grid,y_grid,symmetries,state.basis.sym_qnumbers)
            end
        end
        return Psi2d, x_grid, y_grid
    end
end

"""
    wavefunction(state::BIMEigenstate{K,T,S,Bi}; b::Union{Real,Symbol} = :auto, inside_only::Bool = true, use_float32_bessel::Bool = true, multithreaded::Bool = true) where {K,T,S<:DoubleLayerPotentialSolver,Bi} → (Psi2d::Matrix, x_grid::Vector, y_grid::Vector)

Reconstructs the wavefunction of a [`BIMEigenstate`](@ref) computed with a
[`DoubleLayerPotentialSolver`](@ref) on a Cartesian grid, via the
single-layer-potential Green's-function integral [`ϕ_slp`](@ref) applied to
the boundary normal derivative `u = ∂ₙψ` from [`_dlp_boundary_function`](@ref).

## Description
The Cartesian grid covers the *complete* physical boundary
(`full_boundary(billiard)` when `solver.symmetry !== nothing`, else
`get_boundary_curves(billiard)`), so a symmetry-reduced solver still
reconstructs the wavefunction on the whole billiard, never only its
fundamental domain. Grid points inside the billiard (when `inside_only =
true`) are evaluated with [`ϕ_slp`](@ref); the outer loop over masked grid
points is parallelized (`multithreaded`), with each thread writing to a
disjoint output index.

## Keyword arguments
*  `b::Union{Real,Symbol} = :auto` : Grid sampling density in points per de Broglie wavelength; `:auto` uses `solver.pts_scaling_factor[1]`.
*  `inside_only::Bool = true` : Whether to evaluate only at points inside `state.billiard`.
*  `use_float32_bessel::Bool = true` : Passed to [`ϕ_slp`](@ref).
*  `multithreaded::Bool = true` : Whether the boundary-function/grid evaluation is threaded.

## Returns
*  `Psi2d` : The reconstructed wavefunction values on the grid.
*  `x_grid`,`y_grid` : The Cartesian grid coordinates.
"""
function wavefunction(state::BIMEigenstate{K,T,S,Bi}; b::Union{Real,Symbol}=:auto, inside_only::Bool=true, use_float32_bessel::Bool=true, multithreaded::Bool=true) where {K,T,S<:DoubleLayerPotentialSolver,Bi}
    solver = state.solver
    billiard = state.billiard
    kT = real(state.k)
    u, pts, _ = _dlp_boundary_function(solver, billiard, kT; multithreaded)
    comp = solver.symmetry === nothing ? get_boundary_curves(billiard) : full_boundary(billiard)
    bval = b === :auto ? solver.pts_scaling_factor[1] : T(b)
    Ltot = sum(crv.length for crv in comp)
    xlim, ylim = boundary_limits(comp; grd=max(1000, round(Int, kT*Ltot*bval/(2*pi))))
    dx = xlim[2]-xlim[1]
    dy = ylim[2]-ylim[1]
    nx = max(round(Int, kT*dx*bval/(2*pi)), 512)
    ny = max(round(Int, kT*dy*bval/(2*pi)), 512)
    x_grid::Vector{T} = collect(T, range(xlim..., nx))
    y_grid::Vector{T} = collect(T, range(ylim..., ny))
    pts_grid = collect(SVector(x,y) for y in y_grid for x in x_grid)
    pts_mask = inside_only ? is_inside(billiard, pts_grid) : trues(length(pts_grid))
    Stype = eltype(u) <: Real ? T : Complex{T}
    Psi = zeros(Stype, nx*ny)
    idxs = findall(pts_mask)
    @use_threads multithreading=multithreaded for jj in eachindex(idxs)
        idx = idxs[jj]
        p = pts_grid[idx]
        Psi[idx] = ϕ_slp(p[1], p[2], kT, pts, u; float32_bessel=use_float32_bessel)
    end
    inside_only && (Psi[.!pts_mask] .= convert(Stype, NaN))
    Psi2d::Matrix{Stype} = reshape(Psi, (nx, ny))
    return Psi2d, x_grid, y_grid
end

"""
    wavefunction(state::BasisState; xlim::Tuple = (-2.0, 2.0), ylim::Tuple = (-2.0, 2.0), b::Real = 5.0) → (Psi2d::Matrix, x_grid::Vector, y_grid::Vector)

Computes a single basis function represented by `state` on a regular grid
over the fixed box `xlim × ylim`, without reference to any billiard geometry.

## Description
Grids `x_grid`, `y_grid` are built with `max(round(Int, k*d*b/(2*pi)), 512)`
points along each dimension `d` (`dx` or `dy`, from `xlim`/`ylim`), giving
roughly `b` grid points per de Broglie wavelength. The `idx`-th basis function
of `state.basis` is evaluated directly on the grid with `basis_fun` and
reshaped into a 2D array.

## Arguments
* `state`: The [`BasisState`](@ref) whose basis function is evaluated.

## Keyword arguments
*  `xlim::Tuple = (-2.0, 2.0)` : The `(xmin, xmax)` extent of the evaluation grid.
*  `ylim::Tuple = (-2.0, 2.0)` : The `(ymin, ymax)` extent of the evaluation grid.
*  `b::Real = 5.0` : Oversampling factor controlling the grid resolution; roughly `b` grid points per de Broglie wavelength.

## Returns
*  `Psi2d` : The basis function values on the grid `(x_grid, y_grid)`.
*  `x_grid` : The `x` coordinates of the grid.
*  `y_grid` : The `y` coordinates of the grid.
"""
function wavefunction(state::BasisState; xlim =(-2.0,2.0), ylim=(-2.0,2.0), b=5.0) 
    let k = state.k, basis=state.basis      
        #println(new_basis.dim)
        type = eltype(state.vec)
        #try to find a lazy way to do this
        dx = xlim[2] - xlim[1]
        dy = ylim[2] - ylim[1]
        nx = max(round(Int, k*dx*b/(2*pi)), 512)
        ny = max(round(Int, k*dy*b/(2*pi)), 512)
        x_grid::Vector{type} = collect(type,range(xlim... , nx))
        y_grid::Vector{type} = collect(type,range(ylim... , ny))
        pts_grid = [SVector(x,y) for y in y_grid for x in x_grid]
        Psi::Vector{type} = basis_fun(basis,state.idx,k,pts_grid) 
        #println("Psi type $(eltype(Psi)), $(memory_size(Psi))")
        Psi2d::Array{type,2} = reshape(Psi, (nx,ny))
        return Psi2d, x_grid, y_grid
    end
end


