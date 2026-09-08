
"""
    antisym_vec(x::AbstractVector) → v::Vector

Antisymmetrically extends a one-sided vector `x` (assumed sorted with
`x[1]` closest to zero) to negative values, returning `[-reverse(x[2:end]);
x]`.

## Description
This is used to build symmetric grids of coherent-state evaluation points (in
`q` or `p`) around zero from a one-sided grid `x`, by mirroring and negating
all but the first entry and prepending the result to `x`.

## Arguments
* `x`: The one-sided vector to extend, typically non-negative and increasing.

## Returns
*  `v` : The antisymmetrically extended vector `[-reverse(x[2:end]); x]`, of length `2*length(x) - 1`.
"""
function antisym_vec(x)
    v = reverse(-x[2:end])
    return append!(v,x)
end

"""
    husimi_function(k, u::AbstractVector, s::AbstractVector, L::Real; c::Real = 10.0, w::Real = 7.0) → (H::Matrix, qs::Vector, ps::Vector)

Computes the boundary Husimi function of the normal-derivative boundary
function `u(s)` (sampled at equally spaced arc-length points `s`) at
wavenumber `k`, on a boundary of total length `L`, returning the Husimi
density `H` on a grid of arc-length coordinates `qs` and momenta `ps`.

## Description
`u` is projected onto (approximately) minimal-uncertainty coherent-state
wavepackets of Gaussian width `sig = 1/sqrt(k)` in the arc-length coordinate,
truncated to `w` widths (`x = s[s .<= w*sig]`) and evaluated with `c` points
per width along the momentum direction. Because the boundary is periodic with
period `L`, each coherent state is periodized by adding its two nearest
periodic images at `s ± L` (`gauss_l`, `gauss_r`) before contraction with
`u`, giving the overlap

```math
h(q,p) = \\sum_{s} u(s)\\, \\big[g(s-q) + g(s-q+L) + g(s-q-L)\\big]\\, e^{i k p (s-q)},
```

with `g(\\cdot) = \\exp(-k\\, (\\cdot)^2/2)` the Gaussian envelope. The Husimi
density is `H(q,p) = a\\, |h(q,p)|^2`, with normalization constant
`a = 1/(2\\pi\\sqrt{\\pi k})` (not normalized so that `H` integrates to `1`).
Momenta are sampled on `ps ∈ [0,1]` in steps of `sig/c` and then mirrored to
`[-1,1]` via [`antisym_vec`](@ref), and arc-length points `qs` are subsampled
from `s` with the same step `sig/c`.

## Arguments
* `k`: The wavenumber of the eigenstate.
* `u`: The (real) boundary normal-derivative function, typically from [`boundary_function`](@ref).
* `s`: The arc-length coordinates of `u`, assumed equally spaced.
* `L`: The total length of the (periodized) boundary.

## Keyword arguments
*  `c::Real = 10.0` : Number of coherent-state evaluation points per Gaussian width `sig`, controlling the resolution of `qs` and `ps`.
*  `w::Real = 7.0` : Truncation width of the coherent-state Gaussian envelope, in units of `sig`.

## Returns
*  `H` : The Husimi density on the grid `(qs, ps)`.
*  `qs` : Arc-length coordinates at which `H` is sampled.
*  `ps` : Momenta (in units of `k`, ranging over `[-1,1]`) at which `H` is sampled.
"""
function husimi_function(k,u,s,L; c = 10.0, w = 7.0)
    #c density of points in coherent state peak, w width in units of sigma
    #L is the boundary length for periodization
    #compute coherrent state weights
    N = length(s)
    sig = one(k)/sqrt(k) #width of the gaussian
    x = s[s.<=w*sig]
    idx = length(x) #do not change order here
    x = antisym_vec(x)
    a = one(k)/(2*pi*sqrt(pi*k)) #normalization factor in this version Hsimi is not noramlized to 1
    ds = (x[end]-x[1])/length(x) #integration weigth
    uc = CircularVector(u) #allows circular indexing
    gauss = @. exp(-k/2*x^2)*ds
    gauss_l = @. exp(-k/2*(x+L)^2)*ds
    gauss_r = @. exp(-k/2*(x-L)^2)*ds
    #construct evaluation points in p coordinate
    ps = collect(range(0.0,1.0,step = sig/c))
    #construct evaluation points in q coordinate
    q_stride = length(s[s.<=sig/c])
    q_idx = collect(1:q_stride:N)
    push!(q_idx,N) #add last point
    qs = s[q_idx]
    #println(length(qs))
    H = zeros(typeof(k),length(qs),length(ps))
    for i in eachindex(ps)   
        cs = @. exp(im*ps[i]*k*x)*gauss + exp(im*ps[i]*k*(x+L))*gauss_l + exp(im*ps[i]*k*(x-L))*gauss_r#imag part of coherent state
        for j in eachindex(q_idx)
            u_w = uc[q_idx[j]-idx+1:q_idx[j]+idx-1] #window with relevant values of u
            h = sum(cs.*u_w)
            #hi = sum(ci.*u_w)
            H[j,i] = a*abs2(h)
        end
    end

    ps = antisym_vec(ps)
    H_ref = reverse(H[:, 2:end]; dims=2)
    H = hcat(H_ref,H)
     
    return H, qs, ps    
end

"""
    husimi_function(state::S; b::Real = 5.0, c::Real = 10.0, w::Real = 7.0) where {S<:AbsState} → (H::Matrix, qs::Vector, ps::Vector)

Computes the boundary Husimi function of an eigenstate directly from `state`,
by combining [`boundary_function`](@ref) and
[`husimi_function(k, u, s, L)`](@ref husimi_function).

## Arguments
* `state`: The eigenstate for which the boundary Husimi function is computed.

## Keyword arguments
*  `b::Real = 5.0` : Oversampling factor passed to [`boundary_function`](@ref) controlling the boundary point density.
*  `c::Real = 10.0` : Number of coherent-state evaluation points per Gaussian width, passed to [`husimi_function(k, u, s, L)`](@ref husimi_function).
*  `w::Real = 7.0` : Truncation width of the coherent-state Gaussian envelope, passed to [`husimi_function(k, u, s, L)`](@ref husimi_function).

## Returns
*  `H` : The Husimi density on the grid `(qs, ps)`.
*  `qs` : Arc-length coordinates at which `H` is sampled.
*  `ps` : Momenta at which `H` is sampled.
"""
function husimi_function(state::S;  b = 5.0, c = 10.0, w = 7.0) where {S<:AbsState}
    L = sum(crv.length for crv in full_boundary(state.billiard))
    k = state.k
    u, s, norm = boundary_function(state; b=b)
    return husimi_function(k,u,s,L; c = c, w = w)
end

# Checks (to within `rtol`) whether `s` is uniformly spaced and `ds` is
# uniform, i.e. whether the fast sliding-window `husimi_function(k,u,s,L)`
# kernel's assumption of a uniform arc-length grid actually holds. A BIM
# solver using `GlobalCornerGrading` intentionally clusters nodes near
# corners, so `s` is not uniformly spaced there.
@inline function _bim_boundary_uniformly_spaced(s::AbstractVector{T}, ds::AbstractVector{T}; rtol::Real=1e-6) where {T<:Real}
    N = length(s)
    N < 3 && return true
    Δs = s[2]-s[1]
    return all(isapprox(s[i+1]-s[i], Δs; rtol=rtol) for i in 2:N-1) && all(isapprox(ds[i], ds[1]; rtol=rtol) for i in 2:N)
end

"""
    husimi_function(state::BIMEigenstate{K,T,S,Bi}; b::Real = 5.0, c::Real = 10.0, w::Real = 7.0, multithreaded::Bool = true) where {K,T,S<:DoubleLayerPotentialSolver,Bi} → (H::Matrix, qs::Vector, ps::Vector)

Computes the boundary Husimi function of a [`BIMEigenstate`](@ref) computed
with a [`DoubleLayerPotentialSolver`](@ref), by combining
[`_dlp_boundary_function`](@ref) and
[`husimi_function(k, u, s, L)`](@ref husimi_function).

!!! warning "Non-uniform boundary sampling"
    [`husimi_function(k, u, s, L)`](@ref husimi_function) assumes `s` is
    uniformly spaced in arc-length. A [`DoubleLayerPotentialSolver`](@ref)
    using [`GlobalCornerGrading`](@ref) intentionally clusters nodes near
    corners, so `s` is *not* uniformly spaced there, and a warning is emitted
    in that case since the resulting Husimi function may be inaccurate.

## Keyword arguments
*  `b::Real = 5.0` : Unused, accepted for interface compatibility (see [`boundary_function(state::BIMEigenstate)`](@ref)).
*  `c::Real = 10.0`, `w::Real = 7.0` : Passed to [`husimi_function(k, u, s, L)`](@ref husimi_function).
*  `multithreaded::Bool = true` : Whether the adjoint Fredholm matrix assembly is multithreaded.
"""
function husimi_function(state::BIMEigenstate{K,T,S,Bi}; b=5.0, c=10.0, w=7.0, multithreaded=true) where {K,T,S<:DoubleLayerPotentialSolver,Bi}
    u, pts, norm = _dlp_boundary_function(state.solver, state.billiard, state.k; multithreaded)
    comp = state.solver.symmetry === nothing ? get_boundary_curves(state.billiard) : full_boundary(state.billiard)
    L = sum(crv.length for crv in comp)
    _bim_boundary_uniformly_spaced(pts.s, pts.ds) || @warn "BIMEigenstate boundary discretization is not uniformly spaced in arclength (e.g. GlobalCornerGrading clusters nodes near corners); husimi_function assumes uniform spacing and may be inaccurate."
    return husimi_function(real(state.k), u, pts.s, L; c=c, w=w)
end

