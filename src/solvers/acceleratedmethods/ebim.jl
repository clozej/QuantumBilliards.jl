"""
    ExpandedBIMSolver{T,K} <: AcceleratedBIMSolver

`ExpandedBIMSolver` is a concrete [`AcceleratedBIMSolver`](@ref) implementing
the expanded boundary integral method (EBIM): a single locally-corrected root
of the nonlinear eigenproblem `A(k)v = 0` obtained from a second-order local
Taylor expansion of `A(k)`.

## Description
`ExpandedBIMSolver` wraps an inner `kernel::SweepBIMSolver` (a
[`DoubleLayerPotentialSolver`](@ref), [`CombinedFieldIntegralEquationSolver`](@ref)
or [`CompositeBIMSolver`](@ref)) supplying the Fredholm operator and its first
two `k`-derivatives,

    A(k+ε) = A(k) + ε A'(k) + (1/2) ε² A''(k) + O(ε³).

The generalized eigenproblem `A(k)v = λ A'(k)v` gives the first-order root
correction `ε₁ = -λ`; with the corresponding left generalized eigenvector `u`,
the second-order correction is `ε₂ = -(1/2) ε₁² [u'A''(k)v]/[u'A'(k)v]`, giving
the corrected wavenumber `k_corr = k + ε₁ + ε₂` (see [`construct_matrices`](@ref),
[`solve`](@ref)). Unlike [`BeynSolver`](@ref), EBIM produces a single locally
refined root per call rather than every root in a window.

## Attributes
* `kernel`: The wrapped [`SweepBIMSolver`](@ref) supplying `A(k)`, `A'(k)`, `A''(k)`.
* `use_chebyshev`: Whether Chebyshev-accelerated kernel evaluation is used.
* `n_panels_h`: Hankel-function Chebyshev panel count.
* `M_h`: Hankel-function Chebyshev polynomial degree.
* `n_panels_j`: Bessel-J-function Chebyshev panel count.
* `M_j`: Bessel-J-function Chebyshev polynomial degree.

## API
The following functions can be evaluated for this type:
- [`evaluate_points`](@ref)
- [`construct_matrices`](@ref)
- [`solve`](@ref)
- [`solve_wavenumber`](@ref)
- [`solve_spectrum`](@ref)

!!! note "Migration status"
    `use_chebyshev` is not yet wired to an accelerated evaluation path (Step
    10 of the migration plan); `construct_matrices` always uses direct
    Bessels.jl/SpecialFunctions.jl evaluation via the wrapped kernel
    regardless of `solver.use_chebyshev`. Supports
    [`DoubleLayerPotentialSolver`](@ref), [`CombinedFieldIntegralEquationSolver`](@ref)
    and [`CompositeBIMSolver`](@ref) kernels.
"""
struct ExpandedBIMSolver{T<:Real,K<:SweepBIMSolver} <: AcceleratedBIMSolver
    kernel::K
    use_chebyshev::Bool
    n_panels_h::Int
    M_h::Int
    n_panels_j::Int
    M_j::Int
end

"""
    ExpandedBIMSolver(kernel::K; use_chebyshev::Bool = true, n_panels_h::Int = 15000, M_h::Int = 5, n_panels_j::Int = 10000, M_j::Int = 5) where {K<:SweepBIMSolver} → solver::ExpandedBIMSolver

Constructs an [`ExpandedBIMSolver`](@ref) wrapping the boundary-integral
`kernel`.

## Arguments
* `kernel`: The [`SweepBIMSolver`](@ref) supplying the Fredholm operator `A(k)` and its `k`-derivatives.

## Keyword arguments
* `use_chebyshev::Bool = true`: Whether to use Chebyshev-accelerated kernel evaluation.
* `n_panels_h::Int = 15000`: Hankel-function Chebyshev panel count.
* `M_h::Int = 5`: Hankel-function Chebyshev polynomial degree.
* `n_panels_j::Int = 10000`: Bessel-J-function Chebyshev panel count.
* `M_j::Int = 5`: Bessel-J-function Chebyshev polynomial degree.

## Returns
* `solver`: An [`ExpandedBIMSolver`](@ref) instance.
"""
function ExpandedBIMSolver(kernel::K; use_chebyshev::Bool=true,
                            n_panels_h::Int=15000, M_h::Int=5, n_panels_j::Int=10000, M_j::Int=5) where {K<:SweepBIMSolver}
    T = _bim_numeric_type(kernel)
    return ExpandedBIMSolver{T,K}(kernel, use_chebyshev, n_panels_h, M_h, n_panels_j, M_j)
end

_bim_numeric_type(::ExpandedBIMSolver{T}) where {T} = T

################################################################################
###################### DERIVATIVE-OF-HANKEL-KERNEL HELPERS ###################
################################################################################

# Wavenumber-derivatives of a "linear-in-k prefactor" order-1 kernel term
# value(k) = α(k)*invr*Z1(kr), with α(k) = c*k linear in k (as every DLP/CFIE
# double-layer prefactor αL1=-k/(2π), αL2=ik/2 is) and Z the Bessel-J or
# Hankel-H1 family (Zν satisfies Z1'(z)=Z0(z)-Z1(z)/z). Using α(k)=c*k and
# dα/dk=α/k, the chain rule collapses to
#
#   dvalue/dk  = α*Z0(kr),
#   d²value/dk² = (α/k)*Z0(kr) - α*r*Z1(kr),
#
# with no leftover invr/r factors (they cancel exactly against the r
# introduced by d/dk[Z1(kr)]=r*Z1'(kr)). `z0`/`z1` are the already-evaluated
# `Z0(kr)`/`Z1(kr)`.
@inline function _ebim_lin1_deriv(α, r::T, invr::T, z0, z1, k) where {T<:Real}
    val = α*invr*z1
    dval = α*z0
    ddval = (α/k)*z0-α*r*z1
    return val, dval, ddval
end

# Wavenumber-derivatives of a "k-independent prefactor" order-0 kernel term
# value(k) = β*Z0(kr), with β constant in k (as the CFIE single-layer
# prefactors αM1=-1/(2π), αM2=i/2 are) and Z0'(z)=-Z1(z), Z1'(z)=Z0(z)-Z1(z)/z:
#
#   dvalue/dk   = -β*r*Z1(kr),
#   d²value/dk² = -β*r²*Z0(kr)+(β*r/k)*Z1(kr).
@inline function _ebim_const0_deriv(β, r::T, k, z0, z1) where {T<:Real}
    val = β*z0
    dval = -β*r*z1
    ddval = -β*r*r*z0+(β*r/k)*z1
    return val, dval, ddval
end

################################################################################
######################## DLP KERNEL WITH DERIVATIVES ##########################
################################################################################

# `D(k)` kernel entry (as `_dlp_kernel_entry` in dlp.jl) plus its first two
# `k`-derivatives, at full-boundary indices `(i,j)`.
@inline function _dlp_kernel_entry_with_derivatives(pts::BoundaryPoints{T}, Rmat::AbstractMatrix{T}, G::BoundaryGeomCache{T}, k::Union{T,Complex{T}}, i::Int, j::Int) where {T<:Real}
    if i==j
        return Complex{T}(pts.ws[i]*G.kappa[i], zero(T)), zero(Complex{T}), zero(Complex{T})
    end
    invtwopi = inv(2*T(pi))
    r = G.R[i,j]
    invr = G.invR[i,j]
    lt = G.logterm[i,j]
    inn = G.inner[i,j]
    h1 = _bim_hankelh1(1, k*r)
    j1 = _bim_besselj(1, k*r, h1)
    h0 = _bim_hankelh1(0, k*r)
    j0 = _bim_besselj(0, k*r, h0)
    αL1 = -k*invtwopi
    αL2 = im*k/2
    l1, dl1, ddl1 = _ebim_lin1_deriv(αL1*inn, r, invr, j0, j1, k)
    l2a, dl2a, ddl2a = _ebim_lin1_deriv(αL2*inn, r, invr, h0, h1, k)
    l2 = l2a-l1*lt
    dl2 = dl2a-dl1*lt
    ddl2 = ddl2a-ddl1*lt
    val = Rmat[i,j]*l1+pts.ws[j]*l2
    dval = Rmat[i,j]*dl1+pts.ws[j]*dl2
    ddval = Rmat[i,j]*ddl1+pts.ws[j]*ddl2
    return val, dval, ddval
end

################################################################################
######################## CFIE KERNEL WITH DERIVATIVES #########################
################################################################################

# `(D(k)+ikS(k))` kernel entry (as `_cfie_kernel_entry` in cfie.jl) plus its
# first two `k`-derivatives, at full-boundary indices `(i,j)`. The diagonal
# self-term's `S(k)` piece carries a `log(k²/4·speed²)` singular-kernel
# correction whose `k`-derivatives are `d/dk[log(k²/4·s²)]=2/k`; every
# off-diagonal `D`/`S` piece reuses `_ebim_lin1_deriv`/`_ebim_const0_deriv`.
@inline function _cfie_kernel_entry_with_derivatives(pts::BoundaryPoints{T}, Rmat::AbstractMatrix{T}, G::BoundaryGeomCache{T}, k::Union{T,Complex{T}}, i::Int, j::Int) where {T<:Real}
    invtwopi = inv(2*T(pi))
    ik = im*k
    if i==j
        si = G.speed[i]
        wi = pts.ws[i]
        Dval = Complex{T}(wi*G.kappa[i], zero(T))
        euler_over_pi = T(Base.MathConstants.eulergamma)/T(pi)
        m1 = -invtwopi*si
        m2 = ((Complex{T}(0, one(T)/2)-euler_over_pi)-invtwopi*log((k^2/4)*si^2))*si
        Sval = Complex{T}(Rmat[i,i]*m1, zero(T))+wi*m2
        val = Dval+ik*Sval
        dm2 = -si/(T(pi)*k)
        ddm2 = si/(T(pi)*k^2)
        dSval = wi*dm2
        ddSval = wi*ddm2
        dval = im*Sval+ik*dSval
        ddval = 2*im*dSval+ik*ddSval
        return val, dval, ddval
    end
    r = G.R[i,j]
    invr = G.invR[i,j]
    lt = G.logterm[i,j]
    inn = G.inner[i,j]
    sj = G.speed[j]
    wj = pts.ws[j]
    h0 = _bim_hankelh1(0, k*r)
    h1 = _bim_hankelh1(1, k*r)
    j0 = _bim_besselj(0, k*r, h0)
    j1 = _bim_besselj(1, k*r, h1)
    αL1 = -k*invtwopi
    αL2 = im*k/2
    αM1 = -invtwopi
    αM2 = Complex{T}(0, one(T)/2)
    l1, dl1, ddl1 = _ebim_lin1_deriv(αL1*inn, r, invr, j0, j1, k)
    l2a, dl2a, ddl2a = _ebim_lin1_deriv(αL2*inn, r, invr, h0, h1, k)
    l2 = l2a-l1*lt
    dl2 = dl2a-dl1*lt
    ddl2 = ddl2a-ddl1*lt
    Dval = Rmat[i,j]*l1+wj*l2
    dDval = Rmat[i,j]*dl1+wj*dl2
    ddDval = Rmat[i,j]*ddl1+wj*ddl2
    m1, dm1, ddm1 = _ebim_const0_deriv(αM1*sj, r, k, j0, j1)
    m2a, dm2a, ddm2a = _ebim_const0_deriv(αM2*sj, r, k, h0, h1)
    m2 = m2a-m1*lt
    dm2 = dm2a-dm1*lt
    ddm2 = ddm2a-ddm1*lt
    Sval = Rmat[i,j]*m1+wj*m2
    dSval = Rmat[i,j]*dm1+wj*dm2
    ddSval = Rmat[i,j]*ddm1+wj*ddm2
    val = Dval+ik*Sval
    dval = dDval+im*Sval+ik*dSval
    ddval = ddDval+2*im*dSval+ik*ddSval
    return val, dval, ddval
end

################################################################################
####################### COMPOSITE KERNEL WITH DERIVATIVES #####################
################################################################################

@inline _composite_component_kernel_entry_with_derivatives(::DoubleLayerPotentialSolver, pts::BoundaryPoints{T}, Rmat::AbstractMatrix{T}, G::BoundaryGeomCache{T}, k::Union{T,Complex{T}}, i::Int, j::Int) where {T<:Real} = _dlp_kernel_entry_with_derivatives(pts, Rmat, G, k, i, j)
@inline _composite_component_kernel_entry_with_derivatives(::CombinedFieldIntegralEquationSolver, pts::BoundaryPoints{T}, Rmat::AbstractMatrix{T}, G::BoundaryGeomCache{T}, k::Union{T,Complex{T}}, i::Int, j::Int) where {T<:Real} = _cfie_kernel_entry_with_derivatives(pts, Rmat, G, k, i, j)

# Smooth cross-component double-layer kernel entry (as
# `_composite_cross_kernel_entry` in compositebim.jl) plus its first two
# `k`-derivatives, for a `DoubleLayerPotentialSolver` source component.
@inline function _composite_cross_kernel_entry_with_derivatives(::DoubleLayerPotentialSolver, pb::BoundaryPoints{T}, xi::T, yi::T, k::Union{T,Complex{T}}, j::Int) where {T<:Real}
    xj, yj = pb.xy[j]
    dx = xi-xj
    dy = yi-yj
    r = hypot(dx, dy)
    invr = inv(r)
    tx, ty = pb.tangent[j]
    inn = ty*dx-tx*dy
    h1 = _bim_hankelh1(1, k*r)
    h0 = _bim_hankelh1(0, k*r)
    αL2 = im*k/2
    a, da, dda = _ebim_lin1_deriv(αL2*inn, r, invr, h0, h1, k)
    return pb.ws[j]*a, pb.ws[j]*da, pb.ws[j]*dda
end

# Same as above, combined-field cross-component kernel entry plus derivatives
# for a `CombinedFieldIntegralEquationSolver` source component.
@inline function _composite_cross_kernel_entry_with_derivatives(::CombinedFieldIntegralEquationSolver, pb::BoundaryPoints{T}, xi::T, yi::T, k::Union{T,Complex{T}}, j::Int) where {T<:Real}
    xj, yj = pb.xy[j]
    dx = xi-xj
    dy = yi-yj
    r = hypot(dx, dy)
    invr = inv(r)
    tx, ty = pb.tangent[j]
    inn = ty*dx-tx*dy
    sj = hypot(tx, ty)
    ik = im*k
    h0 = _bim_hankelh1(0, k*r)
    h1 = _bim_hankelh1(1, k*r)
    αL2 = im*k/2
    αM2 = Complex{T}(0, one(T)/2)
    Dval, dDval, ddDval = _ebim_lin1_deriv(αL2*inn, r, invr, h0, h1, k)
    Dval *= pb.ws[j]; dDval *= pb.ws[j]; ddDval *= pb.ws[j]
    Sval, dSval, ddSval = _ebim_const0_deriv(αM2*sj, r, k, h0, h1)
    Sval *= pb.ws[j]; dSval *= pb.ws[j]; ddSval *= pb.ws[j]
    val = Dval+ik*Sval
    dval = dDval+im*Sval+ik*dSval
    ddval = ddDval+2*im*dSval+ik*ddSval
    return val, dval, ddval
end

################################################################################
######################### FREDHOLM ASSEMBLY WITH DERIVATIVES ##################
################################################################################

# Full (unfolded) `A(k)=I-D(k)`/`A(k)=I-(D(k)+ikS(k))` Fredholm matrix and its
# first two `k`-derivatives, dispatched through `entry_fn` (either
# `_dlp_kernel_entry_with_derivatives` or `_cfie_kernel_entry_with_derivatives`).
function _ebim_fredholm_full_with_derivatives!(entry_fn, A::AbstractMatrix{Complex{T}}, dA::AbstractMatrix{Complex{T}}, ddA::AbstractMatrix{Complex{T}}, pts::BoundaryPoints{T}, Rmat::AbstractMatrix{T}, G::BoundaryGeomCache{T}, k::Union{T,Complex{T}}; multithreaded::Bool=true) where {T<:Real}
    N = length(pts)
    @use_threads multithreading=multithreaded for j in 1:N
        @inbounds for i in 1:N
            val, dval, ddval = entry_fn(pts, Rmat, G, k, i, j)
            Iij = i==j ? one(Complex{T}) : zero(Complex{T})
            A[i,j] = Iij-val
            dA[i,j] = -dval
            ddA[i,j] = -ddval
        end
    end
    return A, dA, ddA
end

# Symmetry-reduced `A(k)` Fredholm matrix and its first two `k`-derivatives,
# folding the complete discrete full-boundary kernel over each source
# symmetry orbit (mirrors `_dlp_fredholm_reduced!`/`_cfie_fredholm_reduced!`'s
# image-list folding), dispatched through `entry_fn`.
function _ebim_fredholm_reduced_with_derivatives!(entry_fn, A::AbstractMatrix{Complex{T}}, dA::AbstractMatrix{Complex{T}}, ddA::AbstractMatrix{Complex{T}}, pts::BoundaryPoints{T}, Rmat::AbstractMatrix{T}, G::BoundaryGeomCache{T}, orbits::SymmetryOrbitMap{T}, k::Union{T,Complex{T}}; multithreaded::Bool=true) where {T<:Real}
    m = fundamental_size(orbits)
    N = length(orbits)
    fund = orbits.fundamental_indices
    orbit_of = orbits.orbit_of
    phase = orbits.phase
    images = [Int[] for _ in 1:m]
    @inbounds for j in 1:N
        push!(images[orbit_of[j]], j)
    end
    @use_threads multithreading=(multithreaded && m>=32) for b in 1:m
        @inbounds for a in 1:m
            i = fund[a]
            acc = zero(Complex{T})
            dacc = zero(Complex{T})
            ddacc = zero(Complex{T})
            for j in images[b]
                val, dval, ddval = entry_fn(pts, Rmat, G, k, i, j)
                acc += phase[j]*val
                dacc += phase[j]*dval
                ddacc += phase[j]*ddval
            end
            A[a,b] = -acc
            dA[a,b] = -dacc
            ddA[a,b] = -ddacc
        end
        A[b,b] += one(Complex{T})
    end
    return A, dA, ddA
end

# Full (unfolded) composite `A(k)` Fredholm matrix and its first two
# `k`-derivatives (mirrors `_composite_fredholm_full!` in compositebim.jl):
# same-component diagonal blocks reuse the Kress-corrected DLP/CFIE
# derivative kernels, cross-component blocks use the smooth derivative
# kernels above.
function _composite_fredholm_full_with_derivatives!(A::AbstractMatrix{Complex{T}}, dA::AbstractMatrix{Complex{T}}, ddA::AbstractMatrix{Complex{T}}, solver::CompositeBIMSolver, comp_pts::Vector{BoundaryPoints{T}}, Gs::Vector{BoundaryGeomCache{T}}, Rmats::Vector{Matrix{T}}, offs::Vector{Int}, k::Union{T,Complex{T}}; multithreaded::Bool=true) where {T<:Real}
    nc = length(comp_pts)
    @inbounds for a in 1:nc
        cs = solver.component_solvers[a]
        pa = comp_pts[a]
        Ga = Gs[a]
        Ra = Rmats[a]
        Na = length(pa)
        off = offs[a]
        @use_threads multithreading=(multithreaded && Na>=32) for j in 1:Na
            gj = off+j-1
            @inbounds for i in 1:Na
                gi = off+i-1
                val, dval, ddval = _composite_component_kernel_entry_with_derivatives(cs, pa, Ra, Ga, k, i, j)
                Iij = i==j ? one(Complex{T}) : zero(Complex{T})
                A[gi,gj] = Iij-val
                dA[gi,gj] = -dval
                ddA[gi,gj] = -ddval
            end
        end
    end
    for b in 1:nc
        csb = solver.component_solvers[b]
        pb = comp_pts[b]
        offb = offs[b]
        Nb = length(pb)
        for a in 1:nc
            a==b && continue
            pa = comp_pts[a]
            offa = offs[a]
            Na = length(pa)
            @use_threads multithreading=(multithreaded && Na>=16) for i in 1:Na
                gi = offa+i-1
                xi, yi = pa.xy[i]
                @inbounds for j in 1:Nb
                    gj = offb+j-1
                    val, dval, ddval = _composite_cross_kernel_entry_with_derivatives(csb, pb, xi, yi, k, j)
                    A[gi,gj] = -val
                    dA[gi,gj] = -dval
                    ddA[gi,gj] = -ddval
                end
            end
        end
    end
    return A, dA, ddA
end

# Symmetry-reduced composite `A(k)` Fredholm matrix and its first two
# `k`-derivatives (mirrors `_composite_fredholm_reduced!`).
function _composite_fredholm_reduced_with_derivatives!(A::AbstractMatrix{Complex{T}}, dA::AbstractMatrix{Complex{T}}, ddA::AbstractMatrix{Complex{T}}, solver::CompositeBIMSolver, comp_pts::Vector{BoundaryPoints{T}}, Gs::Vector{BoundaryGeomCache{T}}, Rmats::Vector{Matrix{T}}, offs::Vector{Int}, g2c::Vector{Int}, g2l::Vector{Int}, orbits::SymmetryOrbitMap{T}, k::Union{T,Complex{T}}; multithreaded::Bool=true) where {T<:Real}
    m = fundamental_size(orbits)
    N = length(orbits)
    fund = orbits.fundamental_indices
    orbit_of = orbits.orbit_of
    phase = orbits.phase
    images = [Int[] for _ in 1:m]
    @inbounds for j in 1:N
        push!(images[orbit_of[j]], j)
    end
    @use_threads multithreading=(multithreaded && m>=32) for b in 1:m
        @inbounds for a in 1:m
            gi = fund[a]
            ca = g2c[gi]
            ia = g2l[gi]
            acc = zero(Complex{T})
            dacc = zero(Complex{T})
            ddacc = zero(Complex{T})
            for gj in images[b]
                cb = g2c[gj]
                jb = g2l[gj]
                ph = phase[gj]
                if ca==cb
                    cs = solver.component_solvers[ca]
                    val, dval, ddval = _composite_component_kernel_entry_with_derivatives(cs, comp_pts[ca], Rmats[ca], Gs[ca], k, ia, jb)
                else
                    csb = solver.component_solvers[cb]
                    xi, yi = comp_pts[ca].xy[ia]
                    val, dval, ddval = _composite_cross_kernel_entry_with_derivatives(csb, comp_pts[cb], xi, yi, k, jb)
                end
                acc += ph*val
                dacc += ph*dval
                ddacc += ph*ddval
            end
            A[a,b] = -acc
            dA[a,b] = -dacc
            ddA[a,b] = -ddacc
        end
        A[b,b] += one(Complex{T})
    end
    return A, dA, ddA
end

################################################################################
################### PER-KERNEL construct_matrices DISPATCH ###################
################################################################################

function _ebim_construct_matrices(cs::DoubleLayerPotentialSolver, pts::BoundaryPoints{T}, k; multithreaded::Bool=true) where {T<:Real}
    kT = _bim_widen_k(T, k)
    N = length(pts)
    graded = _is_nontrivial_dlp_grading(pts)
    G = boundary_geom_cache(pts, graded)
    Rmat = zeros(T, N, N)
    kress_R!(Rmat)
    if cs.symmetry===nothing
        A = Matrix{Complex{T}}(undef, N, N)
        dA = similar(A)
        ddA = similar(A)
        _ebim_fredholm_full_with_derivatives!(_dlp_kernel_entry_with_derivatives, A, dA, ddA, pts, Rmat, G, kT; multithreaded)
        return A, dA, ddA
    end
    orbits = symmetry_index_orbits(T, pts.xy, cs.symmetry)
    m = fundamental_size(orbits)
    A = Matrix{Complex{T}}(undef, m, m)
    dA = similar(A)
    ddA = similar(A)
    _ebim_fredholm_reduced_with_derivatives!(_dlp_kernel_entry_with_derivatives, A, dA, ddA, pts, Rmat, G, orbits, kT; multithreaded)
    return A, dA, ddA
end

function _ebim_construct_matrices(cs::CombinedFieldIntegralEquationSolver, pts::BoundaryPoints{T}, k; multithreaded::Bool=true) where {T<:Real}
    kT = _bim_widen_k(T, k)
    N = length(pts)
    graded = _is_nontrivial_dlp_grading(pts)
    G = boundary_geom_cache(pts, graded)
    Rmat = zeros(T, N, N)
    kress_R!(Rmat)
    if cs.symmetry===nothing
        A = Matrix{Complex{T}}(undef, N, N)
        dA = similar(A)
        ddA = similar(A)
        _ebim_fredholm_full_with_derivatives!(_cfie_kernel_entry_with_derivatives, A, dA, ddA, pts, Rmat, G, kT; multithreaded)
        return A, dA, ddA
    end
    orbits = symmetry_index_orbits(T, pts.xy, cs.symmetry)
    m = fundamental_size(orbits)
    A = Matrix{Complex{T}}(undef, m, m)
    dA = similar(A)
    ddA = similar(A)
    _ebim_fredholm_reduced_with_derivatives!(_cfie_kernel_entry_with_derivatives, A, dA, ddA, pts, Rmat, G, orbits, kT; multithreaded)
    return A, dA, ddA
end

function _ebim_construct_matrices(cs::CompositeBIMSolver{T}, pts::BoundaryPoints{T}, k; multithreaded::Bool=true) where {T<:Real}
    kT = _bim_widen_k(T, k)
    nc = length(cs.component_solvers)
    N = length(pts)
    offs = _composite_offsets(pts, nc)
    comp_pts = [_composite_component_slice(pts, offs[a]:offs[a+1]-1, a) for a in 1:nc]
    Gs = Vector{BoundaryGeomCache{T}}(undef, nc)
    Rmats = Vector{Matrix{T}}(undef, nc)
    @inbounds for a in 1:nc
        graded = _is_nontrivial_dlp_grading(comp_pts[a])
        Gs[a] = boundary_geom_cache(comp_pts[a], graded)
        Na = length(comp_pts[a])
        Ra = zeros(T, Na, Na)
        kress_R!(Ra)
        Rmats[a] = Ra
    end
    if cs.symmetry===nothing
        A = Matrix{Complex{T}}(undef, N, N)
        dA = similar(A)
        ddA = similar(A)
        _composite_fredholm_full_with_derivatives!(A, dA, ddA, cs, comp_pts, Gs, Rmats, offs, kT; multithreaded)
        return A, dA, ddA
    end
    orbits = symmetry_index_orbits(T, pts.xy, cs.symmetry)
    m = fundamental_size(orbits)
    g2c, g2l = _composite_global_to_local(offs)
    A = Matrix{Complex{T}}(undef, m, m)
    dA = similar(A)
    ddA = similar(A)
    _composite_fredholm_reduced_with_derivatives!(A, dA, ddA, cs, comp_pts, Gs, Rmats, offs, g2c, g2l, orbits, kT; multithreaded)
    return A, dA, ddA
end

################################################################################
############################## PUBLIC API ######################################
################################################################################

"""
    construct_matrices(solver::ExpandedBIMSolver, pts::BoundaryPoints, k; multithreaded::Bool = true) → (A::Matrix, dA::Matrix, ddA::Matrix)

Assembles the Fredholm matrix `A(k)` and its first two `k`-derivatives
`A'(k)`, `A''(k)` by direct (non-Chebyshev) evaluation of `solver.kernel`'s
Kress-corrected Fredholm kernel and its Hankel/Bessel-function `k`-derivatives
(see [`_dlp_kernel_entry_with_derivatives`](@ref),
[`_cfie_kernel_entry_with_derivatives`](@ref)).
"""
function construct_matrices(solver::ExpandedBIMSolver, pts::BoundaryPoints, k; multithreaded::Bool=true)
    solver.use_chebyshev && error("Chebyshev-accelerated EBIM evaluation not yet implemented, see the QuantumBilliardsTests migration plan step 10. Construct the ExpandedBIMSolver with use_chebyshev=false.")
    return _ebim_construct_matrices(solver.kernel, pts, k; multithreaded)
end

"""
    solve(solver::ExpandedBIMSolver, pts::BoundaryPoints, k; multithreaded::Bool = true) → (k_corr::Real, t0::Real)

Computes the second-order locally-corrected root `k_corr` near `k` and its
tension `t0`.

## Description
Assembles `(A,dA,ddA)=A(k),A'(k),A''(k)` via [`construct_matrices`](@ref),
solves the dense generalized eigenproblem `A v = λ A' v` (`LinearAlgebra.eigen`
on the matrix pencil), and keeps the eigenpair `(λ,v)` of smallest `|λ|` (the
locally dominant root). The corresponding left eigenvector `u` is obtained
from the adjoint pencil `A' u = μ (A')' u` (`eigen(A',dA')`, eigenvalues
`μ≈conj(λ)`). The correction is `ε₁=-λ`,
`ε₂=-(1/2)ε₁²[u'A''(k)v]/[u'A'(k)v]`, giving `k_corr=k+Re(ε₁+ε₂)` and
`t0=|ε₁+ε₂|`.
"""
function solve(solver::ExpandedBIMSolver, pts::BoundaryPoints, k; multithreaded::Bool=true)
    T = _bim_numeric_type(solver)
    A, dA, ddA = construct_matrices(solver, pts, k; multithreaded)
    @blas_multi_then_1 MAX_BLAS_THREADS Fr = eigen(A, dA)
    λ = Fr.values
    V = Fr.vectors
    jr = argmin(abs.(λ))
    λj = λ[jr]
    v = @view V[:,jr]
    @blas_multi_then_1 MAX_BLAS_THREADS Fl = eigen(A', dA')
    μ = Fl.values
    U = Fl.vectors
    jl = argmin(abs.(μ.-conj(λj)))
    u = @view U[:,jl]
    ε1 = -λj
    num = dot(u, ddA*v)
    den = dot(u, dA*v)
    ε2 = abs(den)>eps(T) ? -T(0.5)*ε1^2*(num/den) : zero(ε1)
    corr = ε1+ε2
    return real(k+corr), abs(corr)
end

"""
    solve_wavenumber(solver::ExpandedBIMSolver, billiard::Bi, k, dk; multithreaded::Bool = true) where {Bi<:AbsBilliard} → (k0::Real, t0::Real)

Computes the second-order locally-corrected root nearest `k` (`dk` is retained
for API parity with [`solve_wavenumber(::BeynSolver, ...)`](@ref) but is
unused by the local expansion, which needs no window).
"""
function solve_wavenumber(solver::ExpandedBIMSolver, billiard::Bi, k, dk; multithreaded::Bool=true) where {Bi<:AbsBilliard}
    pts = evaluate_points(solver, billiard, k)
    return solve(solver, pts, k; multithreaded)
end

"""
    solve_spectrum(solver::ExpandedBIMSolver, billiard::Bi, k, dk; multithreaded::Bool = true) where {Bi<:AbsBilliard} → (ks::Vector, ts::Vector)

Computes the second-order locally-corrected roots for every wavenumber in a
sweep, one call per target `k` (`k` is a vector of trial wavenumbers, e.g. a
prior sweep's approximate roots; `dk` is unused, see [`solve_wavenumber`](@ref)).
"""
function solve_spectrum(solver::ExpandedBIMSolver, billiard::Bi, k, dk; multithreaded::Bool=true) where {Bi<:AbsBilliard}
    T = _bim_numeric_type(solver)
    n = length(k)
    ks = Vector{T}(undef, n)
    ts = Vector{T}(undef, n)
    @inbounds for i in 1:n
        ks[i], ts[i] = solve_wavenumber(solver, billiard, k[i], dk; multithreaded)
    end
    return ks, ts
end