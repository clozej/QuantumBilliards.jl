################################################################################
# EXPANDED BOUNDARY INTEGRAL METHOD
#
# For a nonlinear Fredholm matrix A(k), EBIM uses the local expansion
#
#     A(k+ε)=A(k)+εA'(k)+1/2 ε²A''(k)+O(ε³).
#
# The generalized eigenproblem
#
#     A(k)v=λA'(k)v
#
# gives the first-order root correction ε₁=-λ. With the corresponding left
# generalized eigenvector u, the second-order correction is
#
#     ε₂=-1/2 ε₁² [u†A''(k)v]/[u†A'(k)v].
#
# Hence
#
#     k_corr=k+ε₁+ε₂.
#
# Geometry, symmetry reduction, quadrature, Chebyshev plans and solver-specific
# workspaces are delegated completely to the boundary-matrix backends.
################################################################################

"""
    EBIMSolver

Boundary-integral solver backends supporting the production EBIM workflow.

Every solver must provide direct derivative matrix construction and the common
derivative-aware Chebyshev workspace interface used by EBIM.
"""
const EBIMSolver=Union{BoundaryIntegralMethod,DLP_kress,DLP_kress_global_corners,CFIE_kress,CFIE_kress_corners,CFIE_kress_global_corners,CFIE_kress_composite_solver}

################################################################################
######################## CHEBYSHEV PATHWAY #####################################
################################################################################

"""
    EBIMChebBatchCache{W}

Reusable derivative-aware Chebyshev cache for an EBIM segment.

## Attributes
* `ws::W`: Solver-specific derivative Chebyshev workspace.
* `ks::Vector{ComplexF64}`: Complex wavenumbers represented by the workspace.
"""
struct EBIMChebBatchCache{W}
    ws::W
    ks::Vector{ComplexF64}
end

@inline _ebim_complex_ks(ks)=ComplexF64.(ks)

"""
    build_ebim_cheb_cache(solver::EBIMSolver,pts,ks;n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false) → EBIMChebBatchCache

Build the reusable derivative Chebyshev workspace for one EBIM wavenumber
batch. Solver-specific cache construction is delegated to
`build_derivative_chebyshev_workspace`.

## Arguments
* `solver::EBIMSolver`: Active boundary-integral solver.
* `pts`: Boundary discretization used by the solver.
* `ks`: Wavenumbers represented by the cache.

## Keyword Arguments
* `n_panels_h::Int=15000`: Hankel radial Chebyshev panel count.
* `M_h::Int=5`: Hankel Chebyshev polynomial degree.
* `n_panels_j::Int=10000`: Bessel-J radial Chebyshev panel count.
* `M_j::Int=5`: Bessel-J Chebyshev polynomial degree.
* `timeit::Bool=false`: Enable workspace-construction timing.

## Returns
* `cache::EBIMChebBatchCache`: Reusable solver-specific derivative workspace and its wavenumber grid.
"""
function build_ebim_cheb_cache(solver::EBIMSolver,pts,ks;n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,timeit::Bool=false)
    zks=_ebim_complex_ks(ks)
    ws=build_derivative_chebyshev_workspace(solver,pts,zks;n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,timeit=timeit)
    return EBIMChebBatchCache(ws,zks)
end

@inline function allocate_ebim_matrices(solver::EBIMSolver,pts)
    N=boundary_matrix_size(solver,pts)
    return Matrix{ComplexF64}(undef,N,N),Matrix{ComplexF64}(undef,N,N),Matrix{ComplexF64}(undef,N,N)
end

@inline function allocate_ebim_cheb_matrices(solver::EBIMSolver,pts,cache::EBIMChebBatchCache)
    N=boundary_matrix_size(solver,pts)
    Mk=length(cache.ks)
    As=[Matrix{ComplexF64}(undef,N,N) for _ in 1:Mk]
    dAs=[Matrix{ComplexF64}(undef,N,N) for _ in 1:Mk]
    ddAs=[Matrix{ComplexF64}(undef,N,N) for _ in 1:Mk]
    return As,dAs,ddAs
end

function construct_ebim_cheb_matrices!(As::Vector{Matrix{ComplexF64}},dAs::Vector{Matrix{ComplexF64}},ddAs::Vector{Matrix{ComplexF64}},solver::EBIMSolver,pts,cache::EBIMChebBatchCache;multithreaded::Bool=true)
    construct_matrices_chebyshev_with_derivatives!(As,dAs,ddAs,solver,pts,cache.ws;multithreaded=multithreaded)
    return nothing
end

function construct_ebim_cheb_matrices(solver::EBIMSolver,pts,cache::EBIMChebBatchCache;multithreaded::Bool=true)
    As,dAs,ddAs=allocate_ebim_cheb_matrices(solver,pts,cache)
    construct_ebim_cheb_matrices!(As,dAs,ddAs,solver,pts,cache;multithreaded=multithreaded)
    return As,dAs,ddAs
end

function construct_ebim_cheb_matrix_at!(A::Matrix{ComplexF64},dA::Matrix{ComplexF64},ddA::Matrix{ComplexF64},solver::EBIMSolver,pts,cache::EBIMChebBatchCache,idx::Int;multithreaded::Bool=true)
    1<=idx<=length(cache.ks)||throw(BoundsError(cache.ks,idx))
    construct_matrix_chebyshev_with_derivatives_at!(A,dA,ddA,solver,pts,cache.ws,idx;multithreaded=multithreaded)
    return nothing
end

################################################################################
######################## EBIM CORRECTION #######################################
################################################################################

@inline function _ebim_second_order_correction!(buf,dA,ddA,u,v,λ)
    mul!(buf,ddA,v)
    num=dot(u,buf)
    mul!(buf,dA,v)
    den=dot(u,buf)
    c1=-λ
    c2=abs(den)>1e-15 ? -0.5*c1^2*(num/den) : zero(c1)
    return c1,c2
end

################################################################################
######################## DENSE GENERALIZED EVP #################################
################################################################################

"""
    solve_full!(solver::EBIMSolver,A::AbstractMatrix{Complex{T}},dA::AbstractMatrix{Complex{T}},ddA::AbstractMatrix{Complex{T}},pts,k,dk;use_lapack_raw::Bool=false,multithreaded::Bool=true,return_imag_part::Bool=false) where {T<:Real} → Tuple

Solve the dense EBIM generalized eigenproblem

    A(k)v=λA'(k)v

and apply the second-order correction

    ε₁=-λ,
    ε₂=-1/2 ε₁² [u†A''(k)v]/[u†A'(k)v].

Only generalized eigenvalues satisfying `|Re λ|<dk` and `|Im λ|<dk` are
retained.

## Arguments
* `solver::EBIMSolver`: Active solver, retained for the common EBIM interface.
* `A::AbstractMatrix{Complex{T}}`: Fredholm matrix `A(k)`.
* `dA::AbstractMatrix{Complex{T}}`: First derivative `A'(k)`.
* `ddA::AbstractMatrix{Complex{T}}`: Second derivative `A''(k)`.
* `pts`: Boundary discretization, retained for the common EBIM interface.
* `k`: Trial wavenumber.
* `dk`: Local correction half-width.

## Keyword Arguments
* `use_lapack_raw::Bool=false`: Use the low-level generalized eigensolver backend.
* `multithreaded::Bool=true`: Retained for interface compatibility.
* `return_imag_part::Bool=false`: Return complex corrected roots instead of only their real parts.

## Returns
* `λout::Vector`: Corrected local wavenumber estimates.
* `tens::Vector`: Absolute total corrections `|ε₁+ε₂|`.
"""
function solve_full!(solver::EBIMSolver,A::AbstractMatrix{Complex{T}},dA::AbstractMatrix{Complex{T}},ddA::AbstractMatrix{Complex{T}},pts,k,dk;use_lapack_raw::Bool=false,multithreaded::Bool=true,return_imag_part::Bool=false) where {T<:Real}
    if use_lapack_raw
        @blas_multi MAX_BLAS_THREADS λ,VR,VL=generalized_eigen_all_LAPACK_LEGACY(A,dA)
    else
        @blas_multi MAX_BLAS_THREADS λ,VR,VL=generalized_eigen_all(A,dA)
    end
    CT=eltype(λ)
    RT=real(CT)
    KT=return_imag_part ? CT : RT
    valid=isfinite.(λ).&(abs.(real.(λ)).<dk).&(abs.(imag.(λ)).<dk)
    any(valid)||return KT[],RT[]
    λ=λ[valid]
    VR=VR[:,valid]
    VL=VL[:,valid]
    nk=length(λ)
    λout=Vector{KT}(undef,nk)
    tens=Vector{RT}(undef,nk)
    buf=Vector{CT}(undef,size(A,1))
    @inbounds for j in 1:nk
        c1,c2=_ebim_second_order_correction!(buf,dA,ddA,@view(VL[:,j]),@view(VR[:,j]),λ[j])
        t=c1+c2
        kc=complex(k)+t
        λout[j]=return_imag_part ? kc : real(kc)
        tens[j]=abs(t)
    end
    return λout,tens
end

################################################################################
######################## KRYLOV GENERALIZED EVP ################################
################################################################################

"""
    solve_krylov!(solver::EBIMSolver,A::AbstractMatrix{Complex{T}},dA::AbstractMatrix{Complex{T}},ddA::AbstractMatrix{Complex{T}},pts,k,dk;multithreaded::Bool=true,nev::Int=5,tol=1e-14,maxiter::Int=5000,krylovdim::Int=max(40,2*nev+1),return_imag_part::Bool=false) where {T<:Real} → Tuple

Solve the EBIM generalized eigenproblem with shift-invert Krylov iteration.

The right operator is

    C=A⁻¹A',

with eigenvalues `μ=1/λ`, so the smallest generalized corrections `|λ|` become
the largest-magnitude eigenvalues of `C`.

The corresponding left generalized eigenvectors are obtained from

    C_L=(A†)⁻¹A'†,

whose eigenvalues are `conj(μ)`. The left and right Krylov spectra are therefore
paired through `μ_l≈conj(μ_r)`.

`A` is overwritten by its LU factorization.

## Arguments
* `solver::EBIMSolver`: Active solver, retained for the common EBIM interface.
* `A::AbstractMatrix{Complex{T}}`: Fredholm matrix `A(k)`, overwritten by `lu!`.
* `dA::AbstractMatrix{Complex{T}}`: First derivative `A'(k)`.
* `ddA::AbstractMatrix{Complex{T}}`: Second derivative `A''(k)`.
* `pts`: Boundary discretization, retained for the common EBIM interface.
* `k`: Trial wavenumber.
* `dk`: Local correction half-width.

## Keyword Arguments
* `multithreaded::Bool=true`: Retained for interface compatibility.
* `nev::Int=5`: Number of right and left Krylov eigenpairs requested.
* `tol=1e-14`: Krylov convergence tolerance.
* `maxiter::Int=5000`: Maximum Krylov iterations.
* `krylovdim::Int=max(40,2*nev+1)`: Krylov subspace dimension.
* `return_imag_part::Bool=false`: Return complex corrected roots instead of only their real parts.

## Returns
* `λout::Vector`: Accepted corrected local wavenumbers.
* `tens::Vector`: Absolute total corrections `|ε₁+ε₂|`.
"""
function solve_krylov!(solver::EBIMSolver,A::AbstractMatrix{Complex{T}},dA::AbstractMatrix{Complex{T}},ddA::AbstractMatrix{Complex{T}},pts,k,dk;multithreaded::Bool=true,nev::Int=5,tol=1e-14,maxiter::Int=5000,krylovdim::Int=max(40,2*nev+1),return_imag_part::Bool=false) where {T<:Real}
    CT=eltype(A)
    RT=real(CT)
    KT=return_imag_part ? CT : RT
    n=size(A,1)
    nev=min(nev,n)
    @blas_multi MAX_BLAS_THREADS F=lu!(A)
    Ft=adjoint(F)
    dAt=adjoint(dA)
    function op_r!(y,x)
        mul!(y,dA,x)
        ldiv!(F,y)
        return y
    end
    function op_l!(y,x)
        mul!(y,dAt,x)
        ldiv!(Ft,y)
        return y
    end
    C=LinearMaps.LinearMap{CT}(op_r!,n,n;ismutating=true)
    Cl=LinearMaps.LinearMap{CT}(op_l!,n,n;ismutating=true)
    μr,VR,_=eigsolve(C,n,nev,:LM;tol=tol,maxiter=maxiter,krylovdim=krylovdim)
    μl,UL,_=eigsolve(Cl,n,nev,:LM;tol=tol,maxiter=maxiter,krylovdim=krylovdim)
    λ=inv.(μr)
    ord=sortperm(abs.(λ))
    λ=λ[ord]
    μr=μr[ord]
    VR=VR[ord]
    perm=[argmin(abs.(μl.-conj(μ))) for μ in μr]
    UL=UL[perm]
    λout=Vector{KT}(undef,nev)
    tens=Vector{RT}(undef,nev)
    buf=Vector{CT}(undef,n)
    m=0
    @inbounds for j in 1:nev
        λj=λ[j]
        abs(real(λj))<dk&&abs(imag(λj))<dk||continue
        c1,c2=_ebim_second_order_correction!(buf,dA,ddA,UL[j],VR[j],λj)
        t=c1+c2
        kc=complex(k)+t
        abs(real(kc)-k)<dk||continue
        m+=1
        λout[m]=return_imag_part ? kc : real(kc)
        tens[m]=abs(t)
    end
    resize!(λout,m)
    resize!(tens,m)
    return λout,tens
end

################################################################################
######################## SOLVE INTERFACE #######################################
################################################################################

"""
    solve!(solver::EBIMSolver,A::AbstractMatrix{Complex{T}},dA::AbstractMatrix{Complex{T}},ddA::AbstractMatrix{Complex{T}},pts,k,dk;use_lapack_raw::Bool=false,multithreaded::Bool=true,use_krylov::Bool=true,nev::Int=5,return_imag_part::Bool=false) where {T<:Real} → Tuple

Construct `A(k)`, `A'(k)` and `A''(k)` directly and perform one local EBIM
solve.

## Arguments
* `solver::EBIMSolver`: Active boundary-integral solver.
* `A::AbstractMatrix{Complex{T}}`: Reusable Fredholm matrix buffer.
* `dA::AbstractMatrix{Complex{T}}`: Reusable first-derivative buffer.
* `ddA::AbstractMatrix{Complex{T}}`: Reusable second-derivative buffer.
* `pts`: Boundary discretization.
* `k`: Trial wavenumber.
* `dk`: Local correction half-width.

## Keyword Arguments
* `use_lapack_raw::Bool=false`: Use the low-level dense generalized eigensolver.
* `multithreaded::Bool=true`: Enable threaded matrix assembly.
* `use_krylov::Bool=true`: Select Krylov or dense generalized eigensolution.
* `nev::Int=5`: Number of Krylov eigenpairs requested.
* `return_imag_part::Bool=false`: Return complex corrected roots when enabled.

## Returns
* `λs::Vector`: Corrected local roots.
* `tensions::Vector`: Corresponding absolute EBIM corrections.
"""
function solve!(solver::EBIMSolver,A::AbstractMatrix{Complex{T}},dA::AbstractMatrix{Complex{T}},ddA::AbstractMatrix{Complex{T}},pts,k,dk;use_lapack_raw::Bool=false,multithreaded::Bool=true,use_krylov::Bool=true,nev::Int=5,return_imag_part::Bool=false) where {T<:Real}
    basis=AbstractHankelBasis()
    @blas_1 construct_matrices!(solver,basis,A,dA,ddA,pts,k;multithreaded=multithreaded)
    if use_krylov
        return solve_krylov!(solver,A,dA,ddA,pts,k,dk;multithreaded=multithreaded,nev=nev,return_imag_part=return_imag_part)
    end
    return solve_full!(solver,A,dA,ddA,pts,k,dk;use_lapack_raw=use_lapack_raw,multithreaded=multithreaded,return_imag_part=return_imag_part)
end

"""
    solve!(solver::EBIMSolver,A::Matrix{ComplexF64},dA::Matrix{ComplexF64},ddA::Matrix{ComplexF64},pts,k,dk,cache::EBIMChebBatchCache,idx::Int;use_lapack_raw::Bool=false,multithreaded::Bool=true,use_krylov::Bool=true,nev::Int=5,return_imag_part::Bool=false) → Tuple

Construct one matrix triple from a reusable derivative Chebyshev cache and
perform the local EBIM solve without rebuilding geometry or interpolation
plans.

## Arguments
* `solver::EBIMSolver`: Active boundary-integral solver.
* `A::Matrix{ComplexF64}`: Reusable Fredholm matrix buffer.
* `dA::Matrix{ComplexF64}`: Reusable first-derivative buffer.
* `ddA::Matrix{ComplexF64}`: Reusable second-derivative buffer.
* `pts`: Boundary discretization.
* `k`: Trial wavenumber.
* `dk`: Local correction half-width.
* `cache::EBIMChebBatchCache`: Reusable derivative Chebyshev cache.
* `idx::Int`: Wavenumber index inside `cache`.

## Keyword Arguments
* `use_lapack_raw::Bool=false`: Use the low-level dense generalized eigensolver.
* `multithreaded::Bool=true`: Enable threaded cached matrix assembly.
* `use_krylov::Bool=true`: Select Krylov or dense generalized eigensolution.
* `nev::Int=5`: Number of Krylov eigenpairs requested.
* `return_imag_part::Bool=false`: Return complex corrected roots when enabled.

## Returns
* `λs::Vector`: Corrected local roots.
* `tensions::Vector`: Corresponding absolute EBIM corrections.
"""
function solve!(solver::EBIMSolver,A::Matrix{ComplexF64},dA::Matrix{ComplexF64},ddA::Matrix{ComplexF64},pts,k,dk,cache::EBIMChebBatchCache,idx::Int;use_lapack_raw::Bool=false,multithreaded::Bool=true,use_krylov::Bool=true,nev::Int=5,return_imag_part::Bool=false)
    construct_ebim_cheb_matrix_at!(A,dA,ddA,solver,pts,cache,idx;multithreaded=multithreaded)
    if use_krylov
        return solve_krylov!(solver,A,dA,ddA,pts,k,dk;multithreaded=multithreaded,nev=nev,return_imag_part=return_imag_part)
    end
    return solve_full!(solver,A,dA,ddA,pts,k,dk;use_lapack_raw=use_lapack_raw,multithreaded=multithreaded,return_imag_part=return_imag_part)
end

################################################################################
######################## OVERLAP MERGING #######################################
################################################################################

@inline function _ebim_scoring_logic(k,t)
    return log10(abs(imag(ComplexF64(k)))+eps(Float64))+log10(abs(Float64(t))+eps(Float64))
end

function _local_gap(xs,i;w::Int=4)
    i1=max(1,i-w)
    i2=min(length(xs),i+w)
    gaps=diff(@view xs[i1:i2])
    gaps=gaps[gaps.>0]
    return isempty(gaps) ? Inf : median(gaps)
end

"""
    overlap_and_merge_ebim!(k_left::Vector{K},ten_left::Vector{T},k_right::Vector{K},ten_right::Vector{T},control_left::Vector{Bool};tol::T=T(1e-5),spacing_frac::T=T(0.02),tolmax::T=T(5e-3),local_window::Int=4) where {K<:Number,T<:Real} → Nothing

Merge duplicate root candidates generated by overlapping EBIM windows.

Neighboring roots are clustered using an adaptive tolerance based on the local
median spectral spacing. Inside each cluster the candidate minimizing

    log10(|Im k|+eps)+log10(|tension|+eps)

is retained.

## Arguments
* `k_left::Vector{K}`: Accumulated roots, modified in place.
* `ten_left::Vector{T}`: Accumulated tensions, modified in place.
* `k_right::Vector{K}`: New roots to merge.
* `ten_right::Vector{T}`: New root tensions.
* `control_left::Vector{Bool}`: Merge-control flags, modified in place.

## Keyword Arguments
* `tol::T=T(1e-5)`: Minimum clustering tolerance.
* `spacing_frac::T=T(0.02)`: Fraction of the local spectral spacing used for clustering.
* `tolmax::T=T(5e-3)`: Maximum clustering tolerance.
* `local_window::Int=4`: Local spacing window half-width.

## Returns
* `nothing`.
"""
function overlap_and_merge_ebim!(k_left::Vector{K},ten_left::Vector{T},k_right::Vector{K},ten_right::Vector{T},control_left::Vector{Bool};tol::T=T(1e-5),spacing_frac::T=T(0.02),tolmax::T=T(5e-3),local_window::Int=4) where {K<:Number,T<:Real}
    isempty(k_right)&&return nothing
    append!(k_left,k_right)
    append!(ten_left,ten_right)
    append!(control_left,fill(false,length(k_right)))
    p=sortperm(real.(k_left))
    k_all=k_left[p]
    ten_all=ten_left[p]
    ctrl_all=control_left[p]
    xs=real.(k_all)
    new_k=K[]
    new_t=T[]
    new_c=Bool[]
    i=1
    while i<=length(k_all)
        j=i
        while j<length(k_all)
            gap=xs[j+1]-xs[j]
            lgap=min(_local_gap(xs,j;w=local_window),_local_gap(xs,j+1;w=local_window))
            local_tol=min(tolmax,max(tol,spacing_frac*T(lgap)))
            gap<=local_tol||break
            j+=1
        end
        block=i:j
        best=block[argmin(_ebim_scoring_logic.(k_all[block],ten_all[block]))]
        push!(new_k,k_all[best])
        push!(new_t,ten_all[best])
        push!(new_c,any(ctrl_all[block])||length(block)>1)
        i=j+1
    end
    empty!(k_left)
    empty!(ten_left)
    empty!(control_left)
    append!(k_left,new_k)
    append!(ten_left,new_t)
    append!(control_left,new_c)
    return nothing
end

################################################################################
######################## SPECTRUM SWEEP ########################################
################################################################################

"""
    solve_spectrum_ebim(solver::EBIMSolver,billiard::Bi,k1::T,k2::T;dk::Function=(k->0.05*k^(-1/3)),tol::T=T(1e-5),use_lapack_raw::Bool=false,multithreaded_matrices::Bool=false,use_krylov::Bool=true,seg_reuse_frac::T=T(0.95),use_chebyshev::Bool=false,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,cheb_param_strategy::Symbol=:global,cheb_tol::Real=1e-13,max_iter::Int=20,sampling_points::Int=50000,grow_panels::Real=1.5,grow_M::Int=2,verbose_cheb_panelization::Bool=false,return_imag_part::Bool=false) where {T<:Real,Bi<:AbsBilliard} → Tuple

Compute the EBIM spectrum on `[k1,k2]`.

Trial wavenumbers are grouped into segments. Each segment uses the boundary
geometry evaluated at its highest wavenumber and reuses that discretization for
all lower trial points in the segment.

With `use_chebyshev=true`, one derivative Chebyshev workspace is built per
segment and reused for all local solves. Chebyshev parameters may be chosen
globally, independently for each segment, or supplied manually.

The number of Krylov eigenpairs requested in each local window is estimated
from the leading Weyl density with additional padding. Overlapping local root
sets are merged using the adaptive EBIM clustering criterion.

## Arguments
* `solver::EBIMSolver`: Boundary-integral solver.
* `billiard::Bi`: Billiard geometry.
* `k1::T`: Lower wavenumber bound.
* `k2::T`: Upper wavenumber bound.

## Keyword Arguments
* `dk::Function=(k->0.05*k^(-1/3))`: Positive local grid-spacing function.
* `tol::T=T(1e-5)`: Minimum overlap-merging tolerance.
* `use_lapack_raw::Bool=false`: Use the low-level dense generalized eigensolver.
* `multithreaded_matrices::Bool=false`: Enable threaded boundary-matrix assembly.
* `use_krylov::Bool=true`: Use shift-invert Krylov instead of the dense generalized EVP.
* `seg_reuse_frac::T=T(0.95)`: Geometry-reuse fraction controlling segment size.
* `use_chebyshev::Bool=false`: Use reusable derivative Chebyshev matrix construction.
* `n_panels_h::Int=15000`: Initial/manual Hankel Chebyshev panel count.
* `M_h::Int=5`: Initial/manual Hankel Chebyshev degree.
* `n_panels_j::Int=10000`: Initial/manual Bessel-J Chebyshev panel count.
* `M_j::Int=5`: Initial/manual Bessel-J Chebyshev degree.
* `cheb_param_strategy::Symbol=:global`: Chebyshev parameter strategy `:global`, `:segment`, or `:manual`.
* `cheb_tol::Real=1e-13`: Target Chebyshev interpolation error.
* `max_iter::Int=20`: Maximum Chebyshev tuning iterations.
* `sampling_points::Int=50000`: Validation radii used during Chebyshev tuning.
* `grow_panels::Real=1.5`: Panel-count growth factor during tuning.
* `grow_M::Int=2`: Polynomial-degree increment during tuning.
* `verbose_cheb_panelization::Bool=false`: Print Chebyshev tuning diagnostics.
* `return_imag_part::Bool=false`: Return complex corrected roots rather than only their real parts.

## Returns
* `λs::Vector{T}` or `Vector{Complex{T}}`: Merged corrected roots inside `[k1,k2]`.
* `tensions::Vector{T}`: Corresponding EBIM tensions.
"""
function solve_spectrum_ebim(solver::EBIMSolver,billiard::Bi,k1::T,k2::T;dk::Function=(k->0.05*k^(-1/3)),tol::T=T(1e-5),use_lapack_raw::Bool=false,multithreaded_matrices::Bool=false,use_krylov::Bool=true,seg_reuse_frac::T=T(0.95),use_chebyshev::Bool=false,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,cheb_param_strategy::Symbol=:global,cheb_tol::Real=1e-13,max_iter::Int=20,sampling_points::Int=50_000,grow_panels::Real=1.5,grow_M::Int=2,verbose_cheb_panelization::Bool=false,return_imag_part::Bool=false) where {T<:Real,Bi<:AbsBilliard}
    k1<k2||throw(ArgumentError("require k1<k2"))
    0<seg_reuse_frac<=1||throw(ArgumentError("seg_reuse_frac must satisfy 0<seg_reuse_frac<=1"))
    cheb_param_strategy in (:global,:segment,:manual)||throw(ArgumentError("cheb_param_strategy must be :global, :segment or :manual"))
    ks=T[]
    dks=T[]
    k=k1
    while k<k2
        Δk=T(dk(k))
        Δk>0||throw(ArgumentError("dk(k) must be positive; received dk($k)=$Δk"))
        push!(ks,k)
        push!(dks,Δk)
        k+=Δk
    end
    K=return_imag_part ? Complex{T} : T
    isempty(ks)&&return K[],T[]
    nevs=Vector{Int}(undef,length(ks))
    @inbounds for i in eachindex(ks)
        nevs[i]=max(1,ceil(Int,(fundamental_area(billiard)*ks[i]/(2*pi))*dks[i])+10)
    end
    if use_chebyshev&&cheb_param_strategy==:global
        kref=ks[end]
        ptsref=evaluate_points(solver,billiard,kref)
        out=chebyshev_params(solver,ptsref,ComplexF64[kref];tol=cheb_tol,npanels_h_init=n_panels_h,M_h_init=M_h,npanels_j_init=n_panels_j,M_j_init=M_j,sampling_points=sampling_points,max_iter=max_iter,grow_panels=grow_panels,grow_M=grow_M,verbose=verbose_cheb_panelization)
        n_panels_h,M_h,n_panels_j,M_j=out[1],out[2],out[3],out[4]
    end
    results=Vector{Tuple{Vector{K},Vector{T}}}(undef,length(ks))
    p=Progress(length(ks),1)
    pts0=evaluate_points(solver,billiard,ks[1])
    seg_first=1
    while seg_first<=length(ks)
        seg_last=seg_first
        while seg_last<length(ks)&&ks[seg_last+1]<=ks[seg_first]/seg_reuse_frac
            seg_last+=1
        end
        pts=seg_last==1 ? pts0 : evaluate_points(solver,billiard,ks[seg_last])
        A,dA,ddA=allocate_ebim_matrices(solver,pts)
        if use_chebyshev
            segks=@view ks[seg_first:seg_last]
            if cheb_param_strategy==:segment
                out=chebyshev_params(solver,pts,ComplexF64[segks[end]];tol=cheb_tol,npanels_h_init=n_panels_h,M_h_init=M_h,npanels_j_init=n_panels_j,M_j_init=M_j,sampling_points=sampling_points,max_iter=max_iter,grow_panels=grow_panels,grow_M=grow_M,verbose=verbose_cheb_panelization)
                n_panels_h,M_h,n_panels_j,M_j=out[1],out[2],out[3],out[4]
            end
            cache=build_ebim_cheb_cache(solver,pts,segks;n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j)
            @inbounds for (loc,i) in enumerate(seg_first:seg_last)
                results[i]=solve!(solver,A,dA,ddA,pts,ks[i],dks[i],cache,loc;use_lapack_raw=use_lapack_raw,multithreaded=multithreaded_matrices,use_krylov=use_krylov,nev=nevs[i],return_imag_part=return_imag_part)
                next!(p)
            end
        else
            @inbounds for i in seg_first:seg_last
                results[i]=solve!(solver,A,dA,ddA,pts,ks[i],dks[i];use_lapack_raw=use_lapack_raw,multithreaded=multithreaded_matrices,use_krylov=use_krylov,nev=nevs[i],return_imag_part=return_imag_part)
                next!(p)
            end
        end
        seg_first=seg_last+1
    end
    λs_all=K[]
    tensions_all=T[]
    control=Bool[]
    @inbounds for i in eachindex(results)
        λs,tens=results[i]
        isempty(λs)&&continue
        overlap_and_merge_ebim!(λs_all,tensions_all,λs,tens,control;tol=tol)
    end
    keep=[k1<=real(λ)<=k2 for λ in λs_all]
    return λs_all[keep],tensions_all[keep]
end

################################################################################
######################## DIAGNOSTIC UTILITY ####################################
################################################################################

function ebim_inv_diff(kvals::Vector{K},tens::Vector{C}) where {K<:Number,C<:Number}
    length(kvals)==length(tens)||throw(DimensionMismatch("kvals and tens must have equal length"))
    length(kvals)>=2||return K[],Float64[],C[]
    p=sortperm(real.(kvals))
    kvals=kvals[p]
    tens=tens[p]
    dr=diff(real.(kvals))
    invspacing=one(eltype(dr))./dr
    return kvals[1:end-1],invspacing,tens[1:end-1]
end