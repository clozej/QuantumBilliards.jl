#=
#############################################
########### BEYN CONTOUR METHOD #############
#############################################

MAIN REFERENCE: Beyn, Wolf-Jurgen, An integral method for solving nonlinear eigenvalue problems, 2012

- A Beyn-type contour-integral method for nonlinear eigenproblems T(k)φ=0 arising from BIM.
- On each disk Γ: center k0, radius R, we build A0 = (1/2πi)∮ T(z)^{-1}V dz,  A1 = (1/2πi)∮ z T(z)^{-1}V dz
then project with the rank-revealing SVD of A0 to a small dense B, and solve eigen!(B).
- Returned eigenpairs are filtered: (i) |k−k0|≤R and (ii) residual ‖T(k)φ‖ below a tolerance.

Practical guidance

- If you increase m (more levels per disk), increase both r (≥ m by a margin) and usually nq.
- Non-analytic boundaries converge slower with nq; expect to use higher nq and/or smaller R.
- Use solve_INFO logs to diagnose: Σ(A0) gaps, rank rk, counts kept/dropped, and residuals.
- Typical robust defaults: m≈100, Rmax≈0.5, nq≈40–70, r≈m+50, svd_tol≈1e-13, res_tol≈1e-10.
- For very high k or intricate geometries, start conservative (smaller R, larger nq) and relax if safe.

#TODO HSS with lu!
#TODO Backer's idea (Numerical details of wavefunction computation) of using the real Green's function Y0 or a combo with the beta param to avoid spurious ols associated with Y0 -> this would halve the contour nodes since Fredholm matrix would have conjugation symmetry. This is similar toe the FEAST algorithm / Zoloterov filter where we can halve the nodes of real symmetric matrices
#TODO For a smaller number of levels per contour we could dable with NLFEAST (Algorithm II): D. Kressner, Y. Liu, J. E. Roman, M. Shao, and N. Shao, "Linear convergence of iterative contour integral-based eigensolvers for nonlinear eigenvalue problems," arXiv:2606.13357 (2026).
=#

# when adding new ones just put them here and make sure they have construct construct_boundary_matrices! dispatch
const BeynSolver{T}=Union{BoundaryIntegralMethod{T},CFIE_kress{T},CFIE_kress_corners{T},CFIE_kress_global_corners{T},DLP_kress{T},DLP_kress_global_corners{T},CFIE_kress_composite_solver{T}}

#################
#### HELPERS ####
#################

# Width Δk containing approximately m levels from the leading Weyl estimate A*((k+Δk)^2-k^2)/(4π)=m.
@inline function weyl_window_width(billiard::Bi,k::T,m::Int;fundamental::Bool=true) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    A=fundamental ? fundamental_area(billiard) : area(billiard)
    return sqrt(k^2+T(4*pi*m/A))-k
end

# Cover [k1,k2] with consecutive Weyl-balanced windows. Each width is capped at 2Rmax so that the corresponding Beyn disk has radius R≤Rmax.
function plan_weyl_windows(billiard::Bi,k1::T,k2::T;m::Int=10,Rmax::Real=1.0,fundamental::Bool=true) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    k2>k1||return Tuple{T,T}[]
    m>0||throw(ArgumentError("m must be positive; received m=$m"))
    Rmax>0||throw(ArgumentError("Rmax must be positive; received Rmax=$Rmax"))
    iv=Tuple{T,T}[]
    maxwidth=T(2Rmax)
    k=k1
    while k<k2
        Δk=min(weyl_window_width(billiard,k,m;fundamental=fundamental),maxwidth,k2-k)
        Δk>zero(T)||throw(ArgumentError("Weyl window width vanished at k=$k"))
        kR=k+Δk
        push!(iv,(k,kR))
        k=kR
    end
    return iv
end

# Convert real windows [kL,kR] to circular Beyn contours with midpoint center and half-width radius.
function beyn_disks_from_windows(iv::Vector{Tuple{T,T}}) where {T<:Real}
    k0=Vector{Complex{T}}(undef,length(iv))
    R=Vector{T}(undef,length(iv))
    @inbounds for (i,(kL,kR)) in pairs(iv)
        k0[i]=complex((kL+kR)/2)
        R[i]=(kR-kL)/2
    end
    return k0,R
end

"""
    beyn_buffer_matrices(::Type{T},N::Int,r::Int,rng::G) where {T<:Real,G}

Allocate the random probing matrix and working matrices used to form the two Beyn contour moments

    A0=(1/2πi)∮ T(z)⁻¹V dz,
    A1=(1/2πi)∮ zT(z)⁻¹V dz.

The probing matrix is kept complex even for real-valued problems because the contour solves are generally complex.

## Arguments
- `T::Type{<:Real}`: Real scalar type used by the problem.
-- `N::Int`: Dimension of the boundary matrix.
-- `r::Int`: Number of random probing vectors.
- `rng::G`: Random-number generator used to construct `V`.

## Returns
A tuple `(V,X,A0,A1)` where:
- `V::Matrix{Complex{T}}` is the random `N×r` probing matrix.
- `X::Matrix{Complex{T}}` is the solve workspace for `T(z)X=V`.
- `A0::Matrix{Complex{T}}` stores the zeroth contour moment.
- `A1::Matrix{Complex{T}}` stores the first contour moment.
"""
function beyn_buffer_matrices(::Type{T},N::Int,r::Int,rng::G) where {T<:Real,G}
    V=randn(rng,Complex{T},N,r) # best leave as Complex even for Real problems to avoid issues in ldiv!
    X=similar(V)
    A0=zeros(Complex{T},N,r)
    A1=zeros(Complex{T},N,r)
    return V,X,A0,A1
end

"""
    construct_B_matrix(solver::BeynSolver{T},pts::Union{BoundaryPoints{T},Vector{BoundaryPoints{T}}},N::Int,k0::Complex{T},R::T;kwargs...) where {T<:Real}

Construct the reduced Beyn matrix for one circular contour centered at `k0` with radius `R`.
For trapezoidal contour nodes `zj` and weights `wj`, the method approximates

    A0=(1/2πi)∮ T(z)⁻¹V dz ≈ ∑ wj[j] T(zj[j])⁻¹V,
    A1=(1/2πi)∮ zT(z)⁻¹V dz ≈ ∑ wj[j] zj[j] T(zj[j])⁻¹V.

After the thin SVD

    A0=UΣW*,

the numerical rank `rk` is determined from the singular values satisfying `Σ[i]≥svd_tol`. With the retained factors `Uk`, `Wk`, and `Σk`, the reduced matrix is

    B=Uk'*A1*Wk*Σk⁻¹.

The eigenvalues of `B` approximate the nonlinear eigenvalues enclosed by the contour.

If all `r` singular values survive the cutoff, the probing dimension is increased and the moments are recomputed using the already factorized contour matrices.

## Arguments
- `solver::BeynSolver{T}`: Beyn-compatible boundary-integral solver.
- `pts::Union{BoundaryPoints{T},Vector{BoundaryPoints{T}}}`: Boundary discretization used to construct the nonlinear matrix `T(k)`.
- `N::Int`: Dimension of the boundary matrix.
- `k0::Complex{T}`: Center of the circular Beyn contour.
- `R::T`: Radius of the contour.

## Keyword Arguments
- `nq::Int=64`: Number of trapezoidal quadrature nodes on the contour.
- `r::Int=48`: Initial number of random probing vectors.
- `svd_tol::Real=1e-14`: Absolute singular-value threshold used for numerical rank detection.
- `rng=MersenneTwister(0)`: Random-number generator used to construct the probing matrix.
- `use_chebyshev::Bool=true`: Use Chebyshev-accelerated Hankel/Bessel evaluation in matrix construction.
- `n_panels_h::Int=15000`: Number of radial Chebyshev panels for Hankel functions.
- `M_h::Int=5`: Polynomial degree of the Hankel Chebyshev approximation.
- `n_panels_j::Int=10000`: Number of radial Chebyshev panels for Bessel-J functions.
- `M_j::Int=5`: Polynomial degree of the Bessel-J Chebyshev approximation.
- `info::Bool=false`: Print timing information.
- `multithreaded::Bool=true`: Enable multithreaded boundary-matrix construction.

## Returns
A tuple `(B,Uk)` where:
- `B::Matrix{Complex{T}}` is the reduced Beyn matrix.
- `Uk::AbstractMatrix{Complex{T}}` contains the retained left singular vectors of `A0`.

If numerical rank zero is detected, both returned matrices have zero columns.
"""
function construct_B_matrix(solver::BeynSolver{T},pts::Union{BoundaryPoints{T},Vector{BoundaryPoints{T}}},N::Int,k0::Complex{T},R::T;nq::Int=64,r::Int=48,svd_tol::Real=1e-14,rng=MersenneTwister(0),use_chebyshev::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,info::Bool=false,multithreaded::Bool=true) where {T<:Real}
    θ=range(zero(T),TWO_PI;length=nq+1)[1:end-1] # remove last point
    ej=cis.(θ) # unit circle points
    zj=k0.+R.*ej # contour points
    wj=(R/nq).*ej # contour weights
    #TODO Make the Fredholm matrices working buffers from outside to prevent large allocations in a loop (only for RAM critical applications since this is a small part of the actual execution time)
    Tbufs1=[zeros(Complex{T},N,N) for _ in 1:nq] 
    construct_boundary_matrices!(Tbufs1,solver,pts,zj;multithreaded=multithreaded,use_chebyshev=use_chebyshev,n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,timeit=info) # construct the T(zj) matrices for each contour point zj.
    # Allocate the buffers for the Beyn method. These are used in the matrix construction and then in the contour integrations to avoid repeated allocations. The matrices are sized according to the expected number of eigenvalues r and the size of the Fredholm matrices N.
    V,X,A0,A1=beyn_buffer_matrices(T,N,r,rng)
    # Now perform the Beyn contour integrations to form A0 and A1. To do this we need to solve T(zj) X = V for each zj and accumulate A0 += wj[j] * X, A1 += wj[j] * zj[j] * X. So as the first step we LU factor all T(zj) matrices to get the Fj factors which are used for ldiv! to efficiently solve the systems.
    @blas_multi MAX_BLAS_THREADS F1=lu!(Tbufs1[1];check=false) # just to get the type
    Fs=Vector{typeof(F1)}(undef,nq)
    Fs[1]=F1
    @benchit timeit=info "LU factorization" begin
        @blas_multi_then_1 MAX_BLAS_THREADS @inbounds for j in 2:nq # LU factor all T(zj) matrices
            Fs[j]=lu!(Tbufs1[j];check=false)
        end
    end
    xv=reshape(X,:);a0v=reshape(A0,:);a1v=reshape(A1,:) # vector views for BLAS.axpy! operations, to avoid allocations in the loop via reshaping the matrices each time in the loop
    @benchit timeit=info "Contour integration - ldiv! + axpy!" begin
        @blas_multi_then_1 MAX_BLAS_THREADS @inbounds for j in eachindex(zj)
            ldiv!(X,Fs[j],V) # make efficient inverse
            BLAS.axpy!(wj[j],xv,a0v) # A0 += wj[j] * X
            BLAS.axpy!(wj[j]*zj[j],xv,a1v) # A1 += wj[j] * zj[j] * X
        end
    end
    @blas_multi_then_1 MAX_BLAS_THREADS U,Σ,W=svd!(A0;full=false) # thin SVD of A0, revealing rank. The singular values > svd_tol correspond to eigenvalues. If all sv > svd_tol then maybe increase r (expected eigenvalue count) or reduce R (contour around k0), but if increasing r careful with nq. Check ref. section 3 eq. 22
    rk=count(>=(svd_tol),Σ) # filter out those that correspond to actual eigenvalues
    if rk==0 # if nothing found early return
        return Matrix{Complex{T}}(undef,0,0),Matrix{Complex{T}}(undef,N,0)
    end
    if rk==r
        r_tmp=min(r+r,N)
        while true # do again the ldiv + axpy accumulation with larger r until some sv < svd_tol. This does not require another Fredholm matrix construction since the same T(zj) can be used for larger r.
            V,X,A0,A1=beyn_buffer_matrices(T,N,r_tmp,rng)
            xv=reshape(X,:);a0v=reshape(A0,:);a1v=reshape(A1,:)
            @blas_multi_then_1 MAX_BLAS_THREADS @inbounds for j in eachindex(zj)
                ldiv!(X,Fs[j],V)
                BLAS.axpy!(wj[j],xv,a0v)
                BLAS.axpy!(wj[j]*zj[j],xv,a1v)
            end
            U,Σ,W=svd!(A0;full=false)
            rk=count(>=(svd_tol),Σ)
            rk<r_tmp&&break
            r_tmp==N&&throw(ArgumentError("Beyn moment remains rank-saturated at the maximum probe rank N=$N"))
            r_tmp=min(r_tmp+r,N)
        end
    end
    Uk=@view U[:,1:rk] # take the relevant ones corresponding to eigenvalues as in Integral algorithm 1 on p14 of ref
    Wk=@view W[:,1:rk] # take the relevant ones corresponding to eigenvalues as in Integral algorithm 1 on p14 of ref
    Σk=@view Σ[1:rk] # take the relevant ones corresponding to eigenvalues as in Integral algorithm 1 on p14 of ref
    # form B = adjoint(U) * A1 * W * Σ^{-1} as in the reference, p14, integral algorithm 1
    tmp=Matrix{Complex{T}}(undef,N,rk)
    @blas_multi_then_1 MAX_BLAS_THREADS mul!(tmp,A1,Wk)  # tmp := A1 * Wk, not weighted by inverse diagonal Σk
    @inbounds @simd for j in 1:rk # right-divide by diagonal Σk
        @views tmp[:,j]./=Σk[j]
    end
    B=Matrix{Complex{T}}(undef,rk,rk)
    @blas_multi_then_1 MAX_BLAS_THREADS mul!(B,adjoint(Uk),tmp) # B := Uk'*tmp, the final step
    return B,Uk
end

"""
    solve_vect(solver::BeynSolver{T},basis::Ba,pts::Union{BoundaryPoints{T},Vector{BoundaryPoints{T}}},k::Complex{T},dk::T;kwargs...) where {Ba<:AbstractHankelBasis,T<:Real}

Perform one Beyn contour solve on the circle centered at `k` with radius `dk`.
The reduced matrix is constructed with `construct_B_matrix` and diagonalized as

    B*Y=Y*Λ.

If `Uk` is the retained left singular-vector basis of the zeroth moment, the approximate boundary eigenvectors are represented by

    Φ=Uk*Y.

This function returns the provisional Beyn eigenvalues and the factors needed to reconstruct their boundary vectors. It does not perform contour-membership or residual filtering.

## Arguments
- `solver::BeynSolver{T}`: Beyn-compatible boundary-integral solver.
- `basis::Ba`: Hankel basis object used by the BIM interface.
- `pts::Union{BoundaryPoints{T},Vector{BoundaryPoints{T}}}`: Boundary discretization.
- `k::Complex{T}`: Center of the circular contour.
- `dk::T`: Radius of the circular contour.

## Keyword Arguments
- `multithreaded::Bool=true`: Enable multithreaded boundary-matrix construction.
- `nq::Int=32`: Number of contour quadrature nodes.
- `r::Int=48`: Initial random probe rank.
- `svd_tol::Real=1e-14`: Absolute singular-value threshold used for Beyn rank detection.
- `rng=MersenneTwister(0)`: Random-number generator used for the probing matrix.
- `use_chebyshev::Bool=true`: Enable Chebyshev-accelerated special-function evaluation.
- `n_panels_h::Int=15000`: Hankel Chebyshev panel count.
- `M_h::Int=5`: Hankel Chebyshev polynomial degree.
- `n_panels_j::Int=10000`: Bessel-J Chebyshev panel count.
- `M_j::Int=5`: Bessel-J Chebyshev polynomial degree.
- `info::Bool=false`: Print timing diagnostics.

## Returns
A tuple `(λ,Uk,Y)` where:
- `λ::Vector{Complex{T}}` contains the provisional complex eigenvalues.
- `Uk::AbstractMatrix{Complex{T}}` is the retained left singular-vector basis of the zeroth moment.
- `Y::Matrix{Complex{T}}` contains the eigenvectors of the reduced Beyn matrix.

The corresponding approximate boundary eigenvectors satisfy `Φ=Uk*Y`.
"""
function solve_vect(solver::BeynSolver{T},basis::Ba,pts::Union{BoundaryPoints{T},Vector{BoundaryPoints{T}}},k::Complex{T},dk::T;multithreaded::Bool=true,nq::Int=32,r::Int=48,svd_tol::Real=1e-14,rng=MersenneTwister(0),use_chebyshev::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,info::Bool=false) where {Ba<:AbstractHankelBasis} where {T<:Real}
    N=boundary_matrix_size(solver,pts) # Get the boundary-matrix dimension for the active solver and discretization.
    B,Uk=construct_B_matrix(solver,pts,N,k,dk,nq=nq,r=r,svd_tol=svd_tol,rng=rng,use_chebyshev=use_chebyshev,n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,multithreaded=multithreaded,info=info) # here is where the core of the algorithm is found. Constructs B from step 5 in ref p.14
    if isempty(B) # rk==0
        @info "no_roots_in_window" k0=k R=dk nq=nq svd_tol=svd_tol
        return Complex{T}[],Uk,Matrix{Complex{T}}(undef,0,0)
    end
    @blas_multi_then_1 MAX_BLAS_THREADS λ,Y=eigen!(B) # small dense eigendecomposition to get eigenvalues λ are the eigenvalues and v(λ) are the eigenvectors
    # Now form only relevant cols of Φ = U * Y since A0 = U Σ W*, we have A0 * W Σ^{-1} Y = U Y. Each column is now an eigenvector of of T(λ)v(λ) = 0. This is the second layer potential boundary operator now!
    #println("Eigenvalues found in window k0=$(k), R=$(dk): ",λ)
    return λ,Uk,Y
end

"""
    solve(solver::BeynSolver{T},basis::Ba,pts::Union{BoundaryPoints{T},Vector{BoundaryPoints{T}}},k::Complex{T},dk::T;kwargs...) where {Ba<:AbstractHankelBasis,T<:Real}

Return only the provisional Beyn eigenvalues obtained from one circular contour.
This is a lightweight interface to `construct_B_matrix` followed by the eigendecomposition of the reduced Beyn matrix. It does not perform contour-membership or residual filtering.

## Arguments
- `solver::BeynSolver{T}`: Beyn-compatible boundary-integral solver.
- `basis::Ba`: Hankel basis object used by the BIM interface.
- `pts::Union{BoundaryPoints{T},Vector{BoundaryPoints{T}}}`: Boundary discretization.
- `k::Complex{T}`: Center of the circular contour.
- `dk::T`: Radius of the contour.

## Keyword Arguments
- `multithreaded::Bool=true`: Enable multithreaded boundary-matrix construction.
- `nq::Int=32`: Number of contour quadrature nodes.
- `r::Int=48`: Initial random probe rank.
- `svd_tol::Real=1e-14`: Absolute singular-value threshold used for Beyn rank detection.
- `res_tol::Real=1e-8`: Compatibility keyword; residual filtering is not performed here.
- `rng=MersenneTwister(0)`: Random-number generator used for the probing matrix.
- `auto_discard_spurious::Bool=true`: Compatibility keyword; spurious-root filtering is not performed here.
- `use_chebyshev::Bool=true`: Enable Chebyshev-accelerated special-function evaluation.
- `n_panels_h::Int=15000`: Hankel Chebyshev panel count.
- `M_h::Int=5`: Hankel Chebyshev polynomial degree.
- `n_panels_j::Int=10000`: Bessel-J Chebyshev panel count.
- `M_j::Int=5`: Bessel-J Chebyshev polynomial degree.
- `info::Bool=false`: Print timing diagnostics.

## Returns
`Vector{Complex{T}}` containing the provisional Beyn eigenvalues.

## Notes
Use `solve_vect` together with `residual_and_norm_select`, or use `compute_spectrum_beyn`, when validated eigenpairs are required.
"""
function solve(solver::BeynSolver{T},basis::Ba,pts::Union{BoundaryPoints{T},Vector{BoundaryPoints{T}}},k::Complex{T},dk::T;multithreaded::Bool=true,nq::Int=32,r::Int=48,svd_tol::Real=1e-14,res_tol::Real=1e-8,rng=MersenneTwister(0),auto_discard_spurious::Bool=true,use_chebyshev::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,info::Bool=false) where {Ba<:AbstractHankelBasis} where {T<:Real}
    N=boundary_matrix_size(solver,pts) # get the size of the boundary matrix based on the type of pts (BoundaryPoints or Vector{BoundaryPointsCFIE})
    B,_=construct_B_matrix(solver,pts,N,k,dk,nq=nq,r=r,svd_tol=svd_tol,rng=rng,use_chebyshev=use_chebyshev,n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,multithreaded=multithreaded,info=info) # here is where the core of the algorithm is found. Constructs B from step 5 in ref p.14
    if isempty(B) # rk==0
        @info "no_roots_in_window" k0=k R=dk nq=nq svd_tol=svd_tol
        return Complex{T}[]
    end
    @blas_multi_then_1 MAX_BLAS_THREADS λ=eigvals!(B) # small dense eigendecomposition to get eigenvalues 
    return λ
end

"""
    solve_INFO(solver::BeynSolver{T},basis::Ba,pts::Union{BoundaryPoints{T},Vector{BoundaryPoints{T}}},k0::Complex{T},R::T;kwargs...) where {Ba<:AbstractHankelBasis,T<:Real}

Run one diagnostic Beyn solve and report the internal numerical behavior of the contour method.
The routine constructs and factorizes all contour matrices, forms the Beyn moments, prints the singular values of `A0`, determines the numerical rank, solves the reduced eigenproblem, reconstructs

    Φ=Uk*Y,

and evaluates the raw nonlinear residual

    ||T(λj)φj||

for each candidate inside the contour.

This routine is intended for choosing suitable values of `nq`, `r`, and `svd_tol` before a production spectral sweep.

## Arguments
- `solver::BeynSolver{T}`: Beyn-compatible boundary-integral solver.
- `basis::Ba`: Hankel basis object used by the BIM interface.
- `pts::Union{BoundaryPoints{T},Vector{BoundaryPoints{T}}}`: Boundary discretization.
- `k0::Complex{T}`: Center of the circular contour.
- `R::T`: Radius of the contour.

## Keyword Arguments
- `multithreaded::Bool=true`: Enable multithreaded boundary-matrix construction.
- `nq::Int=48`: Number of contour quadrature nodes.
- `r::Int=48`: Number of random probing vectors.
- `svd_tol::Real=1e-10`: Singular-value threshold used for numerical rank detection.
- `res_tol::Real=1e-10`: Raw residual threshold used when rejecting candidates.
- `rng=MersenneTwister(0)`: Random-number generator used for the probing matrix.
- `use_adaptive_svd_tol::Bool=false`: Replace `svd_tol` by `maximum(Σ)*1e-15`.
- `auto_discard_spurious::Bool=false`: Discard contour-inside candidates whose residual is at least `res_tol`.
- `use_chebyshev::Bool=true`: Enable Chebyshev-accelerated special-function evaluation.
- `n_panels_h::Int=15000`: Hankel Chebyshev panel count.
- `M_h::Int=5`: Hankel Chebyshev polynomial degree.
- `n_panels_j::Int=10000`: Bessel-J Chebyshev panel count.
- `M_j::Int=5`: Bessel-J Chebyshev polynomial degree.

## Returns
A tuple `(λ,Phi,tens)` where:
- `λ::Vector{Complex{T}}` contains the retained contour-inside eigenvalues.
- `Phi::Matrix{Complex{T}}` contains their reconstructed boundary eigenvectors.
- `tens::Vector{T}` contains the corresponding raw residual norms `||T(λ)φ||`.
"""
function solve_INFO(solver::BeynSolver{T},basis::Ba,pts::Union{BoundaryPoints{T},Vector{BoundaryPoints{T}}},k0::Complex{T},R::T;multithreaded::Bool=true,nq::Int=48,r::Int=48,svd_tol::Real=1e-10,res_tol::Real=1e-10,rng=MersenneTwister(0),use_adaptive_svd_tol=false,auto_discard_spurious=false,use_chebyshev=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5) where {Ba<:AbstractHankelBasis,T<:Real}
    N=boundary_matrix_size(solver,pts) # get the size of the boundary matrix based on the type of pts (BoundaryPoints or Vector{BoundaryPointsCFIE})
    θ=range(zero(T),TWO_PI;length=nq+1);θ=θ[1:end-1];ej=cis.(θ);zj=k0.+R.*ej;wj=(R/nq).*ej # contour points and weights
    V,X,A0,A1=beyn_buffer_matrices(T,N,r,rng)
    @info "beyn:start" k0=k0 R=R nq=nq N=N r=r
    Tbufs1=[zeros(Complex{T},N,N) for _ in 1:nq] 
    construct_boundary_matrices!(Tbufs1,solver,pts,zj;multithreaded=multithreaded,use_chebyshev=use_chebyshev,n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,timeit=true) # construct the T(zj) matrices for each contour point zj.
    @blas_multi MAX_BLAS_THREADS F1=lu!(Tbufs1[1];check=false) # just to get the type
    Fs=Vector{typeof(F1)}(undef,nq)
    Fs[1]=F1
    @blas_multi_then_1 MAX_BLAS_THREADS @inbounds begin
        @showprogress desc="lu!" for j in 2:nq # LU factor all T(zj) matrices
            Fs[j]=lu!(Tbufs1[j];check=false)
        end
    end
    xv=reshape(X,:);a0v=reshape(A0,:);a1v=reshape(A1,:) # vector views for BLAS.axpy! operations, to avoid allocations in the loop via reshaping the matrices each time in the loop
    @blas_multi_then_1 MAX_BLAS_THREADS @inbounds begin
        @showprogress desc="ldiv! + axpy!" for j in eachindex(zj)
            ldiv!(X,Fs[j],V) # make efficient inverse
            BLAS.axpy!(wj[j],xv,a0v) # A0 += wj[j] * X
            BLAS.axpy!(wj[j]*zj[j],xv,a1v) # A1 += wj[j] * zj[j] * X
        end
    end
    @time "SVD" @blas_multi_then_1 MAX_BLAS_THREADS U,Σ,W=svd!(A0;full=false)
    println("Singular values (<1e-10 tail inspection): ",Σ)
    rk=0
    svd_tol=use_adaptive_svd_tol ? maximum(Σ)*1e-15 : svd_tol
    @inbounds for i in eachindex(Σ)
        if Σ[i]≥svd_tol
            rk+=1
        else
            break
        end
    end
    rk==r && @warn "All singular values are above svd_tol = $(svd_tol), r = $(r) needs to be increased" # in the actual implementation where B matrix is constructed this will increase r by a fixed amount and do the procedure again until we have some singular values under tolerance!
    rk==0 && return Complex{T}[],Matrix{Complex{T}}(undef,N,0),T[]
    Uk=@view U[:,1:rk]
    Wk=@view W[:,1:rk]
    Σk=@view Σ[1:rk]
    tmp=Matrix{Complex{T}}(undef,N,rk)
    @blas_multi MAX_BLAS_THREADS mul!(tmp,A1,Wk)
    @inbounds for j in 1:rk
        @views tmp[:,j]./=Σk[j]
    end
    B=Matrix{Complex{T}}(undef,rk,rk)
    @blas_multi MAX_BLAS_THREADS mul!(B,adjoint(Uk),tmp)
    @time "eigen" @blas_multi_then_1 MAX_BLAS_THREADS ev=eigen!(B)
    λ=ev.values;Y=ev.vectors;Phi=Uk*Y
    keep=trues(length(λ))
    tens=Vector{T}()
    ybuf=Vector{Complex{T}}(undef,size(Phi,1))
    dropped_out=0
    dropped_res=0
    res_keep=T[]
    Tbuf_check=[zeros(Complex{T},N,N)]
    begin 
        @inbounds for j in eachindex(λ)
            d=abs(λ[j]-k0)
            if d>R
                keep[j]=false
                dropped_out+=1
                continue
            end
            fill!(Tbuf_check[1],0.0+0.0im)
            construct_boundary_matrices!(Tbuf_check,solver,pts,[λ[j]];multithreaded=multithreaded,use_chebyshev=use_chebyshev,n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,timeit=false) # construct the T(λ[j]) matrix for the eigenvalue λ[j] to check the residual
            @blas_multi_then_1 MAX_BLAS_THREADS mul!(ybuf,Tbuf_check[1],@view(Phi[:,j]))
            ybuf_norm=norm(ybuf)
            @info "k=$((λ[j])) ||A(k)v(k)|| = $(ybuf_norm)"
            if auto_discard_spurious
                if ybuf_norm≥res_tol
                    keep[j]=false
                    dropped_res+=1
                    if ybuf_norm>1e-8
                        if ybuf_norm>1e-6 # heuristic for when usually it is spurious sqrt(eps())
                            @warn "k=$((λ[j])) ||A(k)v(k)|| = $(ybuf_norm) > $res_tol , definitely spurious" 
                        else # gray zone
                            @warn "k=$((λ[j])) ||A(k)v(k)|| = $(ybuf_norm) > $res_tol , most probably eigenvalue but too low nq" 
                        end
                    else
                        @warn "k=$((λ[j])) ||A(k)v(k)|| = $(ybuf_norm) > $res_tol , could be spurious or try increasing nq (usually spurious) or lowering residual tolerance" 
                    end
                    continue
                end
            end
            push!(tens,ybuf_norm)
            push!(res_keep,ybuf_norm)
        end
        kept=count(keep)
        if kept>0
            @info "STATUS: " kept=kept dropped_outside=dropped_out dropped_residual=dropped_res max_residual=maximum(res_keep)
        else
            @info "STATUS: " kept=0 dropped_outside=dropped_out dropped_residual=dropped_res
        end
    end
    return λ[keep],Phi[:,keep],tens
end

"""
    residual_and_norm_select(solver::BeynSolver{T},λ::AbstractVector{Complex{T}},Uk::AbstractMatrix{Complex{T}},Y::AbstractMatrix{Complex{T}},k0::Complex{T},R::T,pts::Union{BoundaryPoints{T},Vector{BoundaryPoints{T}}};kwargs...) where {T<:Real}

Validate provisional Beyn eigenpairs by contour membership and nonlinear residual.
For each candidate `λ[j]` inside the contour, the corresponding boundary vector is reconstructed as

    φj=Uk*Y[:,j],

and the raw residual

    rj=||T(λj)φj||

is evaluated. A scale-independent normalized residual is also computed as

    rj_norm=||T(λj)φj||/(||T(λj)||*||φj||),

using the norm selected by `matnorm`.

Candidates outside the contour are discarded. If `auto_discard_spurious=true`, candidates with raw residual `rj≥res_tol` are also discarded.

## Arguments
- `solver::BeynSolver{T}`: Beyn-compatible boundary-integral solver.
- `λ::AbstractVector{Complex{T}}`: Provisional complex Beyn eigenvalues.
- `Uk::AbstractMatrix{Complex{T}}`: Retained left singular vectors of the zeroth Beyn moment.
- `Y::AbstractMatrix{Complex{T}}`: Eigenvectors of the reduced Beyn matrix.
- `k0::Complex{T}`: Center of the Beyn contour.
- `R::T`: Radius of the contour.
- `pts::Union{BoundaryPoints{T},Vector{BoundaryPoints{T}}}`: Boundary discretization associated with the solve.

## Keyword Arguments
- `res_tol::T`: Raw residual threshold used for candidate rejection.
- `matnorm::Symbol=:one`: Norm used for residual normalization; one of `:one`, `:two`, or `:inf`.
- `epss::Real=1e-15`: Small positive number protecting the normalization against division by zero.
- `auto_discard_spurious::Bool=true`: Reject candidates whose raw residual is at least `res_tol`.
- `collect_logs::Bool=false`: Collect textual keep/drop diagnostics.
- `use_chebyshev::Bool=true`: Enable Chebyshev-accelerated special-function evaluation.
- `n_panels_h::Int=15000`: Hankel Chebyshev panel count.
- `M_h::Int=5`: Hankel Chebyshev polynomial degree.
- `n_panels_j::Int=10000`: Bessel-J Chebyshev panel count.
- `M_j::Int=5`: Bessel-J Chebyshev polynomial degree.
- `multithreaded::Bool=true`: Enable multithreaded boundary-matrix construction.

## Returns
A tuple `(idx,Φ_kept,tens,tensN,logs)` where:
- `idx::Vector{Int}` contains the retained local eigenvalue indices.
- `Φ_kept::Matrix{Complex{T}}` contains the retained reconstructed boundary eigenvectors.
- `tens::Vector{T}` contains the raw residual norms.
- `tensN::Vector{T}` contains the normalized residuals.
- `logs::Vector{String}` contains optional selection diagnostics.
"""
function residual_and_norm_select(solver::BeynSolver{T},λ::AbstractVector{Complex{T}},Uk::AbstractMatrix{Complex{T}},Y::AbstractMatrix{Complex{T}},k0::Complex{T},R::T,pts::Union{BoundaryPoints{T},Vector{BoundaryPoints{T}}};res_tol::T,matnorm::Symbol=:one,epss::Real=1e-15,auto_discard_spurious::Bool=true,collect_logs::Bool=false,use_chebyshev::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,multithreaded::Bool=true) where {T<:Real}
    N,rk=size(Uk)
    Φtmp=Matrix{Complex{T}}(undef,N,rk)
    y=Vector{Complex{T}}(undef,N)
    keep=falses(rk)
    tens=Vector{T}(undef,rk)
    tensN=Vector{T}(undef,rk)
    logs=collect_logs ? String[] : nothing
    Tbufs=[zeros(Complex{T},N,N) for _ in eachindex(λ)]
    construct_boundary_matrices!(Tbufs,solver,pts,λ;multithreaded=multithreaded,use_chebyshev=use_chebyshev,n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,timeit=false)
    vecnorm= matnorm===:one ? (v->norm(v,1)) : matnorm===:two ? (v->norm(v)) : (v->norm(v,Inf))
    @inbounds for j in 1:rk
        λj=λ[j]
        if abs(λj-k0)>R
            tens[j]=T(NaN)
            tensN[j]=T(NaN)
            continue
        end
        @blas_multi_then_1 MAX_BLAS_THREADS mul!(@view(Φtmp[:,j]),Uk,@view(Y[:,j]))
        @blas_multi_then_1 MAX_BLAS_THREADS mul!(y,Tbufs[j],@view(Φtmp[:,j]))
        rj=norm(y)
        tens[j]=rj
        nA= matnorm===:one ? opnorm(Tbufs[j],1) : matnorm===:two ? opnorm(Tbufs[j],2) : opnorm(Tbufs[j],Inf)
        φn=vecnorm(@view(Φtmp[:,j]))
        yn=vecnorm(y)
        tensN[j]=yn/(nA*(φn+epss)+epss)
        if auto_discard_spurious && rj>=res_tol
            collect_logs && push!(logs,"λ=$(λj) ||Aφ||=$(rj) > $res_tol → DROP")
        else
            keep[j]=true
            collect_logs && push!(logs,"λ=$(λj) ||Aφ||=$(rj) < $res_tol ← KEEP")
        end
    end
    idx=findall(keep)
    Φ_kept=isempty(idx) ? Matrix{Complex{T}}(undef,N,0) : Φtmp[:,idx]
    return idx,Φ_kept,tens[idx],tensN[idx],(collect_logs ? logs : String[])
end

"""
    imag_k_check(solver::BeynSolver{T},λs::Vector{Vector{Complex{T}}},Uks::Vector{Matrix{Complex{T}}},Ys::Vector{Matrix{Complex{T}}},k0s::Vector{Complex{T}},Rs::Vector{T},all_pts;kwargs...) where {T<:Real}

Apply an experimental fast spurious-root filter based on the imaginary parts of provisional Beyn eigenvalues.
For problems with a theoretically real spectrum, discretization and contour-integration errors generally move approximate eigenvalues away from the real axis. The candidates from all Beyn windows are therefore sorted by decreasing `|Im λ|`, and residual checks are performed first on this suspicious tail.

A checked candidate is discarded when

    ||T(λ)φ||≥res_tol.

The scan terminates after `pad` consecutive residual-good candidates, after which the remaining smaller-`|Im λ|` candidates are accepted without explicit residual evaluation.

Residuals are always evaluated using the boundary discretization of the Beyn window that generated the corresponding eigenpair.

## Arguments
- `solver::BeynSolver{T}`: Beyn-compatible boundary-integral solver.
- `λs::Vector{Vector{Complex{T}}}`: Provisional eigenvalues for each Beyn window.
- `Uks::Vector{Matrix{Complex{T}}}`: Retained left singular-vector bases for each window.
- `Ys::Vector{Matrix{Complex{T}}}`: Reduced Beyn eigenvectors for each window.
- `k0s::Vector{Complex{T}}`: Beyn contour centers.
- `Rs::Vector{T}`: Beyn contour radii.
- `all_pts`: Boundary discretization associated with each window.

## Keyword Arguments
- `res_tol::T`: Raw residual threshold used to reject checked candidates.
- `pad::Int=20`: Number of consecutive residual-good candidates required before terminating the scan.
- `group_size::Int=100`: Maximum number of suspicious candidates processed in one batch.
- `use_chebyshev::Bool=true`: Enable Chebyshev-accelerated special-function evaluation.
- `n_panels_h::Int=15000`: Hankel Chebyshev panel count.
- `M_h::Int=5`: Hankel Chebyshev polynomial degree.
- `n_panels_j::Int=10000`: Bessel-J Chebyshev panel count.
- `M_j::Int=5`: Bessel-J Chebyshev polynomial degree.
- `multithreaded::Bool=true`: Enable multithreaded boundary-matrix construction.
- `verbose::Bool=true`: Print candidate and filtering diagnostics.

## Returns
A tuple `(idx_keep,residuals)` where:
- `idx_keep::Vector{Vector{Int}}` contains the surviving local eigenvalue indices for each window.
- `residuals::Vector{Vector{T}}` contains residuals aligned with the surviving candidates.

An accepted candidate that was not explicitly residual-checked has residual `NaN`.

## Notes
This is a heuristic acceleration intended for problems whose exact spectrum is real. It should be validated against complete residual checking when changing solver type, geometry, discretization, contour quadrature, or spectral regime.
"""
function imag_k_check(solver::BeynSolver{T},λs::Vector{Vector{Complex{T}}},Uks::Vector{Matrix{Complex{T}}},Ys::Vector{Matrix{Complex{T}}},k0s::Vector{Complex{T}},Rs::Vector{T},all_pts;res_tol::T,pad::Int=20,group_size::Int=100,use_chebyshev::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,multithreaded::Bool=true,verbose::Bool=true) where {T<:Real}
    nw=length(λs)
    idx_inside=Vector{Vector{Int}}(undef,nw)
    idx_keep=Vector{Vector{Int}}(undef,nw)
    residuals=Vector{Vector{T}}(undef,nw)
    local_pos=Dict{Tuple{Int,Int},Int}()
    candidates=Tuple{Int,Int,T,T}[] # (window index, local eigenvalue index, |Im λ|, Re λ)
    # 1. Collect all contour-inside roots and rank by suspiciousness.
    @time "collect candidates and sort by |Im λ|" begin
        @inbounds for i in 1:nw
            λi=λs[i]
            idx_inside[i]=findall(j->abs(λi[j]-k0s[i])<=Rs[i],eachindex(λi))
            residuals[i]=fill(T(NaN),length(idx_inside[i]))
            for (lp,j) in pairs(idx_inside[i])
                local_pos[(i,j)]=lp
                push!(candidates,(i,j,abs(imag(λi[j])),real(λi[j])))
            end
        end
    end
    # Sort by |Im λ| descending so the most suspicious candidates are checked first.
    sort!(candidates;by=c->c[3],rev=true)
    verbose && begin
        @info "top imag candidates" first_candidates=candidates[1:min(10,length(candidates))]
    end
    drop=Dict{Tuple{Int,Int},Bool}()
    checked=0
    dropped=0
    good_streak=0
    stop_early=false
    pos=1
    # 2. Check only the suspicious large-|Im λ| tail in batches.
    while pos<=length(candidates) && !stop_early
        stop=min(pos+group_size-1,length(candidates))
        # Suspicious roots, globally sorted by descending |Im λ|.
        group=candidates[pos:stop]
        # Residual lookup keyed by (window index, local eigenvalue index).
        rdict=Dict{Tuple{Int,Int},T}()
        # Residuals must be evaluated using the same discretization that produced
        # the corresponding Beyn eigenpair, so split mixed batches by window.
        @time "residual check for candidates with |Im λ| in [$(group[end][3]), $(group[1][3])]" begin
            for iwin in unique(c[1] for c in group)
                sub=[c for c in group if c[1]==iwin]
                pts_ref=all_pts[iwin]
                Nref=boundary_matrix_size(solver,pts_ref)
                # Batched A(λ) assembly for suspicious roots from this window.
                λ_group=Complex{T}[λs[i][j] for (i,j,_,_) in sub]
                Tbufs=[zeros(Complex{T},Nref,Nref) for _ in eachindex(λ_group)]
                construct_boundary_matrices!(Tbufs,solver,pts_ref,λ_group;multithreaded=multithreaded,use_chebyshev=use_chebyshev,n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,timeit=false)
                ybuf=Vector{Complex{T}}(undef,Nref)
                φbuf=Vector{Complex{T}}(undef,Nref)
                @blas_multi_then_1 MAX_BLAS_THREADS begin 
                    @inbounds for q in eachindex(sub)
                        i,j,_,_=sub[q]
                        # Reconstruct Beyn eigenvector φ = Uk*Y[:,j].
                        mul!(φbuf,Uks[i],@view(Ys[i][:,j]))
                        # True residual ||A(λ)φ||.
                        mul!(ybuf,Tbufs[q],φbuf)
                        rdict[(i,j)]=norm(ybuf)
                    end
                end
            end
        end
        @time "filter candidates with |Im λ| in [$(group[end][3]), $(group[1][3])]" begin
            # Consume in original suspiciousness order.
            @inbounds for c in group
                i,j,imj,_=c
                rj=rdict[(i,j)]
                checked+=1
                residuals[i][local_pos[(i,j)]]=rj
                if rj>=res_tol
                    verbose && @info "DROP candidate" i=i j=j k=λs[i][j] abs_imag=imj residual=rj
                    drop[(i,j)]=true
                    dropped+=1
                    good_streak=0
                else
                    good_streak+=1
                    # Once enough consecutive good roots appear, assume the smaller-|Im λ| tail is clean.
                    if good_streak>=pad
                        verbose && @info "grouped tail residual check stopped" checked=checked dropped=dropped good_streak=good_streak last_imag=imj last_residual=rj
                        stop_early=true
                        break
                    end
                end
            end
            pos=stop+1
        end
    end
    # 3. Apply drops and keep residuals aligned with surviving roots.
    @inbounds for i in 1:nw
        old=idx_inside[i]
        mask=[!get(drop,(i,j),false) for j in old]
        idx_keep[i]=old[mask]
        residuals[i]=residuals[i][mask]
    end
    verbose && @info "grouped tail residual check summary" checked=checked dropped=dropped total_candidates=length(candidates)
    return idx_keep,residuals
end

########################
#### HIGH LEVEL API ####
########################

"""
    solve_spectrum_beyn(solver::BeynSolver{T},billiard::Bi,k1::T,k2::T;kwargs...) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}

Compute the boundary-integral spectrum in `[k1,k2]` using Weyl-balanced circular Beyn contours.
The requested interval is first partitioned into windows containing approximately `m` levels according to the leading Weyl estimate

    N(k)≈A*k^2/(4π).

Each window is converted to a circular contour with radius at most `Rmax`. For every contour, the algorithm constructs the boundary discretization, forms the Beyn moments, solves the reduced nonlinear eigenproblem, and reconstructs the corresponding boundary eigenvectors. Candidate eigenvalues are subsequently filtered either by complete residual evaluation or by the experimental `|Im λ|` prefilter. When `use_chebyshev=true`, the Hankel/Bessel interpolation parameters are tuned once using the highest-frequency contour and reused for the complete spectral sweep.

## Arguments
- `solver::BeynSolver{T}`: Beyn-compatible boundary-integral solver.
- `billiard::Bi`: Billiard geometry.
- `k1::T`: Lower endpoint of the requested wavenumber interval.
- `k2::T`: Upper endpoint of the requested wavenumber interval.

## Keyword Arguments
- `m::Int=50`: Target number of eigenvalues per Weyl window.
- `Rmax::T=T(0.5)`: Maximum radius of each Beyn contour.
- `nq::Int=48`: Number of trapezoidal quadrature nodes on each contour.
- `r::Int=m+15`: Initial number of random probing vectors.
- `svd_tol::Real=1e-12`: Absolute singular-value threshold used for Beyn rank detection.
- `res_tol::Real=1e-9`: Raw nonlinear residual threshold used for candidate rejection.
- `auto_discard_spurious::Bool=true`: Enable residual-based rejection of spurious candidates.
- `multithreaded_matrix::Bool=true`: Enable multithreaded boundary-matrix construction.
- `use_adaptive_svd_tol::Bool=false`: Use an adaptive SVD threshold in the initial diagnostic solve.
- `use_chebyshev::Bool=true`: Enable Chebyshev interpolation of Hankel/Bessel functions.
- `n_panels_h::Int=15000`: Initial Hankel Chebyshev panel count.
- `M_h::Int=5`: Initial Hankel Chebyshev polynomial degree.
- `n_panels_j::Int=10000`: Initial Bessel-J Chebyshev panel count.
- `M_j::Int=5`: Initial Bessel-J Chebyshev polynomial degree.
- `do_INFO_init::Bool=true`: Run `solve_INFO` on a representative contour before the production sweep.
- `do_per_solve_INFO::Bool=true`: Enable detailed timing and diagnostic output during the sweep.
- `cheb_tol::Real=1e-13`: Target absolute error for Chebyshev parameter tuning.
- `max_iter::Int=20`: Maximum number of Chebyshev tuning iterations.
- `sampling_points::Int=50_000`: Number of radial validation points used during Chebyshev tuning.
- `grow_panels::Real=1.5`: Multiplicative growth factor for Chebyshev panel counts.
- `grow_M::Int=2`: Additive growth of the Chebyshev polynomial degree.
- `return_imag_part::Bool=false`: Return complex approximate eigenvalues instead of only their real parts.
- `use_imag_check_EXPERIMENTAL::Bool=true`: Use the `|Im λ|` spurious-root filter instead of residual-checking every candidate. Recommended for production when the exact spectrum is known to be real and the heuristic has been validated for the problem class.

## Returns
A tuple `(ks,tens,us,pts,tensN)` where:
- `ks::Union{Vector{T},Vector{Complex{T}}}` contains the retained wavenumbers.
- `tens::Vector{T}` contains raw residuals with complete residual checking, or `|Im λ|` with the experimental filter.
- `us::Vector{Vector{Complex{T}}}` contains the retained boundary densities/eigenvectors.
- `pts::Vector{pts_type}` contains the boundary discretization associated with each returned eigenpair.
- `tensN::Vector{T}` contains normalized residuals with complete residual checking, or explicitly evaluated raw residuals with the experimental filter.

## Notes
- Increasing `m` generally requires increasing both `r` and `nq`.
- Increasing `nq` is the primary way to reduce contour-quadrature error when genuine roots have poor residuals.
- Decreasing `Rmax` can improve robustness when `T(k)` varies rapidly (high k sometimes).
"""
function solve_spectrum_beyn(solver::BeynSolver{T},billiard::Bi,k1::T,k2::T;m::Int=50,Rmax::T=T(0.5),nq::Int=48,r::Int=m+15,svd_tol::Real=1e-12,res_tol::Real=1e-9,auto_discard_spurious::Bool=true,multithreaded_matrix::Bool=true,use_adaptive_svd_tol::Bool=false,use_chebyshev::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,do_INFO_init::Bool=true,do_per_solve_INFO::Bool=true,cheb_tol::Real=1e-13,max_iter::Int=20,sampling_points::Int=50_000,grow_panels::Real=1.5,grow_M::Int=2,return_imag_part::Bool=false,use_imag_check_EXPERIMENTAL::Bool=true) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    fundamental=!isnothing(solver.symmetry)
    basis=AbstractHankelBasis()
    intervals=plan_weyl_windows(billiard,k1,k2;m=m,fundamental=fundamental,Rmax=Rmax)
    if length(intervals)>=2
        kL2,kR2=intervals[end-1]
        kL3,kR3=intervals[end]
        len3=kR3-kL3
        if len3<=max(100*eps(k2),1e-9*max(k2,one(T)))
            if (kR3-kL2)<=2*Rmax+10*eps(k2)
                intervals[end-1]=(kL2,kR3)
                pop!(intervals)
            else
                pop!(intervals)
            end
        end
    end
    k0,R=beyn_disks_from_windows(intervals)
    # get the type of pts
    pts_type=typeof(evaluate_points(solver,billiard,k1))
    isempty(k0)&&return (return_imag_part ? Complex{T}[] : T[]),T[],Vector{Vector{Complex{T}}}(),Vector{pts_type}(),T[]
    do_INFO_init && @info "Weyl windows planned" intervals=intervals k0=k0 R=R
    all_pts=Vector{pts_type}(undef,length(k0))
    p=Progress(length(k0),desc="Computing points for each disk...")
    @benchit timeit=do_per_solve_INFO "Point evaluation" begin
        @use_threads multithreading=true for i in eachindex(k0)
            all_pts[i]=evaluate_points(solver,billiard,real(k0[i]))
            next!(p)
        end
    end
    λs=Vector{Vector{Complex{T}}}(undef,length(k0))
    Uks=Vector{Matrix{Complex{T}}}(undef,length(k0))
    Ys=Vector{Matrix{Complex{T}}}(undef,length(k0))
    if use_chebyshev
        imax=argmax(real.(k0).+R)
        θref=range(zero(T),TWO_PI;length=nq+1)[1:end-1]
        zj_ref=k0[imax].+R[imax].*cis.(θref)
        cheb_out=chebyshev_params(solver,all_pts[imax],zj_ref;tol=cheb_tol,npanels_h_init=n_panels_h,M_h_init=M_h,npanels_j_init=n_panels_j,M_j_init=M_j,sampling_points=sampling_points,max_iter=max_iter,grow_panels=grow_panels,grow_M=grow_M,verbose=do_per_solve_INFO)
        n_panels_h,M_h,n_panels_j,M_j=cheb_out[1],cheb_out[2],cheb_out[3],cheb_out[4]
    end
    if do_INFO_init
        mid=cld(length(k0),2)
        @benchit timeit=do_INFO_init "solve_INFO representative disk" begin
            _=solve_INFO(solver,basis,all_pts[mid],complex(k0[mid]),R[mid];nq=nq,r=r,svd_tol=svd_tol,res_tol=res_tol,use_adaptive_svd_tol=use_adaptive_svd_tol,multithreaded=multithreaded_matrix,use_chebyshev=use_chebyshev,n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j)
        end
    end
    p=Progress(length(k0),1)
    @time "Beyn pass (all disks)" begin
        @inbounds for i in eachindex(k0)
            λ,Uk,Y=solve_vect(solver,basis,all_pts[i],complex(k0[i]),R[i];nq=nq,r=r,svd_tol=svd_tol,rng=Random.MersenneTwister(0),multithreaded=multithreaded_matrix,use_chebyshev=use_chebyshev,n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,info=do_per_solve_INFO)
            λs[i]=λ
            Uks[i]=Uk
            Ys[i]=Y
            next!(p)
        end
    end
    ks_list=return_imag_part ? Vector{Vector{Complex{T}}}(undef,length(k0)) : Vector{Vector{T}}(undef,length(k0))
    tens_list=Vector{Vector{T}}(undef,length(k0))
    tensN_list=Vector{Vector{T}}(undef,length(k0))
    phi_list=Vector{Matrix{Complex{T}}}(undef,length(k0))
    idx_keep=Vector{Vector{Int}}(undef,length(k0))
    if use_imag_check_EXPERIMENTAL
        idx_keep,residuals=imag_k_check(solver,λs,Uks,Ys,k0,R,all_pts;res_tol=T(res_tol),pad=20,group_size=100,use_chebyshev=use_chebyshev,n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,multithreaded=multithreaded_matrix,verbose=do_per_solve_INFO)
        # now idx_keep gives us those we keep
        @benchit timeit=do_per_solve_INFO "Imag-check selection pass" begin
            @inbounds @showprogress for i in eachindex(k0) # for each batch center k0
                N=boundary_matrix_size(solver,all_pts[i]) # matrix size for this batch 
                if isempty(idx_keep[i])
                    ks_list[i]=return_imag_part ? Complex{T}[] : T[]
                    tens_list[i]=T[]
                    tensN_list[i]=T[]
                    phi_list[i]=Matrix{Complex{T}}(undef,N,0)
                    continue
                end
                λi=λs[i];Uk=Uks[i];Y=Ys[i];idx=idx_keep[i]
                Φ_kept=Matrix{Complex{T}}(undef,N,length(idx))
                for (jj,j) in pairs(idx)
                    @QuantumBilliards.blas_multi_then_1 MAX_BLAS_THREADS mul!(@view(Φ_kept[:,jj]),Uk,@view(Y[:,j]))
                end
                ks_list[i]=return_imag_part ? λi[idx] : real.(λi[idx])
                tens_list[i]=abs.(imag.(λi[idx])) # solver quadrature error proxy. Not a true residual norm
                tensN_list[i]=residuals[i] # for those spurios that were checked and kept, this is the true residual norm ||A(λ)φ||. For unchecked accepted roots, this remains NaN.
                phi_list[i]=Φ_kept
            end
        end
    else
        @benchit timeit=do_per_solve_INFO "Residuals/tensions pass" begin
            @inbounds @showprogress for i in eachindex(k0)
                if isempty(λs[i])
                    ks_list[i]=return_imag_part ? Complex{T}[] : T[]
                    tens_list[i]=T[]
                    tensN_list[i]=T[]
                    phi_list[i]=Matrix{Complex{T}}(undef,boundary_matrix_size(solver,all_pts[i]),0)
                    continue
                end
                idx,Φ_kept,traw,tnorm,_=residual_and_norm_select(solver,λs[i],Uks[i],Ys[i],k0[i],R[i],all_pts[i];res_tol=T(res_tol),matnorm=:one,epss=1e-15,auto_discard_spurious=auto_discard_spurious,collect_logs=false,use_chebyshev=use_chebyshev,n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,multithreaded=multithreaded_matrix)
                ks_list[i]=return_imag_part ? λs[i][idx] : real.(λs[i][idx])
                tens_list[i]=traw
                tensN_list[i]=tnorm
                phi_list[i]=Matrix(Φ_kept)
            end
        end
    end
    nw=length(phi_list)
    n_by_win=Vector{Int}(undef,nw)
    @inbounds for i in 1:nw
        n_by_win[i]=size(phi_list[i],2)
    end
    offs=zeros(Int,nw)
    @inbounds for i in 2:nw
        offs[i]=offs[i-1]+n_by_win[i-1]
    end
    ntot=offs[end]+n_by_win[end]
    ks_all=return_imag_part ? Vector{Complex{T}}(undef,ntot) : Vector{T}(undef,ntot)
    tens_all=Vector{T}(undef,ntot)
    tensN_all=Vector{T}(undef,ntot)
    us_all=Vector{Vector{Complex{T}}}(undef,ntot)
    pts_all=Vector{pts_type}(undef,ntot)
    Threads.@threads for i in 1:nw
        n=n_by_win[i]
        n==0 && continue
        off=offs[i]
        ksi=ks_list[i]
        tr=tens_list[i]
        tn=tensN_list[i]
        Φ=phi_list[i]
        ptsi=all_pts[i]
        @inbounds for j in 1:n
            ks_all[off+j]=ksi[j]
            tens_all[off+j]=tr[j]
            tensN_all[off+j]=tn[j]
            us_all[off+j]=vec(@view Φ[:,j])
            pts_all[off+j]=ptsi
        end
    end
    return ks_all,tens_all,us_all,pts_all,tensN_all
end