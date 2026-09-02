"""
    solve_wavenumber(solver::VerginiSaracenoSolver{T},basis::Ba,billiard::Bi,k,dk;multithreaded::Bool=true) where {T<:Real,Ba<:AbsBasis,Bi<:BilliardGeometry.AbsBilliard} → (k0,t0)

Find the Vergini–Saraceno eigenvalue candidate closest to the target wavenumber
`k` and return it together with its tension.

## Description
The basis dimension is scaled from the boundary length and `k`, bounded below by
`solver.min_dim`. The basis is resized with [`resize_basis`](@ref), boundary
points are generated with [`evaluate_points`](@ref), and one
Vergini–Saraceno generalized eigenproblem is solved with [`solve`](@ref).
The candidate closest to `k` is returned.

## Arguments
* `solver::VerginiSaracenoSolver{T}`: Vergini–Saraceno solver.
* `basis::Ba`: Basis used to approximate the eigenstates.
* `billiard::Bi`: Billiard geometry.
* `k`: Target wavenumber.
* `dk`: Half-width of the local wavenumber window.

## Keyword Arguments
* `multithreaded::Bool=true`: Enable multithreaded matrix construction.

## Returns
* `k0`: Candidate wavenumber closest to `k`.
* `t0`: Vergini–Saraceno tension associated with `k0`.
"""
function solve_wavenumber(solver::VerginiSaracenoSolver{T},basis::Ba,billiard::Bi,k,dk;multithreaded::Bool=true) where {T<:Real,Ba<:AbsBasis,Bi<:BilliardGeometry.AbsBilliard}
    L=BilliardGeometry.CompositeCurve(BilliardGeometry.get_boundary_curves(billiard)).length
    dim=max(solver.min_dim,round(Int,L*k*solver.dim_scaling_factor/(2*pi)))
    new_basis=resize_basis(basis,billiard,dim,k)
    pts=evaluate_points(solver,billiard,k)
    ks,ts=solve(solver,new_basis,pts,k,dk;multithreaded)
    idx=findmin(abs.(ks.-k))[2]
    return ks[idx],ts[idx]
end
"""
    solve_wavenumber_beyn(solver::BeynSolver{T},billiard::Bi,k::T,dk::T;nq::Int=48,r::Int=48,svd_tol::Real=1e-12,res_tol::Real=1e-9,auto_discard_spurious::Bool=true,multithreaded_matrix::Bool=true,use_chebyshev::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,return_imag_part::Bool=false) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard} → (k0,t0)

Perform one validated Beyn contour solve centered at the target wavenumber `k`.

## Description
A circular contour of radius `dk` centered at `complex(k)` is solved with
[`solve_vect`](@ref). Provisional roots are filtered with
[`residual_and_norm_select`](@ref), and the retained root whose real part is
closest to `k` is returned.

## Arguments
* `solver::BeynSolver{T}`: Beyn-compatible boundary-integral solver.
* `billiard::Bi`: Billiard geometry.
* `k::T`: Real center of the Beyn contour.
* `dk::T`: Radius of the circular contour.

## Keyword Arguments
* `nq::Int=48`: Number of contour quadrature nodes.
* `r::Int=48`: Initial random probe rank.
* `svd_tol::Real=1e-12`: Singular-value threshold for Beyn rank detection.
* `res_tol::Real=1e-9`: Raw nonlinear residual threshold.
* `auto_discard_spurious::Bool=true`: Reject candidates whose residual is at least `res_tol`.
* `multithreaded_matrix::Bool=true`: Enable multithreaded boundary-matrix construction.
* `use_chebyshev::Bool=true`: Enable Chebyshev-accelerated matrix construction.
* `n_panels_h::Int=15000`: Hankel Chebyshev panel count.
* `M_h::Int=5`: Hankel Chebyshev polynomial degree.
* `n_panels_j::Int=10000`: Bessel-J Chebyshev panel count.
* `M_j::Int=5`: Bessel-J Chebyshev polynomial degree.
* `return_imag_part::Bool=false`: Return the complex approximate root instead of only its real part.

## Returns
* `k0`: Retained Beyn root closest to `k`, or `NaN` if no validated root is found.
* `t0::T`: Raw nonlinear residual norm, or `Inf` if no validated root is found.
"""
function solve_wavenumber_beyn(solver::BeynSolver{T},billiard::Bi,k::T,dk::T;nq::Int=48,r::Int=48,svd_tol::Real=1e-12,res_tol::Real=1e-9,auto_discard_spurious::Bool=true,multithreaded_matrix::Bool=true,use_chebyshev::Bool=true,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,return_imag_part::Bool=false) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    basis=AbstractHankelBasis()
    pts=evaluate_points(solver,billiard,k)
    k0=complex(k)
    λ,Uk,Y=solve_vect(solver,basis,pts,k0,dk;nq=nq,r=r,svd_tol=svd_tol,rng=MersenneTwister(0),multithreaded=multithreaded_matrix,use_chebyshev=use_chebyshev,n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j)
    isempty(λ)&&return (return_imag_part ? Complex{T}(NaN,NaN) : T(NaN)),T(Inf)
    idx,_,tens,_,_=residual_and_norm_select(solver,λ,Uk,Y,k0,dk,pts;res_tol=T(res_tol),auto_discard_spurious=auto_discard_spurious,use_chebyshev=use_chebyshev,n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j,multithreaded=multithreaded_matrix)
    isempty(idx)&&return (return_imag_part ? Complex{T}(NaN,NaN) : T(NaN)),T(Inf)
    λkeep=λ[idx]
    j=argmin(abs.(real.(λkeep).-k))
    return return_imag_part ? λkeep[j] : real(λkeep[j]),tens[j]
end
"""
    solve_wavenumber_ebim(solver::EBIMSolver,billiard::Bi,k::T,dk::T;use_lapack_raw::Bool=false,multithreaded_matrices::Bool=false,use_krylov::Bool=true,nev::Int=5,use_chebyshev::Bool=false,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,return_imag_part::Bool=false) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard} → (k0,t0)

Perform one local EBIM expansion around the target wavenumber `k`.

## Description
The boundary geometry is evaluated at `k`, the matrices `A(k)`, `A'(k)` and
`A''(k)` are constructed, and one dense or Krylov EBIM generalized eigenproblem
is solved. The accepted corrected root whose real part is closest to `k` is
returned.

## Arguments
* `solver::EBIMSolver`: EBIM-compatible boundary-integral solver.
* `billiard::Bi`: Billiard geometry.
* `k::T`: Expansion wavenumber.
* `dk::T`: Local correction half-width.

## Keyword Arguments
* `use_lapack_raw::Bool=false`: Use the low-level dense generalized eigensolver.
* `multithreaded_matrices::Bool=false`: Enable multithreaded boundary-matrix construction.
* `use_krylov::Bool=true`: Use the shift-invert Krylov eigensolver instead of the dense generalized EVP.
* `nev::Int=5`: Number of Krylov eigenpairs requested.
* `use_chebyshev::Bool=false`: Use derivative-aware Chebyshev matrix construction.
* `n_panels_h::Int=15000`: Hankel Chebyshev panel count.
* `M_h::Int=5`: Hankel Chebyshev polynomial degree.
* `n_panels_j::Int=10000`: Bessel-J Chebyshev panel count.
* `M_j::Int=5`: Bessel-J Chebyshev polynomial degree.
* `return_imag_part::Bool=false`: Return the complex corrected root instead of only its real part.

## Returns
* `k0`: Corrected EBIM root closest to `k`, or `NaN` if no root is found.
* `t0::T`: Absolute EBIM correction magnitude `|ε₁+ε₂|`, or `Inf` if no root is found.
"""
function solve_wavenumber_ebim(solver::EBIMSolver,billiard::Bi,k::T,dk::T;use_lapack_raw::Bool=false,multithreaded_matrices::Bool=false,use_krylov::Bool=true,nev::Int=5,use_chebyshev::Bool=false,n_panels_h::Int=15000,M_h::Int=5,n_panels_j::Int=10000,M_j::Int=5,return_imag_part::Bool=false) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    pts=evaluate_points(solver,billiard,k)
    A,dA,ddA=allocate_ebim_matrices(solver,pts)
    if use_chebyshev
        cache=build_ebim_cheb_cache(solver,pts,T[k];n_panels_h=n_panels_h,M_h=M_h,n_panels_j=n_panels_j,M_j=M_j)
        λs,tens=solve!(solver,A,dA,ddA,pts,k,dk,cache,1;use_lapack_raw=use_lapack_raw,multithreaded=multithreaded_matrices,use_krylov=use_krylov,nev=nev,return_imag_part=return_imag_part)
    else
        λs,tens=solve!(solver,A,dA,ddA,pts,k,dk;use_lapack_raw=use_lapack_raw,multithreaded=multithreaded_matrices,use_krylov=use_krylov,nev=nev,return_imag_part=return_imag_part)
    end
    isempty(λs)&&return (return_imag_part ? Complex{T}(NaN,NaN) : T(NaN)),T(Inf)
    j=argmin(abs.(real.(λs).-k))
    return λs[j],tens[j]
end
"""
    solve_wavenumber(method::Symbol,solver,billiard::Bi,k::T,dk::T;kwargs...) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard} → (k0,t0)

Dispatch one local boundary-integral wavenumber solve to the Beyn or EBIM
backend.

## Arguments
* `method::Symbol`: Spectral method, either `:beyn` or `:ebim`.
* `solver`: Boundary-integral solver compatible with the selected method.
* `billiard::Bi`: Billiard geometry.
* `k::T`: Target wavenumber.
* `dk::T`: Beyn contour radius or EBIM local correction half-width.

## Keyword Arguments
* `kwargs...`: Keyword arguments forwarded to the selected backend.

## Returns
* `k0`: Candidate closest to `k`.
* `t0`: Method-specific quality measure: the raw nonlinear residual for Beyn or the absolute local correction for EBIM.
"""
function solve_wavenumber(method::Symbol,solver,billiard::Bi,k::T,dk::T;kwargs...) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    method===:beyn&&return solve_wavenumber_beyn(solver,billiard,k,dk;kwargs...)
    method===:ebim&&return solve_wavenumber_ebim(solver,billiard,k,dk;kwargs...)
    throw(ArgumentError("Unknown wavenumber method $method; expected :beyn or :ebim"))
end
"""
    solve_spectrum(solver::VerginiSaracenoSolver{T},basis::Ba,billiard::Bi,k,dk;multithreaded::Bool=true) where {T<:Real,Ba<:AbsBasis,Bi<:BilliardGeometry.AbsBilliard} → (ks,ts)

Compute all Vergini–Saraceno eigenvalue candidates and tensions found in one
diagonalization near `k`.

## Description
The basis dimension is scaled from the boundary length and `k`, bounded below by
`solver.min_dim`. The basis is resized with [`resize_basis`](@ref), boundary
points are generated with [`evaluate_points`](@ref), and one
Vergini–Saraceno generalized eigenproblem is solved with [`solve`](@ref).

## Arguments
* `solver::VerginiSaracenoSolver{T}`: Vergini–Saraceno solver.
* `basis::Ba`: Basis used to approximate the eigenstates.
* `billiard::Bi`: Billiard geometry.
* `k`: Wavenumber around which the solve is performed.
* `dk`: Half-width of the local wavenumber window.

## Keyword Arguments
* `multithreaded::Bool=true`: Enable multithreaded matrix construction.

## Returns
* `ks`: Candidate wavenumbers found in the window.
* `ts`: Vergini–Saraceno tensions associated with `ks`.
"""
function solve_spectrum(solver::VerginiSaracenoSolver{T},basis::Ba,billiard::Bi,k,dk;multithreaded::Bool=true) where {T<:Real,Ba<:AbsBasis,Bi<:BilliardGeometry.AbsBilliard}
    L=BilliardGeometry.CompositeCurve(BilliardGeometry.get_boundary_curves(billiard)).length
    dim=max(solver.min_dim,round(Int,L*k*solver.dim_scaling_factor/(2*pi)))
    new_basis=resize_basis(basis,billiard,dim,k)
    pts=evaluate_points(solver,billiard,k)
    ks,ts=solve(solver,new_basis,pts,k,dk;multithreaded)
    return ks,ts
end
"""
    solve_spectrum(method::Symbol,solver,billiard::Bi,k1::T,k2::T;kwargs...) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}

Dispatch a boundary-integral spectral sweep over `[k1,k2]` to the Beyn or EBIM
backend.

## Arguments
* `method::Symbol`: Spectral method, either `:beyn` or `:ebim`.
* `solver`: Boundary-integral solver compatible with the selected method.
* `billiard::Bi`: Billiard geometry.
* `k1::T`: Lower wavenumber bound.
* `k2::T`: Upper wavenumber bound.

## Keyword Arguments
* `kwargs...`: Keyword arguments forwarded to the selected spectral backend.

## Returns
The return value of [`solve_spectrum_beyn`](@ref) for `:beyn` or
[`solve_spectrum_ebim`](@ref) for `:ebim`.
"""
function solve_spectrum(method::Symbol,solver,billiard::Bi,k1::T,k2::T;kwargs...) where {T<:Real,Bi<:BilliardGeometry.AbsBilliard}
    method===:beyn&&return solve_spectrum_beyn(solver,billiard,k1,k2;kwargs...)
    method===:ebim&&return solve_spectrum_ebim(solver,billiard,k1,k2;kwargs...)
    throw(ArgumentError("Unknown spectrum method $method; expected :beyn or :ebim"))
end