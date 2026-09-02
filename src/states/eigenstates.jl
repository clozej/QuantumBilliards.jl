"""
    Eigenstate{K,T,S,Bi,Ba}<:StationaryState

`Eigenstate` is a concrete type representing a numerically computed eigenstate
of a quantum billiard at a given wavenumber.

## Description
Eigenstates are produced by [`compute_eigenstate`](@ref), which combines a
sweep or accelerated solver, a basis and a billiard geometry to obtain the
expansion coefficients `vec` in `basis` and an estimate of the tension `ten`,
quantifying how well the boundary condition is satisfied. Coefficients
smaller in magnitude than the numerical precision `eps`
are set to zero when the state is constructed.

## Attributes
* `k`: The wavenumber of the eigenstate, as refined by the solver.
* `k_basis`: The wavenumber at which `basis` was evaluated to obtain `vec` (may differ slightly from `k` for accelerated solvers).
* `vec`: Expansion coefficients of the eigenstate in `basis`.
* `ten`: Tension of the solution, measuring the residual boundary condition violation.
* `dim`: Dimension of `vec` (and of `basis`).
* `eps`: Numerical precision threshold below which coefficients of `vec` are treated as zero.
* `solver`: The solver (`S<:AbsSolver`) used to compute the eigenstate.
* `basis`: The basis (`Ba<:AbsBasis`), resized/evaluated at `k_basis`, in which `vec` is expressed.
* `billiard`: The billiard (`Bi<:AbsBilliard`) the eigenstate is defined on.

## API
The following functions can be evaluated for this type:
- [`compute_eigenstate`](@ref)
- [`boundary_function`](@ref)
- [`momentum_function`](@ref)
- [`wavefunction`](@ref)
- [`husimi_function`](@ref)
"""
struct Eigenstate{K,T,S,Bi,Ba}<:StationaryState
    k::K
    k_basis::K
    vec::Vector{K}
    ten::T
    dim::Int64
    eps::T
    solver::S
    basis::Ba
    billiard::Bi
end

"""
    Eigenstate(k::K,vec::Vector{K},ten::T,solver::S,basis::Ba,billiard::Bi) where {K<:Number,T<:Real,S<:AbsSolver,Ba<:AbsBasis,Bi<:AbsBilliard} → state::Eigenstate

Construct an [`Eigenstate`](@ref) with `k_basis=k`.
"""
Eigenstate(k::K,vec::Vector{K},ten::T,solver::S,basis::Ba,billiard::Bi) where {K<:Number,T<:Real,S<:AbsSolver,Ba<:AbsBasis,Bi<:AbsBilliard}=Eigenstate(k,k,vec,ten,solver,basis,billiard)

"""
    Eigenstate(k::K,k_basis::K,vec::Vector{K},ten::T,solver::S,basis::Ba,billiard::Bi) where {K<:Number,T<:Real,S<:AbsSolver,Ba<:AbsBasis,Bi<:AbsBilliard} → state::Eigenstate

Construct an [`Eigenstate`](@ref) with independently specified eigenstate and
basis wavenumbers. Real coefficients below the numerical threshold are set to
zero.

## Arguments
* `k::K`: Refined eigenstate wavenumber.
* `k_basis::K`: Wavenumber at which the basis was evaluated.
* `vec::Vector{K}`: Expansion coefficients.
* `ten::T`: Solver-specific tension or residual.
* `solver::S`: Solver used to compute the state.
* `basis::Ba`: Basis containing the expansion.
* `billiard::Bi`: Billiard geometry.

## Returns
* `state::Eigenstate`: Constructed eigenstate.
"""
function Eigenstate(k::K,k_basis::K,vec::Vector{K},ten::T,solver::S,basis::Ba,billiard::Bi) where {K<:Number,T<:Real,S<:AbsSolver,Ba<:AbsBasis,Bi<:AbsBilliard}
    eps=T(K===Float32 ? 1e-8 : 1e-16)
    filtered_vec=K<:Real ? K[abs(v)>eps ? v : zero(K) for v in vec] : vec
    return Eigenstate(k,k_basis,filtered_vec,ten,length(vec),eps,solver,basis,billiard)
end

"""
    BasisState{K,T,Ba} <: StationaryState

Stationary state representing one basis function.

## Attributes
* `k::K`: Wavenumber.
* `k_basis::K`: Basis-evaluation wavenumber, equal to `k`.
* `vec::Vector{T}`: Unit coefficient vector.
* `idx::Int64`: Selected basis-function index.
* `dim::Int64`: Basis dimension.
* `eps::T`: Numerical precision threshold.
* `basis::Ba`: Basis containing the represented function.
"""
struct BasisState{K,T,Ba}<:StationaryState
    k::K
    k_basis::K
    vec::Vector{T}
    idx::Int64
    dim::Int64
    eps::T
    basis::Ba
end

"""
    BasisState(basis::Ba,k::T,i::Int) where {T<:Real,Ba<:AbsBasis} → state::BasisState

Construct a [`BasisState`](@ref) representing the `i`-th basis function at wavenumber `k`.

## Arguments
* `basis::Ba`: Basis containing the represented function.
* `k::T`: Wavenumber.
* `i::Int`: Basis-function index.

## Returns
* `state::BasisState`: State with a unit coefficient vector at index `i`.
"""
function BasisState(basis::Ba,k::T,i::Int) where {T<:Real,Ba<:AbsBasis}
    dim=basis.dim
    eps=T===Float32 ? T(1e-8) : T(1e-16)
    vec=zeros(T,dim)
    vec[i]=one(T)
    return BasisState(k,k,vec,i,dim,eps,basis)
end