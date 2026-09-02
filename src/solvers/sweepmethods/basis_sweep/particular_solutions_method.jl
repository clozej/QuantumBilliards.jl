"""
    struct ParticularSolutionsMethod{T} <: SweepSolver where {T<:Real}

A data structure to hold configuration parameters for the "Particular Solutions Method."
This method typically uses a scaled number of boundary and interior points based on
the wavenumber `k`. It also defines sampling strategies and minimal dimension settings.

# Fields
- `dim_scaling_factor::T`: A real scale factor for matrix dimension or basis dimension.
- `pts_scaling_factor::Vector{T}`: A list of scaling factors for boundary points across different curves.
- `int_pts_scaling_factor::T`: Scaling factor for interior points.
- `sampler::Vector`: A list of samplers (e.g. Gauss-Legendre) for boundary integration.
- `eps::T`: A small numeric epsilon, set to `eps(T)`.
- `min_dim::Int64`: The minimal dimension of a matrix or basis.
- `min_pts::Int64`: The minimal number of boundary points.
- `min_int_pts::Int64`: The minimal number of interior points.
"""
struct ParticularSolutionsMethod{T} <: SweepSolver where {T<:Real}
    dim_scaling_factor::T
    pts_scaling_factor::Vector{T}
    int_pts_scaling_factor::T
    sampler::Vector
    eps::T
    min_dim::Int64
    min_pts::Int64
    min_int_pts::Int64
end

"""
    ParticularSolutionsMethod(dim_scaling_factor::T,pts_scaling_factor::Union{T,Vector{T}},int_pts_scaling_factor::T;min_dim::Int=100,min_pts::Int=500,min_int_pts::Int=500) where {T<:Real} → ParticularSolutionsMethod{T}

Construct a `ParticularSolutionsMethod{T}` object with default Gauss-Legendre samplers.

# Arguments
- `dim_scaling_factor::T`: A scaling factor that controls basis dimension.
- `pts_scaling_factor::Union{T,Vector{T}}`: A scaling factor that controls final matrix dimension for boundary point discretization.
- `int_pts_scaling_factor::T`: Scaling factor for interior points. Usually use the same as the boundary point scaling factor.
- `min_dim::Int`: The minimal dimension used in computations (default 100).
- `min_pts::Int`: The minimal number of boundary points (default 500).
- `min_int_pts::Int`: The minimal number of interior points (default 500).

# Returns
- `ParticularSolutionsMethod{T}`: A new solver configuration with a Gauss-Legendre sampler
  in the `sampler` field, and a numerically safe epsilon in the `eps` field.
"""
function ParticularSolutionsMethod(dim_scaling_factor::T,pts_scaling_factor::Union{T,Vector{T}},int_pts_scaling_factor::T;min_dim::Int=100,min_pts::Int=500,min_int_pts::Int=500) where {T<:Real}
    d=dim_scaling_factor
    bs=pts_scaling_factor isa T ? [pts_scaling_factor] : pts_scaling_factor
    sampler=[BilliardGeometry.GaussLegendreNodes()]
    epsilon=max(eps(T),T(1e-14)) # set a floor on eps to avoid numerical issues in the generalized eigenvalue solver
    return ParticularSolutionsMethod(d,bs,int_pts_scaling_factor,sampler,epsilon,min_dim,min_pts,min_int_pts)
end

"""
    ParticularSolutionsMethod(dim_scaling_factor::T,pts_scaling_factor::Union{T,Vector{T}},int_pts_scaling_factor::T,samplers::Vector{<:AbsSampler};min_dim::Int=100,min_pts::Int=500,min_int_pts::Int=500) where {T<:Real} → ParticularSolutionsMethod{T}

Construct a `ParticularSolutionsMethod{T}` object with user-provided samplers.

# Arguments
#- `dim_scaling_factor::T`: A scaling factor that controls basis dimension.
#- `pts_scaling_factor::Union{T,Vector{T}}`: A scaling factor that controls final matrix dimension for boundary point discretization.
#- `int_pts_scaling_factor::T`: Scaling factor for interior points. Usually use the same as the boundary point scaling factor.
#- `samplers::Vector{<:BilliardGeometry.AbsSampler}`: A vector of samplers (e.g. `BilliardGeometry.GaussLegendreNodes()`, `BilliardGeometry.PolarSampler()`, etc.).
- `min_dim::Int`: Minimal dimension (default 100).
- `min_pts::Int`: Minimal boundary points (default 500).
- `min_int_pts::Int`: Minimal interior points (default 500).

# Returns
- `ParticularSolutionsMethod{T}`: A solver configuration with the given `samplers` and other settings.
"""
function ParticularSolutionsMethod(dim_scaling_factor::T,pts_scaling_factor::Union{T,Vector{T}},int_pts_scaling_factor::T,samplers::Vector{<:BilliardGeometry.AbsSampler};min_dim::Int=100,min_pts::Int=500,min_int_pts::Int=500) where {T<:Real}
    d=dim_scaling_factor
    bs=pts_scaling_factor isa T ? [pts_scaling_factor] : pts_scaling_factor
    epsilon=max(eps(T),T(1e-14)) # set a floor on eps to avoid numerical issues in the generalized eigenvalue solver
    return ParticularSolutionsMethod(d,bs,int_pts_scaling_factor,samplers,epsilon,min_dim,min_pts,min_int_pts)
end

"""
    evaluate_points(solver::ParticularSolutionsMethod,billiard::Bi,k) where {Bi<:BilliardGeometry.AbsBilliard} → BoundaryPoints

Generate boundary and interior points for the Particular Solutions Method. Uses the solver's scaling
factors and minimal point constraints to decide how many points to sample on each boundary curve
and how many interior points to sample.

# Arguments
- `solver::ParticularSolutionsMethod`: The PSM solver configuration.
- `billiard::Bi<:BilliardGeometry.AbsBilliard`: A geometry or domain object holding boundary curves.
- `k::Real`: Wavenumber, used to scale the number of sampled points.

# Returns
- `pts::BoundaryPoints`: Boundary discretization containing `xy`, `normal`, `s`, `ds` and `xy_int`.
"""
function evaluate_points(solver::ParticularSolutionsMethod,billiard::Bi,k) where {Bi<:BilliardGeometry.AbsBilliard}
    bs,samplers=adjust_scaling_and_samplers(solver,billiard)
    curves=BilliardGeometry.get_boundary_curves(billiard)
    T=eltype(solver.pts_scaling_factor)
    Ns=_determine_bp_sizes(curves,bs,k)
    M=length(Ns)
    xy_all=Vector{Vector{SVector{2,T}}}(undef,M)
    normal_all=Vector{Vector{SVector{2,T}}}(undef,M)
    s_all=Vector{Vector{T}}(undef,M)
    ds_all=Vector{Vector{T}}(undef,M)
    L0=zero(T)
    @inbounds for i in eachindex(curves)
        xy,normal,s,ds=boundary_coords(curves[i],samplers[i],Ns[i])
        xy_all[i]=xy
        normal_all[i]=normal
        s_all[i]=s.+L0
        ds_all[i]=ds
        L0+=curves[i].length
    end
    M_int=max(solver.min_int_pts,round(Int,k*L0*solver.int_pts_scaling_factor/(2*pi)))
    xy_int=random_interior_points(billiard,M_int)
    return BoundaryPoints(vcat(xy_all...);normal=vcat(normal_all...),s=vcat(s_all...),ds=vcat(ds_all...),xy_int=xy_int)
end

"""
    construct_matrices_benchmark(solver::ParticularSolutionsMethod,basis::Ba,pts::BoundaryPoints,k;multithreaded::Bool=true) where {Ba<:AbsBasis}

Construct the basis matrices for boundary points and interior points, with timing information.

# Arguments
- `solver::ParticularSolutionsMethod`: Solver config, specifying scaling factors.
- `basis::Ba<:AbsBasis`: A basis to evaluate (e.g. FourierBessel basis).
- `pts::BoundaryPoints{T}`: The boundary/interior points to evaluate.
- `k::Real`: Wavenumber for the basis.
- `multithreaded::Bool=true`: If the matrix construction should be multithreaded.

# Returns
- `(B, B_int)`: A pair of matrices:
  - `B::Matrix`: The basis matrix evaluated at boundary points.
  - `B_int::Matrix`: The basis matrix evaluated at interior points.
"""
function construct_matrices_benchmark(solver::ParticularSolutionsMethod,basis::Ba,pts::BoundaryPoints,k;multithreaded::Bool=true) where {Ba<:AbsBasis}
    pts_bd=pts.xy
    pts_int=pts.xy_int
    @time "boundary" B=basis_matrix(basis,k,pts_bd;multithreaded=multithreaded)
    @time "interior" B_int=basis_matrix(basis,k,pts_int;multithreaded=multithreaded)
    return B,B_int
end

"""
    construct_matrices(solver::ParticularSolutionsMethod,basis::Ba,pts::BoundaryPoints,k;multithreaded::Bool=true) where {Ba<:AbsBasis}

Construct two matrices for the Particular Solutions Method: one for boundary points, one for interior points.
These represent the basis functions evaluated at the domain's boundary and interior.

# Arguments
- `solver::ParticularSolutionsMethod`: The PSM solver config.
- `basis::Ba<:AbsBasis`: A basis type implementing `basis_matrix(...)`.
- `pts::BoundaryPoints{T}`: Contains boundary (`xy`) and interior (`xy_int`) points.
- `k::Real`: Wavenumber.
- `multithreaded::Bool=true`: If the matrix construction should be multithreaded.

# Returns
- `(B, B_int)`:
  - `B::Matrix`: Basis matrix at boundary points.
  - `B_int::Matrix`: Basis matrix at interior points.
"""
function construct_matrices(solver::ParticularSolutionsMethod,basis::Ba,pts::BoundaryPoints,k;multithreaded::Bool=true) where {Ba<:AbsBasis}
    pts_bd=pts.xy
    pts_int=pts.xy_int
    @blas_1 begin
        B=basis_matrix(basis,k,pts_bd;multithreaded=multithreaded)
        B_int=basis_matrix(basis,k,pts_int;multithreaded=multithreaded)
    end
    return B,B_int
end

function solve_full(solver::ParticularSolutionsMethod,basis::Ba,pts::BoundaryPoints,k;multithreaded::Bool=true) where {Ba<:AbsBasis}
    B,B_int=construct_matrices(solver,basis,pts,k;multithreaded=multithreaded)
    @blas_multi_then_1 MAX_BLAS_THREADS solution=svdvals(B,B_int)
    return minimum(solution)
end

#### INTERNAL - FIND NUMERICAL RANK FROM FACTORIZED MATRIX ####
@inline function _numerical_rank_from_F(F,tol::Real)
    A=F.factors
    n=min(size(A,1),size(A,2))
    n==0&&return 0
    t=tol*abs(@inbounds A[1,1])
    @inbounds for i=n:-1:1
        if abs(A[i,i])>t
            return i
        end
    end
    return 0
end

function solve_with_rank_reduction(solver::ParticularSolutionsMethod,basis::Ba,pts::BoundaryPoints,k;multithreaded::Bool=true,tol=1e-10) where {Ba<:AbsBasis}
    # tol is adjustable, based on how the R[1,1] scales below with k. Tested up to k=500 it seems fine with 1e-14
    B,C=construct_matrices(solver,basis,pts,k;multithreaded=multithreaded)
    T=eltype(B)
    @blas_multi_then_1 MAX_BLAS_THREADS begin
        F=qr!(C,ColumnNorm()) # rank-revealing QR with column pivoting: C*P = Q*R. Overwrite C since we do not need it anymore
        r=QuantumBilliards._numerical_rank_from_F(F,tol) # numerical rank r from packed factors (no copy)
        r==0&&return Inf # in case degenerate fail
        Rview=@views UpperTriangular(view(F.factors,1:r,1:r)) # triangular view onto R (no copy)
        piv=F.p # permutation vector piv such that C[:,piv] = Q*R
        B=@views B[:,piv[1:r]]/Rview # Br = B[:,piv[1:r]] * R^{-1} via triangular solve (stable, no inv) and overwrite B
        r=size(B,2) # number of columns after reduction
        B_sq=Matrix{T}(undef,r,r) # small matrix B'B preallocate
        BLAS.syrk!('U','T',one(T),B,zero(T),B_sq) # real matrix so transpose is ok
        _symmetrize_from_upper!(B_sq) # fill the lower triangle inplace
        return sqrt(abs(eigmin(Symmetric(B_sq)))) # smallest singular value via eigmin(B'B). This is usually faster than using Krylov'S smallest svd since the matrix rank reduction is significant
    end
end

"""
    solve(solver::ParticularSolutionsMethod,basis::Ba,pts::BoundaryPoints,k;multithreaded::Bool=true,use_rank_reduction::Bool=true,tol=1e-10) where {Ba<:AbsBasis} → Real

Solve the Particular Solutions Method by constructing `(B, B_int)` and computing a measure
(e.g. minimum singular value) that indicates how well the interior and boundary constraints
are satisfied via the `GeneralizedSVD` object. For mode reading it is based on the paper:
https://people.maths.ox.ac.uk/trefethen/publication/PDF/2005_112.pdf.

The idea is to represent the minimization of the boundary tension of the wavefunctions as a generalized eigenvalue problem `F * x = λ * G * x` as per the afformentioned paper. There are 2 approaches implemented here: the full SVD approach and a rank-reduced approach that offers more performance compared with the full SVD approach. If very high basis dimensions are used, this is recommended, but sometimes tolerance needs adjustement.

# Arguments
- `solver::ParticularSolutionsMethod`: PSM solver config.
- `basis::Ba<:AbsBasis`: The basis (e.g. a trigonometric or radial basis).
- `pts::BoundaryPoints{T}`: Boundary points for evaluation.
- `k::Real`: Wavenumber or frequency-like parameter.
- `multithreaded::Bool=true`: If the matrix construction should be multithreaded.
- `use_rank_reduction::Bool=true`: If true, uses a rank-reduced solver variant that is more stable for large bases & offers unparalleled performance compared with the full SVD approach. If very high basis dimensions are used, this is recommended, but sometimes tolerance needs adjustement.
- `tol::Real=1e-10`: Tolerance for numerical rank determination in the rank-reduced solver. Adjust if necessary (you see too much spiking in the tension plot). A good choice is to start with 1e-10 and lower it if necessary. Rule of thumb here is that tol should increase with k since R[1,1] decreases with k.

# Returns
- `::Real`: The minimum generalized singular value (or similar measure). Lower values can
  indicate a better "fit" to the PDE boundary conditions.
"""
function solve(solver::ParticularSolutionsMethod,basis::Ba,pts::BoundaryPoints,k;multithreaded::Bool=true,use_rank_reduction::Bool=true,tol=1e-10) where {Ba<:AbsBasis}
    if use_rank_reduction
        return solve_with_rank_reduction(solver,basis,pts,k;multithreaded=multithreaded,tol=tol)
    else
        return solve_full(solver,basis,pts,k;multithreaded=multithreaded)
    end
end

# INTERNAL
function solve_INFO(solver::ParticularSolutionsMethod,basis::Ba,pts::BoundaryPoints,k;multithreaded::Bool=true,tol=1e-10) where {Ba<:AbsBasis}
    @time "Matrix construction" B,C=construct_matrices(solver,basis,pts,k;multithreaded=multithreaded)
    T=eltype(B)
    @blas_multi_then_1 MAX_BLAS_THREADS begin
        @time "QR" F=qr!(C,ColumnNorm()) # rank-revealing QR with column pivoting: C*P = Q*R. This is the main trick. Overwrite C since we do not need it anymore
        @time "Finding r" r=_numerical_rank_from_F(F,tol) # numerical rank r: keep diagonal entries of R down to a relative threshold and discard near-null interior directions that cause spurious minima
        r==0&&return Inf # in case degenerate fail
        @time "Triangular view" Rview=@views UpperTriangular(view(F.factors,1:r,1:r)) # well-determined r×r block of R as a view (no copy)
        piv=F.p # permutation vector piv such that C[:,piv] = Q*R
        println("rank: ",r)
        println("size decrease: ",size(B,2)-r)
        @time "/ solve" B=@views B[:,piv[1:r]]/Rview # Br = B[:,piv[1:r]] * R^{-1} via triangular solve (stable, no inv) and overwrite B
        @time "eigmin solve" v=sqrt(real(eigmin(B'*B))) # smallest singular value via eigmin(B'B)
        return v
    end
end

"""
    solve_vect(solver::ParticularSolutionsMethod,basis::Ba,pts::BoundaryPoints,k;multithreaded::Bool=true,tol=1e-10) where {Ba<:AbsBasis}

Returns the smallest singular value and the basis expansion coefficient vector for use in `boundary_function`.

# Arguments
- `solver::ParticularSolutionsMethod`: PSM solver configuration.
- `basis::Ba<:AbsBasis`: Basis used to assemble boundary/interior matrices.
- `pts::BoundaryPoints{T}`: Boundary points for evaluation.
- `k::Real`: Wavenumber.
- `multithreaded::Bool=true`: Pass-through for matrix construction.
- `tol::Real=1e-10`: Tolerance for numerical rank determination in the rank-reduced solver. Adjust if necessary (you see too much spiking in the tension plot). A good choice is to start with 1e-10 and lower it if necessary. Rule of thumb here is that tol should increase with k since R[1,1] decreases with k.

# Returns
- `μ::Real`: the stabilized analogue of the smallest generalized singular value.
- `chat::Vector{Real}`: Coefficient vector in the given `basis`.
"""
function solve_vect(solver::ParticularSolutionsMethod,basis::Ba,pts::BoundaryPoints,k;multithreaded::Bool=true,tol=1e-10) where {Ba<:AbsBasis}
    B,C=construct_matrices(solver,basis,pts,k;multithreaded=multithreaded)
    @blas_multi_then_1 MAX_BLAS_THREADS begin
        T=eltype(B)
        F=qr(C,ColumnNorm()) # rank-revealing QR with column pivoting: C*P = Q*R. This is the main trick
        R=UpperTriangular(F.R) # for fast triangular solves, just in case API changes
        piv=F.p # permutation vector piv such that C[:,piv] = Q*R
        r=findlast(i->abs(R[i,i])>tol*abs(R[1,1]),1:min(size(R)...)) # numerical rank r: keep diagonal entries of R down to a relative threshold and discard near-null interior directions that cause spurious minima. #TODO could use _numerical_rank_from_F here which was implemented later, but this function is not in the hot loop
        isnothing(r)&&return (Inf,zeros(T,size(B,2))) # in case degenerate fail
        Rr=R[1:r,1:r] # well-determined r×r block on Q
        Br=B[:,piv[1:r]]/Rr # Br = B[:,piv[1:r]] * Rr^{-1} via triangular solve (stable, no inv)
        _,S,Vt=LAPACK.gesvd!('A','A',Br) # SVD(Br) = U*Diag(S)*transpose(V); the smallest singular value gives the minimum of ‖Br y‖ subject to ‖y‖=1, and its right singular vector is the minimizer y. In principle could use Krylov with :SM but this is not called in the bottleneck eigenvalue search, and it actually also fails! The only way to use Krylov is to form Br'*Br and find its smallest eigenvalue (but for large k the smallest singular value could be below machine precision since we get the eigenvalue of Br'*Br and then take sqrt).
        idx=findmin(S)[2]
        mu=S[idx]
        u_mu=Vt[idx,:]
        y=real.(u_mu)
        chat=zeros(T,size(B,2))
        chat[piv[1:r]]=Rr\y # back-substitute: c[piv[1:r]]=Rr^{-1} y; rest are zeros
        return mu,chat
    end
end