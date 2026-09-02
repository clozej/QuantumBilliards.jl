# Construct an antisymmetric vector by reversing `x[2:end]`, negating the
# reversed values and prepending them to `x`.
function antisym_vec(x::AbstractVector)
    v=reverse(-x[2:end])
    return append!(v,x)
end

# Extract the boundary nodes satisfying `s<=width` and construct the corresponding
# symmetric relative-coordinate window together with its reflected quadrature weights.
function _husimi_symmetric_window(s::AbstractVector{T},ds::AbstractVector{T},width::Real) where {T<:Real}
    mask=s.<=width
    x=collect(s[mask])
    dx=collect(ds[mask])
    idx=length(x)
    x=antisym_vec(x)
    dx=vcat(reverse(dx[2:end]),dx)
    return x,dx,idx
end

# check of teh spacing is equidistant so one can use a sliding window approach for the Husimi function evaluation
function _husimi_uniform_arclength_grid(s::AbstractVector{T},ds::AbstractVector{T};rtol::Real=100*eps(T)) where {T<:Real}
    N=length(s)
    N==length(ds)||throw(DimensionMismatch("s and ds must have equal length"))
    Δs=s[2]-s[1]
    atol_s=rtol*max(one(T),abs(Δs))
    atol_ds=rtol*max(one(T),abs(ds[1]))
    return all(isapprox(s[i+1]-s[i],Δs;rtol=rtol,atol=atol_s) for i in 2:N-1)&&all(isapprox(ds[i],ds[1];rtol=rtol,atol=atol_ds) for i in 2:N)
end

"""
    _husimi_uniform_arclength(k::T,s::AbstractVector{T},ds::AbstractVector{T},u::AbstractVector{Num},L::T;c::Real=10.0,w::Real=7.0,full_p::Bool=false) where {T<:Real,Num<:Number} → Tuple{Matrix{T},Vector{T},Vector{T}}

Evaluate the Poincaré-Husimi function using the fast translating stencil for a
uniform physical-arclength boundary grid.

## Arguments
* `k::T`: Wavenumber.
* `s::AbstractVector{T}`: Uniform physical boundary arclength coordinates.
* `ds::AbstractVector{T}`: Uniform boundary quadrature weights.
* `u::AbstractVector{Num}`: Boundary-function values.
* `L::T`: Total boundary length.

## Keyword Arguments
* `c::Real=10.0`: Phase-space sampling density in units of the coherent-state width.
* `w::Real=7.0`: Gaussian truncation width in units of `σ=1/√k`.
* `full_p::Bool=false`: Whether to evaluate both signs of momentum explicitly.

## Returns
* `H::Matrix{T}`: Husimi-function matrix.
* `qs::Vector{T}`: Boundary-position coordinates.
* `ps::Vector{T}`: Full signed momentum coordinates.
"""
function _husimi_uniform_arclength(k::T,s::AbstractVector{T},ds::AbstractVector{T},u::AbstractVector{Num},L::T;c::Real=10.0,w::Real=7.0,full_p::Bool=false) where {T<:Real,Num<:Number}
    N=length(s)
    length(ds)==N==length(u)||throw(DimensionMismatch("s, ds and u must have equal length"))
    sig=inv(sqrt(k))
    x,dx,idx=_husimi_symmetric_window(s,ds,T(w)*sig)
    a=inv(T(2π)*sqrt(T(π)*k))
    uc=CircularVector(u)
    gauss=@. exp(-k*x^2/2)*dx
    gauss_l=@. exp(-k*(x+L)^2/2)*dx
    gauss_r=@. exp(-k*(x-L)^2/2)*dx
    np=max(1,ceil(Int,T(c)*sqrt(k)))
    ps=full_p ? collect(range(-one(T),one(T);length=2np+1)) : collect(range(zero(T),one(T);length=np+1))
    Δs=s[2]-s[1]
    q_stride=max(1,round(Int,(sig/T(c))/Δs))
    q_idx=collect(1:q_stride:N)
    qs=collect(s[q_idx])
    H=zeros(T,length(qs),length(ps))
    @fastmath for ip in eachindex(ps)
        p=ps[ip]
        cs=@. exp(-im*p*k*x)*gauss + exp(-im*p*k*(x+L))*gauss_l + exp(-im*p*k*(x-L))*gauss_r
        @inbounds for iq in eachindex(q_idx)
            j=q_idx[iq]
            h=sum(cs.*uc[j-idx+1:j+idx-1])
            H[iq,ip]=a*abs2(h)
        end
    end
    if !full_p
        H=hcat(reverse(H[:,2:end];dims=2),H)
        ps=antisym_vec(ps)
    end
    return H,qs,ps
end

"""
    husimi_function(k::T,s::AbstractVector{T},ds::AbstractVector{T},u::AbstractVector{Num},L::T,q::T,p::T;w::Real=4.0) where {T<:Real,Num<:Number} → T

Evaluate the Poincaré-Husimi function at a single phase-space point `(q,p)`.

The boundary integral is evaluated using the physical arclength coordinates
`s` and the supplied quadrature weights `ds`. Periodicity is handled by
extending the boundary data by one period in each direction.

## Arguments
* `k::T`: Wavenumber.
* `s::AbstractVector{T}`: Physical boundary arclength coordinates.
* `ds::AbstractVector{T}`: Boundary quadrature weights.
* `u::AbstractVector{Num}`: Boundary-function values.
* `L::T`: Total boundary length.
* `q::T`: Boundary-position coordinate.
* `p::T`: Momentum coordinate.

## Keyword Arguments
* `w::Real=4.0`: Gaussian truncation width in units of `σ`.

## Returns
* `H::T`: Husimi-function value at `(q,p)`.
"""
function husimi_function(k::T,s::AbstractVector{T},ds::AbstractVector{T},u::AbstractVector{Num},L::T,q::T,p::T;w::Real=4.0) where {T<:Real,Num<:Number}
    # original algorithm by Benjamin Batistić in python (https://github.com/clozej/quantum_billiards/blob/crt_public/src/CoreModules/HusimiFunctionsOld.py)
    length(s)==length(ds)==length(u)||throw(DimensionMismatch("s, ds and u must have equal length"))
    width=T(w)/sqrt(k)
    s_ext=vcat(s.-L,s,s.+L)
    u_ext=vcat(u,u,u)
    ds_ext=vcat(ds,ds,ds)
    q_ext=q+L
    lo=searchsortedfirst(s_ext,q_ext-width)
    hi=searchsortedlast(s_ext,q_ext+width)
    nf=sqrt(sqrt(k/π))
    sracc=zero(T)
    siacc=zero(T)
    @inbounds for j in lo:hi
        si=s_ext[j]-q_ext
        wt=nf*exp(-0.5*k*si*si)*ds_ext[j]
        θ=k*p*si
        s_,c_=sincos(θ)
        uj=u_ext[j]
        a=wt*real(uj)
        b=wt*imag(uj)
        sracc+=a*c_+b*s_
        siacc+=b*c_-a*s_
    end
    return (sracc*sracc+siacc*siacc)/(2*π*k) # not the actual normalization
end

"""
    husimi_function(k::T,s::AbstractVector{T},ds::AbstractVector{T},u::AbstractVector{Num},L::T,qs::AbstractVector{T},ps::AbstractVector{T};w::Real=4.0,full_p::Bool=false) where {T<:Real,Num<:Number} → Tuple{Matrix{T},AbstractVector{T},Vector{T}}

Evaluate the normalized Poincaré-Husimi function on prescribed boundary-position
and momentum grids.

The physical locations of the boundary nodes are supplied by `s`, while `ds`
contains the corresponding quadrature weights. If `full_p=false`, only the
nonnegative momentum grid is explicitly evaluated and the negative half is
reconstructed by reflection.

## Arguments
* `k::T`: Wavenumber.
* `s::AbstractVector{T}`: Physical boundary arclength coordinates.
* `ds::AbstractVector{T}`: Boundary quadrature weights.
* `u::AbstractVector{Num}`: Boundary-function values.
* `L::T`: Total boundary length.
* `qs::AbstractVector{T}`: Boundary-position grid.
* `ps::AbstractVector{T}`: Momentum grid.

## Keyword Arguments
* `w::Real=4.0`: Gaussian truncation width in units of `σ`.
* `full_p::Bool=false`: Whether `ps` already spans the full signed momentum interval.

## Returns
* `H::Matrix{T}`: Normalized Husimi-function matrix.
* `qs::AbstractVector{T}`: Boundary-position grid.
* `ps_out::Vector{T}`: Full signed momentum grid.
"""
function husimi_function(k::T,s::AbstractVector{T},ds::AbstractVector{T},u::AbstractVector{Num},L::T,qs::AbstractVector{T},ps::AbstractVector{T};w::Real=4.0,full_p::Bool=false) where {T<:Real,Num<:Number}
    length(s)==length(ds)==length(u)||throw(DimensionMismatch("s, ds and u must have equal length"))
    s_ext=vcat(s.-L,s,s.+L) # concatenate shifted copies so that any window centered at q can be sliced as a contiguous subarray without computing indices modulo L -> no CircularVector needed
    u_ext=vcat(u,u,u) # same for u
    ds_ext=vcat(ds,ds,ds) # same for ds
    nx=length(qs) # number of q grid points
    ny=length(ps) # number of p grid points
    Hp=zeros(T,ny,nx) # preallocate Husimi matrix (for p x q)
    nf=sqrt(sqrt(k/pi)) # normalization factor
    width=T(w)/sqrt(k) # Gaussian width (±wσ)
    # temporary vectors that will grow to the current window size and be reused for each q
    c_re=Vector{T}(undef,0) # buffer for real parts of coefficients
    c_im=Vector{T}(undef,0) # buffer for imaginary parts of coefficients
    si=Vector{T}(undef,0) # buffer for s differences (shifted arclengths (s−q))
    @inbounds for iq in 1:nx
        q=qs[iq]+L # Add +L so that the center q sits in the middle copy of s_ext for slicing
        lo=searchsortedfirst(s_ext,q-width) # find left index of window
        hi=searchsortedlast(s_ext,q+width) # find right index of window
        W=max(0,hi-lo+1) # size of the window indexwise
        if length(c_re)<W # binary-search the indices in s_ext that fall within [q−width, q+width] and resize buffers if needed (since they are reused for each iq)
            resize!(c_re,W)
            resize!(c_im,W)
            resize!(si,W)
        end
        @inbounds for t=0:W-1 # for each point in the window calculate weights and shifted arclengths
            j=lo+t # the index in the window (shifted)
            sdiff=s_ext[j]-q # shifted arclength with the above index
            si[t+1]=sdiff # store shifted arclength
            wt=nf*exp(-0.5*k*sdiff*sdiff)*ds_ext[j] # gaussian weight with quadrature 
            uj=u_ext[j] # corresponding boundary function value in that window
            if uj isa Real # hack to avoid complex multiplications when u is real-valued
                c_re[t+1]=wt*uj # real part of summand
                c_im[t+1]=zero(T) 
            else
                c_re[t+1]=wt*real(uj) # real & imag part of summand separately
                c_im[t+1]=wt*imag(uj)
            end 
        end
        @inbounds for ip in 1:ny # for each p grid point compute Husimi value at (q,p) using the precomputed buffers for this iq index
            kp=k*ps[ip] # k*p
            sracc=zero(T);siacc=zero(T) # real & imag accumulators for the sum
            @inbounds for t=1:W # sum over the window
                θ=kp*si[t] 
                s_,c_=sincos(θ) # s_=sin,c_=cos, the cheap way to compute these
                # (a+ib)(cos-i*sin)=(a*c+b*s)+i(b*c-a*s) to avoid complex multiplications temporaries
                a=c_re[t];b=c_im[t]
                sracc+=a*c_+b*s_
                siacc+=b*c_-a*s_
            end
            # Final Husimi value at (q,p) is the squared modulus of the sum, scaled by 1/(2πk) (removed scaling and just normalize in the end)
            Hp[ip,iq]=(sracc*sracc+siacc*siacc)
        end
    end
    if full_p
        H=permutedims(Hp)
        ps_out=collect(ps)
    else
        H=vcat(reverse(Hp[2:end,:];dims=1),Hp)|>permutedims
        ps_out=vcat(-reverse(ps[2:end]),ps)
    end
    H./=sum(H) # normalize it in the end since the 1/(2πk) normalization does not work well in practice with finite grids
    return H,qs,ps_out
end

"""
    husimi_function(k::T,pts::BoundaryPoints{T},u::AbstractVector{Num};c::Real=10.0,w::Real=7.0,full_p::Bool=false) where {T<:Real,Num<:Number} → Tuple{Matrix{T},Vector{T},Vector{T}}

Compute an automatically gridded Poincaré-Husimi function from a
[`BoundaryPoints`](@ref) discretization.

If `pts.s` and `pts.ds` define a uniform physical-arclength grid, the fast
translating-stencil implementation [`_husimi_uniform_arclength`](@ref) is used.
Otherwise an approximately equally resolved phase-space grid is constructed and
the general nonuniform-arclength quadrature is used.

The automatic grid spacing is approximately `σ/c`, where `σ=1/√k`. The
boundary-position grid is periodic and therefore contains `q=0` but not the
duplicate endpoint `q=L`.

## Arguments
* `k::T`: Wavenumber.
* `pts::BoundaryPoints{T}`: Full physical boundary discretization containing arclength coordinates `s` and quadrature weights `ds`.
* `u::AbstractVector{Num}`: Boundary-function values corresponding to `pts`.

## Keyword Arguments
* `c::Real=10.0`: Phase-space sampling density in units of the coherent-state width `σ=1/√k`.
* `w::Real=7.0`: Gaussian truncation width in units of `σ`.
* `full_p::Bool=false`: Whether to evaluate both signs of momentum explicitly. If `false`, only `p≥0` is evaluated and the negative half is reconstructed by reflection.

## Returns
* `H::Matrix{T}`: Poincaré-Husimi function on the returned phase-space grid.
* `qs::Vector{T}`: Periodic boundary-position grid.
* `ps::Vector{T}`: Full signed momentum grid.
"""
function husimi_function(k::T,pts::BoundaryPoints{T},u::AbstractVector{Num};c::Real=10.0,w::Real=7.0,full_p::Bool=false) where {T<:Real,Num<:Number}
    length(pts)==length(u)||throw(DimensionMismatch("pts and u must have equal length"))
    L=sum(pts.ds)
    if _husimi_uniform_arclength_grid(pts.s,pts.ds)
        return _husimi_uniform_arclength(k,pts.s,pts.ds,u,L;c=c,w=w,full_p=full_p)
    end
    sig=inv(sqrt(k))
    nq=max(2,ceil(Int,L*T(c)/sig))
    np=max(1,ceil(Int,T(c)/sig))
    qs=collect(range(zero(T),L;length=nq+1))[1:end-1]
    ps=full_p ? collect(range(-one(T),one(T);length=2np+1)) : collect(range(zero(T),one(T);length=np+1))
    return husimi_function(k,pts.s,pts.ds,u,L,qs,ps;w=w,full_p=full_p)
end

"""
    husimi_function(k::T,pts::BoundaryPoints{T},u::AbstractVector{Num},q::T,p::T;w::Real=4.0) where {T<:Real,Num<:Number} → T

Evaluate the Poincaré-Husimi function at one phase-space point `(q,p)` from a
[`BoundaryPoints`](@ref) discretization.

This overload uses the general physical-arclength quadrature and is valid for
both uniform and nonuniform boundary discretizations.

## Arguments
* `k::T`: Wavenumber.
* `pts::BoundaryPoints{T}`: Full physical boundary discretization containing arclength coordinates `s` and quadrature weights `ds`.
* `u::AbstractVector{Num}`: Boundary-function values corresponding to `pts`.
* `q::T`: Boundary-position coordinate.
* `p::T`: Momentum coordinate.

## Keyword Arguments
* `w::Real=4.0`: Gaussian truncation width in units of `σ=1/√k`.

## Returns
* `H::T`: Poincaré-Husimi value at `(q,p)`.
"""
function husimi_function(k::T,pts::BoundaryPoints{T},u::AbstractVector{Num},q::T,p::T;w::Real=4.0) where {T<:Real,Num<:Number}
    length(pts)==length(u)||throw(DimensionMismatch("pts and u must have equal length"))
    L=sum(pts.ds)
    return husimi_function(k,pts.s,pts.ds,u,L,q,p;w=w)
end

"""
    husimi_function(k::T,pts::BoundaryPoints{T},u::AbstractVector{Num},qs::AbstractVector{T},ps::AbstractVector{T};w::Real=4.0,full_p::Bool=false) where {T<:Real,Num<:Number} → Tuple{Matrix{T},AbstractVector{T},Vector{T}}

Evaluate the Poincaré-Husimi function on prescribed phase-space grids from a
[`BoundaryPoints`](@ref) discretization.

The general physical-arclength quadrature is used directly, so this overload is
valid for both uniform and nonuniform boundary discretizations.

## Arguments
* `k::T`: Wavenumber.
* `pts::BoundaryPoints{T}`: Full physical boundary discretization containing arclength coordinates `s` and quadrature weights `ds`.
* `u::AbstractVector{Num}`: Boundary-function values corresponding to `pts`.
* `qs::AbstractVector{T}`: Boundary-position coordinates at which to evaluate the Husimi function.
* `ps::AbstractVector{T}`: Momentum coordinates to evaluate explicitly.

## Keyword Arguments
* `w::Real=4.0`: Gaussian truncation width in units of `σ=1/√k`.
* `full_p::Bool=false`: Whether `ps` already spans the full signed momentum interval. If `false`, `ps` must contain the nonnegative half and the negative half is reconstructed by reflection.

## Returns
* `H::Matrix{T}`: Normalized Poincaré-Husimi function on the returned phase-space grid.
* `qs::AbstractVector{T}`: Input boundary-position grid.
* `ps_out::Vector{T}`: Full signed momentum grid.
"""
function husimi_function(k::T,pts::BoundaryPoints{T},u::AbstractVector{Num},qs::AbstractVector{T},ps::AbstractVector{T};w::Real=4.0,full_p::Bool=false) where {T<:Real,Num<:Number}
    length(pts)==length(u)||throw(DimensionMismatch("pts and u must have equal length"))
    L=sum(pts.ds)
    return husimi_function(k,pts.s,pts.ds,u,L,qs,ps;w=w,full_p=full_p)
end

"""
    husimi_function(ks::AbstractVector{T},vec_us::AbstractVector{<:AbstractVector{<:Number}},vec_pts::AbstractVector{<:BoundaryPoints{T}};c::Real=10.0,w::Real=7.0,full_p::Bool=false) where {T<:Real} → Tuple

Construct automatically gridded Poincaré-Husimi functions for a collection of
states.

Each state is dispatched independently according to its boundary discretization.
Uniform physical-arclength grids use the fast translating-stencil algorithm,
while nonuniform grids use the general physical-arclength quadrature.

## Arguments
* `ks::AbstractVector{T}`: Eigenwavenumbers.
* `vec_us::AbstractVector{<:AbstractVector{<:Number}}`: Boundary functions, one vector per state.
* `vec_pts::AbstractVector{<:BoundaryPoints{T}}`: Boundary discretizations, one per state.

## Keyword Arguments
* `c::Real=10.0`: Phase-space sampling density in units of the coherent-state width `σ=1/√k`.
* `w::Real=7.0`: Gaussian truncation width in units of `σ`.
* `full_p::Bool=false`: Whether to evaluate both signs of momentum explicitly.

## Returns
* `Hs_return::Vector{Matrix}`: Poincaré-Husimi matrices for successfully processed states.
* `ps_return::Vector{Vector}`: Full signed momentum grids corresponding to `Hs_return`.
* `qs_return::Vector{Vector}`: Boundary-position grids corresponding to `Hs_return`.
"""
function husimi_function(ks::AbstractVector{T},vec_us::AbstractVector{<:AbstractVector{<:Number}},vec_pts::AbstractVector{<:BoundaryPoints{T}};c::Real=10.0,w::Real=7.0,full_p::Bool=false) where {T<:Real}
    length(ks)==length(vec_us)==length(vec_pts)||throw(DimensionMismatch("ks, vec_us and vec_pts must have equal length"))
    valid_indices=trues(length(ks))
    Hs_return=Vector{Matrix}(undef,length(ks))
    ps_return=Vector{Vector}(undef,length(ks))
    qs_return=Vector{Vector}(undef,length(ks))
    p=Progress(length(ks);desc="Constructing Husimi matrices, N=$(length(ks))")
    Threads.@threads for i in eachindex(ks)
        try
            H,qs,ps=husimi_function(ks[i],vec_pts[i],vec_us[i];c=c,w=w,full_p=full_p)
            Hs_return[i]=H
            ps_return[i]=ps
            qs_return[i]=qs
        catch e
            @debug "Husimi fail at k=$(ks[i])" exception=(e,catch_backtrace())
            valid_indices[i]=false
        end
        next!(p)
    end
    return Hs_return[valid_indices],ps_return[valid_indices],qs_return[valid_indices]
end

"""
    husimi_function(ks::AbstractVector{T},vec_us::AbstractVector{<:AbstractVector{<:Number}},vec_pts::AbstractVector{<:BoundaryPoints{T}},nx::Integer,ny::Integer;w::Real=4.0,full_p::Bool=false) where {T<:Real} → Tuple

Construct fixed-grid Poincaré-Husimi functions for a collection of states.

All states are evaluated on a common periodic boundary-position grid containing
`nx` points on `[0,L)` and a common momentum grid. The general
physical-arclength quadrature is used, so the underlying boundary
discretizations may be nonuniform.

The perimeter estimate is taken from the boundary discretization corresponding
to the largest wavenumber, which normally has the highest boundary resolution.

## Arguments
* `ks::AbstractVector{T}`: Eigenwavenumbers.
* `vec_us::AbstractVector{<:AbstractVector{<:Number}}`: Boundary functions, one vector per state.
* `vec_pts::AbstractVector{<:BoundaryPoints{T}}`: Boundary discretizations, one per state.
* `nx::Integer`: Number of periodic boundary-position grid points.
* `ny::Integer`: Requested number of momentum grid points. For `full_p=false`, reflection through `p=0` produces an odd signed grid and therefore the returned count is `2cld(ny,2)-1`.

## Keyword Arguments
* `w::Real=4.0`: Gaussian truncation width in units of `σ=1/√k`.
* `full_p::Bool=false`: Whether to evaluate the full signed momentum interval directly.

## Returns
* `Hs::Vector{Matrix{T}}`: Poincaré-Husimi matrices for successfully processed states.
* `ps_return::Vector{Vector{T}}`: Full signed momentum grids corresponding to `Hs`.
* `qs_return::Vector{Vector{T}}`: Common periodic boundary-position grids corresponding to `Hs`.
"""
function husimi_function(ks::AbstractVector{T},vec_us::AbstractVector{<:AbstractVector{<:Number}},vec_pts::AbstractVector{<:BoundaryPoints{T}},nx::Integer,ny::Integer;w::Real=4.0,full_p::Bool=false) where {T<:Real}
    length(ks)==length(vec_us)==length(vec_pts)||throw(DimensionMismatch("ks, vec_us and vec_pts must have equal length"))
    isempty(ks)&&return Matrix{T}[],Vector{Vector{T}}(),Vector{Vector{T}}()
    imax=argmax(ks)
    L=sum(vec_pts[imax].ds)
    qs=collect(range(zero(T),L;length=nx+1))[1:end-1]
    ps=full_p ? collect(range(-one(T),one(T);length=ny)) : collect(range(zero(T),one(T);length=cld(ny,2)))
    Hs=Vector{Matrix{T}}(undef,length(ks))
    ok=trues(length(ks))
    pbar=Progress(length(ks);desc="Husimi N=$(length(ks))")
    Threads.@threads for i in eachindex(ks)
        try
            H,_,_=husimi_function(ks[i],vec_pts[i],vec_us[i],qs,ps;w=w,full_p=full_p)
            Hs[i]=H
        catch e
            @debug "Husimi fail at k=$(ks[i])" exception=(e,catch_backtrace())
            ok[i]=false
        end
        next!(pbar)
    end
    ps_out=full_p ? ps : vcat(-reverse(ps[2:end]),ps)
    n=count(ok)
    return Hs[ok],[copy(ps_out) for _ in 1:n],[copy(qs) for _ in 1:n]
end