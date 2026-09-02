###############################################################################
# Core Chebyshev routines for piecewise spectral approximation.
#
# This file provides the low–level building blocks used throughout the BIM /
# Beyn machinery to approximate Bessel functions
# on many geometric panels with high accuracy and stability.
#
# Contents
# --------
# 1.  _cheb_clenshaw(c,t)
#       Stable evaluation of a Chebyshev series
#           f(t) = ∑_{m=0}^M c[m] T_m(t)
#     using Clenshaw’s reverse recurrence.  This routine is called tens of
#     millions of times inside matrix–vector products, so it is optimized for
#     zero allocations and minimal branches.
#
# 2.  _chebfit!(c,f)
#       Compute Chebyshev coefficients c from function samples f evaluated at
#       Chebyshev–Lobatto nodes t_j = cos(π j / M).  This is a direct O(M²)
#       cosine transform (DCT-I) with endpoint half-weights:
#           c_m = (2/M) * ∑ w_j f_j cos(π j m / M),  w_0=w_M=1/2.
#       After fitting, c[1] is halved so that the normalization matches
#       _cheb_clenshaw.
#
# 3.  _breaks_uniform(rmin,rmax,np)
#       Generate np uniformly sized panels, returning np+1 breakpoints.
#
# Usage Pattern
# -------------
# These routines are intentionally minimal: they know nothing about the
# underlying special function or geometry.  High-level modules call these primitives to:
#
#   • evaluate f(r) at Chebyshev nodes on each panel;
#   • use _chebfit! to obtain Chebyshev coefficients;
#   • evaluate the function later with _cheb_clenshaw;
#   • divide the domain using _breaks_uniform depending on
#     the expected oscillatory behavior.
#
# All functions are allocation-free inside tight loops and suitable for
# multi-threaded vectorized usage.
###############################################################################

# =============================================================================
# Compute Chebyshev expansion coefficients c (in place) from samples f taken at
# Chebyshev–Lobatto nodes on [-1,1]. This is a direct O(M^2) Discrete Cosine
# Transform of type I with endpoint half-weights.
#
# Inputs
#   c :: Vector{ComplexF64}   # output coefficients (length M+1), filled in place
#   f :: Vector{ComplexF64}   # samples at Lobatto nodes, length M+1 with t_j = cos(π j / M), j=0..M
#
# Construction
#   Given f_j = f(t_j), this computes coefficients c_m (m=0..M) such that
#       f(t) ≈ ∑_{m=0}^M c_m T_m(t)
#   via the weighted projection
#       c_m = (2/M) * ∑_{j=0}^M w_j * f_j * cos(π j m / M),
#   with Lobatto endpoint weights
#       w_0 = w_M = 1/2,   and   w_j = 1 for j=1..M-1.
# =============================================================================
@inline function _chebfit!(c::Vector{ComplexF64},f::Vector{ComplexF64})::Vector{ComplexF64}
    M=length(f)-1
    @inbounds for m in 0:M
        s=0.0+0.0im
        s+=0.5*f[1] # endpoint j = 0: cospi(0) = 1, w0 = 0.5
        s+=0.5*((isodd(m) ? -1.0 : 1.0)*f[M+1]) # endpoint j = M: cospi(m) = (-1)^m
        for j in 1:M-1
            s+=f[j+1]*cospi((j*m)/M)
        end
        c[m+1]=(2/M)*s
    end
    c[1]*=0.5
    c[end]*=0.5
    return c
end

# =============================================================================
# Generate uniformly spaced panel breakpoints between rmin and rmax.
# The interval [rmin, rmax] is divided into np equal panels, returning
# np+1 breakpoints suitable for piecewise-Chebyshev interpolation.
#
# Inputs
#   rmin :: Float64    # lower bound (must be positive)
#   rmax :: Float64    # upper bound (must be > rmin)
#   np   :: Int        # number of panels
#
# Output
#   b :: Vector{Float64}   # length np+1, with b[1]=rmin, b[end]=rmax, and uniform spacing h = (rmax - rmin)/np
# =============================================================================
@inline function _breaks_uniform(rmin::Float64,rmax::Float64,np::Int)::Vector{Float64}
    h=(rmax-rmin)/np
    b=Vector{Float64}(undef,np+1)
    @inbounds for i in 0:np
        b[i+1]=muladd(i,h,rmin)
    end
    b[end]=rmax
    return b
end

# =============================================================================
# UNROLLED CLENSHAW
# =============================================================================

#=
Unrolled versions of the Clenshaw recurrence for small fixed M (4 to 8) to avoid loop overhead. These are used in the innermost loops of Hankel evaluations where M is typically small with large panelization.
=#

#4th degree
@inline function _cheb_clenshaw_4(coeffs::AbstractVector{ComplexF64},t::Float64)::ComplexF64
    b1=0.0+0.0im;b2=0.0+0.0im;u=2*t
    b0=muladd(u,b1,-b2)+coeffs[5];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[4];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[3];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[2];b2=b1;b1=b0
    return muladd(t,b1,-b2)+coeffs[1]
end
@inline function _cheb_clenshaw_4(coeffs::AbstractMatrix{ComplexF64},col::Int,t::Float64)::ComplexF64
    b1=0.0+0.0im;b2=0.0+0.0im;u=2*t
    b0=muladd(u,b1,-b2)+coeffs[5,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[4,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[3,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[2,col];b2=b1;b1=b0
    return muladd(t,b1,-b2)+coeffs[1,col]
end

# 5th degree
@inline function _cheb_clenshaw_5(coeffs::AbstractVector{ComplexF64},t::Float64)::ComplexF64
    b1=0.0+0.0im;b2=0.0+0.0im;u=2*t
    b0=muladd(u,b1,-b2)+coeffs[6];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[5];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[4];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[3];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[2];b2=b1;b1=b0
    return muladd(t,b1,-b2)+coeffs[1]
end
@inline function _cheb_clenshaw_5(coeffs::AbstractMatrix{ComplexF64},col::Int,t::Float64)::ComplexF64
    b1=0.0+0.0im;b2=0.0+0.0im;u=2*t
    b0=muladd(u,b1,-b2)+coeffs[6,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[5,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[4,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[3,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[2,col];b2=b1;b1=b0
    return muladd(t,b1,-b2)+coeffs[1,col]
end

# 6th degree
@inline function _cheb_clenshaw_6(coeffs::AbstractVector{ComplexF64},t::Float64)::ComplexF64
    b1=0.0+0.0im;b2=0.0+0.0im;u=2*t
    b0=muladd(u,b1,-b2)+coeffs[7];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[6];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[5];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[4];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[3];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[2];b2=b1;b1=b0
    return muladd(t,b1,-b2)+coeffs[1]
end
@inline function _cheb_clenshaw_6(coeffs::AbstractMatrix{ComplexF64},col::Int,t::Float64)::ComplexF64
    b1=0.0+0.0im;b2=0.0+0.0im;u=2*t
    b0=muladd(u,b1,-b2)+coeffs[7,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[6,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[5,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[4,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[3,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[2,col];b2=b1;b1=b0
    return muladd(t,b1,-b2)+coeffs[1,col]
end

# 7th degree
@inline function _cheb_clenshaw_7(coeffs::AbstractVector{ComplexF64},t::Float64)::ComplexF64
    b1=0.0+0.0im;b2=0.0+0.0im;u=2*t
    b0=muladd(u,b1,-b2)+coeffs[8];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[7];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[6];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[5];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[4];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[3];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[2];b2=b1;b1=b0
    return muladd(t,b1,-b2)+coeffs[1]
end
@inline function _cheb_clenshaw_7(coeffs::AbstractMatrix{ComplexF64},col::Int,t::Float64)::ComplexF64
    b1=0.0+0.0im;b2=0.0+0.0im;u=2*t
    b0=muladd(u,b1,-b2)+coeffs[8,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[7,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[6,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[5,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[4,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[3,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[2,col];b2=b1;b1=b0
    return muladd(t,b1,-b2)+coeffs[1,col]
end

# 8th degree
@inline function _cheb_clenshaw_8(coeffs::AbstractVector{ComplexF64},t::Float64)::ComplexF64
    b1=0.0+0.0im;b2=0.0+0.0im;u=2*t
    b0=muladd(u,b1,-b2)+coeffs[9];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[8];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[7];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[6];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[5];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[4];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[3];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[2];b2=b1;b1=b0
    return muladd(t,b1,-b2)+coeffs[1]
end
@inline function _cheb_clenshaw_8(coeffs::AbstractMatrix{ComplexF64},col::Int,t::Float64)::ComplexF64
    b1=0.0+0.0im;b2=0.0+0.0im;u=2*t
    b0=muladd(u,b1,-b2)+coeffs[9,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[8,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[7,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[6,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[5,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[4,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[3,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[2,col];b2=b1;b1=b0
    return muladd(t,b1,-b2)+coeffs[1,col]
end

# 9th degree
@inline function _cheb_clenshaw_9(coeffs::AbstractVector{ComplexF64},t::Float64)::ComplexF64
    b1=0.0+0.0im;b2=0.0+0.0im;u=2*t
    b0=muladd(u,b1,-b2)+coeffs[10];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[9];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[8];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[7];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[6];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[5];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[4];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[3];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[2];b2=b1;b1=b0
    return muladd(t,b1,-b2)+coeffs[1]
end
@inline function _cheb_clenshaw_9(coeffs::AbstractMatrix{ComplexF64},col::Int,t::Float64)::ComplexF64
    b1=0.0+0.0im;b2=0.0+0.0im;u=2*t
    b0=muladd(u,b1,-b2)+coeffs[10,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[9,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[8,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[7,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[6,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[5,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[4,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[3,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[2,col];b2=b1;b1=b0
    return muladd(t,b1,-b2)+coeffs[1,col]
end

# 10th degree
@inline function _cheb_clenshaw_10(coeffs::AbstractVector{ComplexF64},t::Float64)::ComplexF64
    b1=0.0+0.0im;b2=0.0+0.0im;u=2*t
    b0=muladd(u,b1,-b2)+coeffs[11];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[10];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[9];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[8];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[7];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[6];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[5];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[4];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[3];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[2];b2=b1;b1=b0
    return muladd(t,b1,-b2)+coeffs[1]
end
@inline function _cheb_clenshaw_10(coeffs::AbstractMatrix{ComplexF64},col::Int,t::Float64)::ComplexF64
    b1=0.0+0.0im;b2=0.0+0.0im;u=2*t
    b0=muladd(u,b1,-b2)+coeffs[11,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[10,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[9,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[8,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[7,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[6,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[5,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[4,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[3,col];b2=b1;b1=b0
    b0=muladd(u,b1,-b2)+coeffs[2,col];b2=b1;b1=b0
    return muladd(t,b1,-b2)+coeffs[1,col]
end

# =============================================================================
# TOP CALLING FUNCTIONS FOR CHEBYSHEV
# =============================================================================

# We evaluated the function at Chebyshev–Lobatto nodes at each panel and then fit the Chebyshev coefficients. c is the vector of Chebyshev coefficients and f is the vector of function evaluations at Chebyshev nodes
# This is all done with clenshaw's reverse recurrence algorithm for numerical stability.
@inline function _cheb_clenshaw(c::AbstractVector{ComplexF64},t::Float64)::ComplexF64
    M=length(c)-1
    if M==4
        return _cheb_clenshaw_4(c,t)
    elseif M==5
        return _cheb_clenshaw_5(c,t)
    elseif M==6
        return _cheb_clenshaw_6(c,t)
    elseif M==7
        return _cheb_clenshaw_7(c,t)
    elseif M==8
        return _cheb_clenshaw_8(c,t)
    elseif M==9
        return _cheb_clenshaw_9(c,t)
    elseif M==10
        return _cheb_clenshaw_10(c,t)
    else
        b1=0.0+0.0im
        b2=0.0+0.0im
        u=2*t
        @inbounds for k in M:-1:1
            b=muladd(u,b1,c[k+1]-b2)
            b2=b1
            b1=b
        end
        return muladd(t,b1,c[1])-b2
    end
end

# Matrix of coefficients is stored in column-major order, so we pass the column index to these functions when we have a set of  coefficients for multiple panels. The unrolled versions are used when M is small (4 to 8) to avoid loop overhead in the innermost loops of Hankel evaluations where M is typically small with large panelization.
@inline function _cheb_clenshaw_col(coeffs::AbstractMatrix{ComplexF64},col::Int,M::Int,t::Float64)::ComplexF64
    if M==4
        return _cheb_clenshaw_4(coeffs,col,t)
    elseif M==5
        return _cheb_clenshaw_5(coeffs,col,t)
    elseif M==6
        return _cheb_clenshaw_6(coeffs,col,t)
    elseif M==7
        return _cheb_clenshaw_7(coeffs,col,t)
    elseif M==8
        return _cheb_clenshaw_8(coeffs,col,t)
    elseif M==9
        return _cheb_clenshaw_9(coeffs,col,t)
    elseif M==10
        return _cheb_clenshaw_10(coeffs,col,t)
    else
        b1=0.0+0.0im
        b2=0.0+0.0im
        u=2*t
        @inbounds for k in M:-1:1
            b0=muladd(u,b1,coeffs[k+1,col]-b2)
            b2=b1
            b1=b0
        end
        return muladd(t,b1,coeffs[1,col])-b2
    end
end

