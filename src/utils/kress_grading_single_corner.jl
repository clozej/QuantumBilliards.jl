#  Reference:
#  R. Kress, "Boundary integral equations in time-harmonic acoustic scattering",
#  Math. Comput. Modelling 15 (1991), 229–243.
const TWO_PI=2*pi

# NOTE: This will still probably cluster 2 pts below machine precision to the corner, but no problem since they
# have jacobian weight also below machine precision. When practically doing computations they have not bearing to
# the accuracy of the result. One can even safely take adjoints ds[j]/ds[i] at this point since it gives finite values.

# v(s,q)
@inline function _kress_v(s::T,q::T) where {T<:Real}
    x=(pi-s)/pi
    return (inv(q)-T(0.5))*x^3+inv(q)*((s-pi)/pi)+T(0.5) # (1 / q - 1 / 2) * ((pi - s) / pi)^3 + 1 / q * (s - pi) / pi + 1 / 2
end
# dv(s,q)/ds
@inline function _kress_vprime(s::T, q::T) where {T<:Real}
    x=(pi-s)/pi
    return -(T(3)/pi)*(inv(q)-T(0.5))*x^2+inv(q)/pi # d/ds x^3 = -3 x^2 / pi
end
@inline function _kress_vdoubleprime(s::T,q::T) where {T<:Real}
    x=(pi-s)/pi
    return T((6/(pi^2))*(inv(q)-0.5)*x)
end
# w(s,q) = v(s,q)^q / (v(s,q)^q + v(2π-s,q)^q)
@inline function _kress_w(s::T,q::T) where {T<:Real}
    a=_kress_v(s,q)^q
    b=_kress_v(TWO_PI-s,q)^q
    return TWO_PI*a/(a+b)
end
# dw(s,q)/ds = 2π * (dv(s,q)/ds * (v(s,q)^q + v(2π-s,q)^q) - v(s,q)^q * (dv(s,q)/ds + dv(2π-s,q)/ds)) / (v(s,q)^q + v(2π-s,q)^q)^2
@inline function _kress_wprime(s::T,q::T) where {T<:Real}
    vs=_kress_v(s,q)
    vsp=_kress_vprime(s,q)
    vt=_kress_v(TWO_PI-s,q)
    # d/ds v(2π-s) = -v'(2π-s)
    vtp=-_kress_vprime(TWO_PI-s,q)
    a=vs^q
    b=vt^q
    ap=q*vs^(q-1)*vsp
    bp=q*vt^(q-1)*vtp
    den=a+b
    return TWO_PI*(ap*den-a*(ap+bp))/den^2
end
@inline function _kress_wdoubleprime(s::T,q::T) where {T<:Real}
    vs=_kress_v(s,q)
    vsp=_kress_vprime(s,q)
    vspp=_kress_vdoubleprime(s,q)
    vt=_kress_v(TWO_PI-s,q)
    vtp=-_kress_vprime(TWO_PI-s,q)
    vtpp=_kress_vdoubleprime(TWO_PI-s,q)
    a=vs^q
    b=vt^q
    ap=q*vs^(q-1)*vsp
    bp=q*vt^(q-1)*vtp
    app=q*((q-one(T))*vs^(q-2)*vsp^2+vs^(q-1)*vspp)
    bpp=q*((q-one(T))*vt^(q-2)*vtp^2+vt^(q-1)*vtpp)
    den=a+b
    denp=ap+bp
    num=ap*b-a*bp
    nump=app*b-a*bpp
    return T(TWO_PI*(nump*den-2*num*denp)/den^3)
end

# helper to determine the minimal spacings, if too small it causes issues with solvers and 
# especially desymmetrization. Best to have a threashold for smallest distance, say 1e-12
@inline function _min_periodic_spacing_sorted(xs::Vector{T}) where {T<:Real}
    N=length(xs)
    N<=1 && return typemax(T)
    dmin=typemax(T)
    @inbounds for i in 1:N-1
        dmin=min(dmin,xs[i+1]-xs[i])
    end
    return min(dmin,T(TWO_PI)+xs[1]-xs[end])
end

# Computes the Kress graded nodes and weights for a given number of points `N`
# and grading parameter `q`, using a periodic grid that works for both even and
# odd N.
function kress_graded_nodes_data(::Type{T},N::Int;q=3,minsep_tol=1e-12) where {T<:Real}
    qT=T(q)
    qT>one(T) || error("Require q>1 for Kress grading.")
    h=T(TWO_PI)/T(N)
    δ=h/T(2)
    while qT>one(T)
        σ=Vector{T}(undef,N)
        s=Vector{T}(undef,N)
        a=Vector{T}(undef,N)
        a2=Vector{T}(undef,N)
        wq=Vector{T}(undef,N)
        @inbounds for k in 1:N
            σ[k]=δ+T(k-1)*h
            s[k]=_kress_w(σ[k],qT)
            a[k]=_kress_wprime(σ[k],qT)
            a2[k]=_kress_wdoubleprime(σ[k],qT)
            wq[k]=h*a[k]
        end
        minsep=_min_periodic_spacing_sorted(s)
        minsep>=minsep_tol && return σ,s,a,a2,wq
        qnew=max(one(T),qT*0.9) # multiply by 0.9 each time to reduce it until it gives larger than min separation
        @warn "Kress grading nodes too close; reducing q." q_old=qT q_new=qnew minsep=minsep minsep_tol=minsep_tol N=N
        qT=qnew
    end
    error("Kress grading is impossible: q reached 1 while min periodic spacing stayed below minsep_tol=$(minsep_tol).")
end
