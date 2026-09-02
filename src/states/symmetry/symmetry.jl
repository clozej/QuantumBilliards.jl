# Estimate the minimum nonzero and maximum interaction radii on the complete physical boundary.
function estimate_rmin_rmax(pts::BoundaryPoints{T},::Nothing) where {T<:Real}
    N=length(pts)
    xy=pts.xy
    tol2=eps(T)^2
    rmin=Inf
    rmax=0.0
    @inbounds for j in 2:N
        xj=xy[j]
        for i in 1:j-1
            xi=xy[i]
            dx=xi[1]-xj[1]
            dy=xi[2]-xj[2]
            d2=muladd(dx,dx,dy*dy)
            d2<=tol2&&continue
            r=sqrt(Float64(d2))
            rmin=min(rmin,r)
            rmax=max(rmax,r)
        end
    end
    isfinite(rmin)&&rmax>0.0||throw(ArgumentError("Unable to determine radial interval"))
    return rmin,rmax
end

# Estimate interaction radii for the exact symmetry-reduced operator using its full-boundary image interactions.
function estimate_rmin_rmax(pts::BoundaryPoints{T},symmetry) where {T<:Real}
    isnothing(symmetry)&&return estimate_rmin_rmax(pts,nothing)
    orbits=symmetry_index_orbits(T,pts,symmetry)
    Ifund=orbits.Ifund
    m=fundamental_size(orbits)
    ng=orbit_size(orbits)
    xy=pts.xy
    tol2=eps(T)^2
    rmin=Inf
    rmax=0.0
    @inbounds for b in 1:m
        j=Ifund[b]
        for a in 1:m
            i=Ifund[a]
            xi=xy[i]
            for l in 1:ng
                q=orbits.fund_to_full[l,b]
                l==1&&i==j&&continue
                xq=xy[q]
                dx=xi[1]-xq[1]
                dy=xi[2]-xq[2]
                d2=muladd(dx,dx,dy*dy)
                d2<=tol2&&continue
                r=sqrt(Float64(d2))
                rmin=min(rmin,r)
                rmax=max(rmax,r)
            end
        end
    end
    isfinite(rmin)&&rmax>0.0||throw(ArgumentError("Unable to determine radial interval"))
    return rmin,rmax
end