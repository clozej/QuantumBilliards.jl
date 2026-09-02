@inline midpoints(v::AbstractVector)=(v[1:end-1].+v[2:end])./2
@inline midpoints(r::AbstractRange)=r[1:end-1].+step(r)/2

function median(v::AbstractVector{T}) where {T<:Real}
    n=length(v)
    w=copy(v)
    if isodd(n)
        return partialsort!(w,(n+1)÷2)
    else
        k=n÷2
        m1=partialsort!(w,k)
        m2=partialsort!(w,k+1)
        return (m1+m2)/2
    end
end

function quantile(v::AbstractVector{T},p::Real) where {T<:Real}
    n=length(v);n==1&&return v[firstindex(v)]
    h=1+(n-1)*p;k=floor(Int,h);γ=h-k
    w=copy(v);x0=partialsort!(w,k)
    k==n&&return x0
    x1=partialsort!(w,k+1)
    return x0+γ*(x1-x0)
end

# Batch n objects into batches of size batch_size
@inline function _nbatches(n::Int,batch_size::Int)
    return cld(n,batch_size)
end
# Get the first index of the b-th batch
@inline _batch_first(b::Int,batch_size::Int)=1+(b-1)*batch_size
# Get the last index of the b-th batch, ensuring it does not exceed n
@inline _batch_last(b::Int,batch_size::Int,n::Int)=min(b*batch_size,n)

function pad_limits(xlim,ylim;padding=0.01)
    return (xlim[1]-padding,xlim[2]+padding),(ylim[1]-padding,ylim[2]+padding)
end
function boundary_limits(curves;grd=1000,padding=0.01)
    x_bnd=Vector{Any}()
    y_bnd=Vector{Any}()
    for crv in curves
        L=crv.length
        N_bnd=max(512,round(Int,grd/L))
        t=range(0.0,1.0,N_bnd)[1:end-1]
        pts=BilliardGeometry.curve(crv,t)
        append!(x_bnd,getindex.(pts,1))
        append!(y_bnd,getindex.(pts,2))
    end
    x_bnd[end]=x_bnd[1]
    y_bnd[end]=y_bnd[1]
    xlim=extrema(x_bnd)
    ylim=extrema(y_bnd)
    return pad_limits(xlim,ylim;padding=padding)
end

"""
    random_interior_points(billiard::BilliardGeometry.AbsBilliard, N::Int; grd::Int = 1000)

Retrieves the bounding limits of the billiard’s fundamental boundary (using boundary_limits) and generates points within the bounds.
1st Checks if each point is inside the billiard using is_inside then Continues until N valid interior points are collected.

# Arguments
- `billiard::BilliardGeometry.AbsBilliard`: Instance of the geometry so we can check what is the interior.
- `N::Int`: The number of internal points we want.
- `grd::Int=1000`: Parameter that determines the precision of the limits of the billiard boundary. Usually 1000 is enough and there is no need to change.

# Returns:
- `pts::Vector{SVector{2,Float64}}`: A vector of points inside the billiard.
"""
function random_interior_points(billiard::BilliardGeometry.AbsBilliard,N::Int;grd::Int=1000)
    xlim,ylim=boundary_limits(BilliardGeometry.get_boundary_curves(billiard);grd=grd)
    dx=xlim[2]-xlim[1]
    dy=ylim[2]-ylim[1]
    pts=[]
    while length(pts)<N
        x=(dx.*rand().+xlim[1]) 
        y=(dy.*rand().+ylim[1])
        pt=SVector(x,y)
        if BilliardGeometry.is_inside(billiard,[pt])[1] #TODO Kind of stupid that we have to access 1 element vector b/c there is no single vector implementation of is_inside
            push!(pts,pt)
        end
    end
    return pts
end