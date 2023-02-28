#include("../abstracttypes.jl")

using CircularArrays

function antisym_vec(x)
    v = reverse(-x[2:end])
    return append!(v,x)
end

function husimi_function(k,u,s; c = 10.0, w = 7.0)
    #compute coherrent state weights
    N = length(s)
    sig = one(k)/sqrt(k) #width of the gaussian
    x = s[s.<=w*sig]
    idx = length(x) #do not change order here
    x = antisym_vec(x)
    a = one(k)/(2*pi*sqrt(pi*k)) #normalization factor in this version Hsimi is not noramlized to 1
    ds = (x[end]-x[1])/length(x) #integration weigth
    uc = CircularVector(u) #allows circular indexing
    gauss = @. exp(-k/2*x^2)*ds
    #construct evaluation points in p coordinate
    ps = collect(range(0.0,1.0,step = sig/c))
    #construct evaluation points in q coordinate
    q_stride = length(s[s.<=sig/c])
    q_idx = collect(1:q_stride:N)
    push!(q_idx,N) #add last point
    #println(length(q_idx))
    #println(q_idx)
    #println(N)
    qs = s[q_idx]
    #println(length(qs))
    H = zeros(typeof(k),length(qs),length(ps))
    for i in eachindex(ps)   
        cr = @. cos(ps[i]*k*x)*gauss #real part of coherent state
        ci = @. sin(ps[i]*k*x)*gauss #imag part of coherent state
        for j in eachindex(q_idx)
            u_w = uc[q_idx[j]-idx+1:q_idx[j]+idx-1] #window with relevant values of u
            hr = sum(cr.*u_w)
            hi = sum(ci.*u_w)
            H[j,i] = a*(hr*hr + hi*hi)
        end
    end
    return H, qs, ps    
end

function husimi_function(state::S, basis::Ba, billiard::Bi;  b = 5.0, c = 10.0, w = 7.0, sampler=fourier_nodes, include_virtual=true) where {S<:AbsState,Ba<:AbsBasis,Bi<:AbsBilliard}
    k = state.k
    u, s = boundary_function(state, basis, billiard; b=b, sampler=sampler, include_virtual=include_virtual)
    return husimi_function(k,u,s; c = c, w = w)
end
#=
function coherent(q,p,k,s,L,m::Int)
    let x = s-q+m*L
        a = (k/pi)^0.25
        ft = exp(im*p*k*x) 
        gauss = exp(-k/2*x^2)
        return a*ft*gauss
    end 
end

function coherent(q,p,k,s,L,m::Int,b::Complex)
    let x = s-q+m*L
        a = (k*imag(b)/pi)^0.25
        ft = exp(im*p*k*x) 
        fb = exp(im*real(b)/2*k*x^2)
        gauss = exp(-imag(b)/2*k*x^2)
        return a*ft*fb*gauss
    end 
end
=#