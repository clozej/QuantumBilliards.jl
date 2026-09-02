################# GENERAL ####################

@inline H(n::Int,x::T) where {T<:Real}=Bessels.hankelh1(n,x)
@inline H(n::Int,x::Complex{T}) where {T<:Real}=SpecialFunctions.besselh(n,1,x)
@inline function hankel_pair01(x);h0=H(0,x);h1=H(1,x);return h0,h1;end
@inline Φ_helmholtz(k::T,r::T) where {T<:Real}=(im/4)*Bessels.hankelh1(0,k*r) # used for CFIE solvers to get the normal derivative of the wavefunctions
# this is the slp kernel used in the hypersingular maue kress formula for the normal derivative of the DLP kernel.