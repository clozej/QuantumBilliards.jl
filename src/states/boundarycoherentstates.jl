include("../abstracttypes.jl")

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
