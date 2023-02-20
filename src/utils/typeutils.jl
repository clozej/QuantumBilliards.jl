
function set_precision(a)
    #expand for other types of numbers
    t = typeof(a)
    return t == Float32 ? Float32(1e-8) : convert(t,1e-16) 
end