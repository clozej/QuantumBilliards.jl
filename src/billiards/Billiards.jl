module Billiards

include("geometry.jl")
include("stadium.jl")

include("triangle.jl")
export adapt_basis
#include("limacon.jl")
#include("rectangle.jl")
export Stadium, Triangle
export curve, tangent, normal, arc_length
export tangent_vec, normal_vec
#convenience functions may be moved somewhere else
export make_stadium_and_basis, make_triangle_and_basis 
end