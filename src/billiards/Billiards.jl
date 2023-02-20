module Billiards

include("geometry.jl")
include("stadium.jl")

include("triangle.jl")
#include("limacon.jl")
#include("rectangle.jl")
export Stadium, Triangle

#convenience functions may be moved somewhere else
export make_stadium_and_basis, make_triangle_and_basis 
end