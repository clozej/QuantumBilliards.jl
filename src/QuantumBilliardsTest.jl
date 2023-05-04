
#using Reexport

#abstract types
include("abstracttypes.jl")
#utils must be included here so modules work
#export AbsBasis
include("utils/coordinatesystems.jl")
include("utils/geometryutils.jl")
include("utils/billiardutils.jl")
include("utils/typeutils.jl")
include("utils/gridutils.jl")
include("utils/symmetry.jl")

#solvers
#include("solvers/Solvers.jl")
#@reexport using .Solvers
include("solvers/samplers.jl")


include("billiards/boundarypoints.jl")

#basis
#include("basis/Basis.jl")
#@reexport using .Basis
include("basis/planewaves/realplanewaves.jl")

include("basis/fourierbessel/corneradapted.jl")

include("basis/fundamental/fundamental.jl")
#include("basis/fundamental/fundamentalbessels.jl")


#billiards
#include("billiards/Billiards.jl")
#@reexport using .Billiards

include("billiards/geometry/geometry.jl")

include("billiards/stadium.jl")
include("billiards/lemon.jl")
include("billiards/sinai.jl")
include("billiards/triangle.jl")

#include("limacon.jl")
#include("rectangle.jl")

#convenience functions may be moved somewhere else
#export make_stadium_and_basis, make_triangle_and_basis 


include("solvers/decompositions.jl")
include("solvers/matrixconstructors.jl")


include("solvers/acceleratedmethods/acceleratedmethods.jl")
include("solvers/sweepmethods/sweepmethods.jl")


#spectra
#include("spectra/Spectra.jl")
#@reexport using .Spectra

include("spectra/spectralutils.jl")
include("spectra/unfolding.jl")

#states
#include("states/States.jl")
#@reexport using .States
include("states/eigenstates.jl")
include("states/basisstates.jl")
include("states/randomstates.jl")


include("states/wavefunctions.jl")
include("states/boundaryfunctions.jl")
include("states/husimifunctions.jl")



#plotting functions in Makie
#include("plotting/Plotting.jl")
#@reexport using .Plotting
include("plotting/plottingmakie.jl")

include("plotting/testplotting.jl")


include("utils/benchmarkutils.jl")


include("plotting/benchmarkplotting.jl")


#convenience functions
function make_stadium_and_basis(half_width;radius=1.0,x0=zero(half_width),y0=zero(half_width), rot_angle=zero(half_width))
    billiard = Stadium(half_width; radius=radius,x0=x0,y0=y0)
    basis = CornerAdaptedFourierBessel(1, pi/2.0, SVector(x0,y0),rot_angle) 
    return billiard, basis 
end

function make_lemon_and_basis(half_separation;full_domain=false,radius=1.0,x0=zero(half_separation),y0=zero(half_separation), rot_angle=zero(half_separation))
    billiard = Lemon(half_separation; radius=radius,x0=x0,y0=y0)
    basis = CornerAdaptedFourierBessel(1, pi/2.0, SVector(x0,y0),rot_angle) 
    return billiard, basis 
end

function make_triangle_and_basis(gamma,chi; edge_i=1)
    cor = Triangle(gamma,chi).corners
    x0,y0 = cor[mod1(edge_i+2,3)]
    re = [:Virtua, :Virtual, :Virtual]
    re[edge_i] = :Real 
    tr = Triangle(gamma,chi; curve_types = re, x0 = x0, y0 =y0)
    basis = CornerAdaptedFourierBessel(1, adapt_basis(tr,edge_i+2)...) 
    return tr, basis 
end
