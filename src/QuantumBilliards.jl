module QuantumBilliards
#abstract types
include("abstracttypes.jl")

#billiards
include("billiards/billiard.jl")
include("billiards/curves.jl")
include("billiards/geometry.jl")

include("billiards/limacon.jl")
include("billiards/rectangle.jl")
include("billiards/triangle.jl")
include("billiards/stadium.jl")

#plotting
include("plotting/plottingmakie.jl")
include("plotting/matrixplotting.jl")

#solvers
include("solvers/acceleratedmethods.jl")
include("solvers/boundaryintegralmethod.jl")
include("solvers/decompositionmethod.jl")
include("solvers/decompositions.jl")
include("solvers/matrixconstructors.jl")
include("solvers/particularsolutionsmethod.jl")
include("solvers/samplers.jl")
include("solvers/scalingmethod.jl")
include("solvers/sweepmethods.jl")

#spectra
include("spectra/spectralutils.jl")
include("spectra/unfolding.jl")

#states
include("states/basisstates.jl")
include("states/boundaryfunctions.jl")
include("states/eigenstates.jl")
include("states/husimifunctions.jl")
include("states/randomstates.jl")
include("states/wavefunctions.jl")

#utils
include("utils/angleutils.jl")
include("utils/benchmarkutils.jl")
end
