include("../src/QuantumBilliards.jl")
using Revise
using .QuantumBilliards 
using StaticArrays

basis = CornerAdaptedFourierBessel(1,pi/2,SVector(0.0,0.0), 1/3)
