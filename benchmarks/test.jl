
using MKL
include("../src/QuantumBilliards.jl")
#using Revise
using .QuantumBilliards
#using Revise 
#using GLMakie
using CairoMakie
#using Latexify
using StaticArrays

crv = CircleSegment(2.0, pi/2.0, 0.0, 0.0,0.0)
crv2 = LineSegment(SVector(0.0,0.0),SVector(1.0,1.0))

ts = collect(range(0.0,1.0,100))

kappa = curvature(crv, ts)