# QuantumBilliards

[![Build Status](https://github.com/clozej/QuantumBilliards.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/clozej/QuantumBilliards.jl/actions/workflows/CI.yml?query=branch%3Amain)

A package for the computation, analysis and visualisation of spectra and states of quantum billiards. 

# List of computational methods
The folowing methods have been fully or partialy implemented:
- Scaling method 
- Decomposition method
- Particular solutions method
# To do list
Solvers:
- Boundary integral method
- Vebles accelerated method
- Gague Freedom method
Basis sets:
- Fourier Bessel functions (symmetrized)
- Fundamnetal basis Hankel functions
Billiards:
- Mushroom
- Squash
- Cut circle
- Limacon
- Squrcle
- Circle
- Rectangle
- Barrier
- Isospectral polygons
- Rough billiards
- Annular
- Track billiards
New features:
- Wigner functions
- time evolution
- classical billiards
- new symmetry clases (rotational symmetry)
- neutrino billiards
- optical cavities
- scattering solutions
Documentation:
- all of it


# Dependencies
- Bessels
- CircularArrays
- CoordinateTransformations
- Distributions
- FFTW 
- FastGaussQuadrature
- ForwardDiff
- IntervalArithmetic
- LinearAlgebra
- Makie
- Optim
- Random
- Rotations
- StaticArrays
- TimerOutputs