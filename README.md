# QuantumBilliards

[![Build Status](https://github.com/clozej/QuantumBilliards.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/clozej/QuantumBilliards.jl/actions/workflows/CI.yml?query=branch%3Amain)

A package for the computation, analysis and visualisation of spectra and states of quantum billiards. The package is a work in progres.

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
- add arbitrary precision features

Documentation:
- all of it

Fixes and improvements
- improve wavefunction computation for large memory
- add type conversions
- separate plotting to separate package
- make unit tests
- fix bugs in symmetrized plots
- add more samplers
- fix integration weights for nonlinear boundary functions 


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
