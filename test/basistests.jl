include("../src/QuantumBilliards.jl")
#using Revise
using .QuantumBilliards
using Test

@testset "Basis tests: CornerAdaptedFourierBessel" begin
    basis = CornerAdaptedFourierBessel(10, pi/2,SVector(0.0,0.0),0.0)
    nx = 10
    ny = 10
    x_grid = range(0.0, 1.0, nx)
    y_grid = range(0.0, 1.0, ny)
    pts = [SVector(x,y) for y in y_grid for x in x_grid]
    @test typeof(basis_fun(basis,1,10.0,pts)) <: AbstractArray
    @test typeof(basis_fun(basis,1:10,10.0,pts)) <: AbstractArray
    @test typeof(gradient(basis,1,10.0,pts)) <: Tuple
    @test typeof(gradient(basis,1:10,10.0,pts)) <: Tuple
    @test typeof(basis_and_gradient(basis,1,10.0,pts)) <: Tuple
    @test typeof(basis_and_gradient(basis,1:10,10.0,pts)) <: Tuple
    @test typeof(dk_fun(basis,1,10.0,pts)) <: AbstractArray
    @test typeof(dk_fun(basis,1:10,10.0,pts)) <: AbstractArray
end;