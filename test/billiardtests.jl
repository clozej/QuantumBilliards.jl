include("../src/QuantumBilliards.jl")
#using Revise
using .QuantumBilliards
using Test

@testset "Billiard tests: Stadium" begin
    width = 0.5
    billiard = Stadium(width)
    t = collect(range(0.0,1.0,10))
    for crv in billiard.boundary
        @test typeof(curve(crv,t)) <: Vector
        @test typeof(normal_vec(crv,t)) <: Vector
        @test typeof(arc_length(crv,t)) <: Vector
    end
end;

@testset "Billiard tests: Triangle" begin
    gamma = 2/3*pi #sqrt(2)/2 * pi
    chi  = 2.0
    billiard = Triangle(gamma,chi)
    t = collect(range(0.0,1.0,10))
    for crv in billiard.boundary
        @test typeof(curve(crv,t)) <: Vector
        @test typeof(normal_vec(crv,t)) <: Vector
        @test typeof(arc_length(crv,t)) <: Vector
    end
end;
